import os
import h5py
import allel
import numpy as np
import dask.array as da
import pandas as pd
import zarr

# compute discordant reads based on the first observed allele.
def compute_discords(genotype, allele_depth):
    
    index = genotype[:, 0]
    x = allele_depth.compute()
    n_disc = x[np.arange(0, x.shape[0], dtype=int), index]
    ad = x.sum(axis=1)
    return ad - n_disc, ad

chrom = snakemake.wildcards.chrom

phase2_callset_pass = zarr.open_group(snakemake.input.phase2_callset, mode="r")
pass_pos = allel.SortedIndex(phase2_callset_pass[chrom]["variants/POS"])

x_callset = h5py.File(snakemake.input.cross_callset, mode="r")
xdf = pd.read_table(
    snakemake.input.metadata,
    index_col=0)

x_pos = allel.SortedIndex(x_callset[chrom]['variants/POS'])
x_gt = allel.GenotypeDaskArray(x_callset[chrom]['calldata/genotype'])

x_ad = x_callset[chrom]['calldata/AD']
x_ad = da.from_array(x_ad, chunks=x_ad.chunks)

call_class = ("HOMREF", "HET", "HOMALT")
columns = pd.MultiIndex.from_product(
    (("ALL", "PASS"),
     call_class,
     ("SUCCESS", "N")))

# take sample names from the hdf5 file
sample_list = x_callset[chrom]["samples"][:].astype("U8").tolist()

# Drop samples that were not included in sequencing.
xdf = xdf.set_index("ox_code").reindex(sample_list).reset_index()
xdf = xdf.loc[xdf.cross.notna()]
xids = xdf.cross.unique()
assert len(xids) == 10, xids

# HOM REF, ALT, HET
# ALL / PASS
# sites success / sites called
df = pd.DataFrame(index=xids, columns=columns, dtype=int)

for cross_id in xids:

    print("processing", cross_id, "...")

    par_ix = xdf.query('cross == @cross_id').query('role == "parent"').index.values
    pro_ix = xdf.query('cross == @cross_id').query('role == "progeny"').index.values
    if par_ix.size != 2:
        print("Must be two parents: {0} found".format(par_ix.size))
        continue

    # grab genotypes of cross and AD
    pr_gt = x_gt.take(par_ix, axis=1)
    pg_gt = x_gt.take(pro_ix, axis=1)
    pr_ad = da.take(x_ad, par_ix, axis=1)

    # count hom refs and alts of parents
    hom_alt_sum = pr_gt.is_hom_alt().sum(axis=1).compute()
    hom_ref_sum = pr_gt.is_hom_ref().sum(axis=1).compute()

    # identify discordance
    mat_discords, mat_cov = compute_discords(pr_gt[:, 0], pr_ad[:, 0])
    pat_discords, pat_cov = compute_discords(pr_gt[:, 1], pr_ad[:, 1])
    
    # identify high coverage sites
    loc_high_cov = (mat_cov >= 30) & (pat_cov >= 30)
    loc_no_discords = (mat_discords + pat_discords) == 0
    loc_sufficient_calls = pg_gt.count_called(axis=1).compute() >= 10
    loc_all_filters = (loc_high_cov & loc_sufficient_calls) & loc_no_discords
    
    # where do we expect to find het/hom-alt/hom-ref
    loc_het_expected = (hom_alt_sum == 1) & (hom_ref_sum == 1)
    loc_homref_expected = hom_ref_sum == 2
    loc_homalt_expected = hom_alt_sum == 2

    expected_genotypes = (
        loc_homref_expected, loc_het_expected, loc_homalt_expected)

    count_funcs = (
        lambda x: x.count_hom_ref(axis=1), 
        lambda x: x.count_het(allele=0, axis=1), 
        lambda x : x.count_hom_alt(axis=1))

    for expected_loc, func, label in zip(expected_genotypes, count_funcs, call_class):

        loc_review = expected_loc & loc_high_cov & loc_sufficient_calls & loc_no_discords
        
        gt_exp = pg_gt.compress(loc_review, axis=0)
        pos_exp = x_pos.compress(loc_review, axis=0)

        obs = func(gt_exp).compute()
        n_calls = gt_exp.count_called(axis=1).compute()

        # handle_pass
        loc_pass, _ = pos_exp.locate_intersection(pass_pos)

        obs_pass = np.compress(loc_pass, obs, axis=0)
        n_calls_pass = np.compress(loc_pass, n_calls, axis=0)

        # must have at least 10 calls to count....
        df.loc[cross_id]["ALL"][label] = np.sum(obs == n_calls), n_calls.shape[0]
        df.loc[cross_id]["PASS"][label] = np.sum(obs_pass == n_calls_pass), n_calls_pass.shape[0]

df.to_csv(snakemake.output.csv)

