import scipy
import numpy as np
import pandas as pd
import allel
import h5py
from collections import OrderedDict
import matplotlib.pyplot as plt
import sys
import os
from re import *
import seaborn as sns

#ag1k_dir = '/kwiat/vector/ag1000g/release'
ag1k_dir = '/home/elucas/data'
from ag1k import phase2_ar1
phase2_ar1.init(os.path.join(ag1k_dir, 'phase2.AR1'))
chroms = list(phase2_ar1.callset_phased.keys())

# Give species their full names in the metadata
samples_df = phase2_ar1.df_samples.copy()
samples_df['m_s'] = samples_df['m_s'].fillna('M/S')
samples_df['species'] = [re.sub('M/S', 'Unknown', str(x)) for x in samples_df['m_s']]
samples_df['species'] = [re.sub('M', 'An. coluzzii', str(x)) for x in samples_df['species']]
samples_df['species'] = [re.sub('S', 'An. gambiae', str(x)) for x in samples_df['species']]

# Create columns concatenating two population and sex (or species and sex)
samples_df['pop_sex'] = samples_df['population'] + '_' + samples_df['sex']
samples_df['sp_sex'] = samples_df['species'] + '_' + samples_df['sex']

subpops = samples_df.groupby('pop_sex').indices
subpops_ix = {k: list(v) for k, v in subpops.items()}
species = samples_df.groupby('sp_sex').indices
species_ix = {k: list(v) for k, v in species.items()}

# Set the populations of interest
pop = ['GW', 'GM', 'KE']

gff_fn = ("/home/eric/Liverpool/AR3/geneset/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.2.sorted.gff3.gz")
fasta_fn = ("/home/eric/Liverpool/AR3/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa")
accessibility_fn = "/home/eric/Liverpool/phase2.AR1/accessibility/accessibility.h5"
callset_fn = '/media/eric/Ucalegon/variation/phase2/main/hdf5/pass/ag1000g.phase2.ar1.pass.' + chrom + '.h5'

haplotypes_file = '/media/eric/Ucalegon/haplotypes_phase2/main/hdf5/ag1000g.phase2.ar1.haplotypes.' + chrom + '.h5'
haplotypes_fstem = '/media/eric/Ucalegon/haplotypes_phase2/main/hdf5/ag1000g.phase2.ar1.haplotypes.{chrom}.h5'
meta_file = '/home/eric/Liverpool/phase2.AR1/samples/samples.meta.txt'
hap_X_metadata = '/media/eric/Ucalegon/haplotypes_phase2/haplotypes.X.meta_manual.txt'

window_size = 100000

chrom = '2R'

GET THE CHROM RIGHT HERE
ag1k_hap_data = selection.data.HaplotypeData(hap_X_metadata, fasta_fn, haplotypes_fstem, gff_fn, accessibility_fn, gene_labels, populations, pop_labels, subsets, pop_colors)

# Get sample metadata
meta = pd.read_csv(meta_file, sep = '\t')
hap_X_meta = pd.read_csv(hap_X_metadata, sep = '\t', index_col = 0)

# Get the metadata for the population of interest
meta_pop = meta.loc[meta['population'].isin(pop)]
hap_X_meta_pop = hap_X_meta.loc[hap_X_meta['population'].isin(pop)]
# Get the ox codes for the population of interest
focal_samples = meta_pop['ox_code']
focal_females = meta_pop['ox_code'].iloc[np.where(meta_pop['sex'] == 'F')[0]]

# Load the haplotypes file
all_haplotypes = h5py.File(haplotypes_file, 'r')

# Get the indices of the samples of interest
sample_list = pd.Series([x.decode('utf-8') for x in all_haplotypes[chrom]['samples'].value])
indices_list = np.where(sample_list.isin(focal_samples))[0]

accessibility = h5py.File(accessibility_fn, "r")
pos = allel.SortedIndex(all_haplotypes[chrom]['variants']['POS'])

downstream_hap_matrix = allel.GenotypeArray(all_haplotypes[chrom]['calldata']['genotype'][range(downstream_start_index, downstream_end_index)][:,indices_list,:]).to_haplotypes()
dist = allel.stats.pairwise_distance(downstream_hap_matrix, metric = 'hamming')
is_accessible = accessibility[chrom]['is_accessible'][pos[downstream_start_index]:pos[downstream_end_index]]
n_bases = np.count_nonzero(is_accessible)
dist_dxy = dist * downstream_hap_matrix.n_variants / n_bases


