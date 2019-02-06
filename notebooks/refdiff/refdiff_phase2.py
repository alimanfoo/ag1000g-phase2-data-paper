import sys
sys.path.insert(0, '../../agam-report-base/src/python')
from util import *
import zcache
import veff
import os
import numpy as np
import pandas as pd
import allel
import re

#ag1k_dir = '/kwiat/vector/ag1000g/release'
# On Liverpool server:
ag1k_dir = '/home/elucas/data'
from ag1k import phase2_ar1
phase2_ar1.init(os.path.join(ag1k_dir, 'phase2.AR1'))
chroms = list(phase2_ar1.callset_phased.keys())
colours = phase2_ar1.pop_colors
colours['An. gambiae'] = 'red'
colours['An. coluzzii'] = 'blue'

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

def compute_divergence(allele_freqs):
	sum_alt = allele_freqs.sum(axis=0)
	return (sum_alt[1:].sum())

window_size = 100000
dxy_by_window = dict()

for chrom in chroms:
	print('\nChromosome ' + chrom)
	gt = allel.GenotypeDaskArray(phase2_ar1.callset_pass[chrom]['calldata/genotype'])
	ac = gt.count_alleles_subpops(subpops_ix)
	ac_species = gt.count_alleles_subpops(species_ix)

	pos = allel.SortedIndex(phase2_ar1.callset_pass[chrom]['variants/POS'])

	eqa = allel.equally_accessible_windows(phase2_ar1.accessibility[chrom]['is_accessible'], window_size)

	# Use the middle of the window as the index
	window_middle = np.sum(eqa, axis=1)/2
	dxy_by_window[chrom] = pd.DataFrame(index=window_middle.astype(int), columns=list(subpops.keys()) + list(species.keys()))
	dxy_by_window[chrom] = pd.DataFrame(index=window_middle.astype(int), columns= species.keys())

	for pop in subpops.keys():
		print('processing', pop)

		# Faster if we drop non variant loci first, and load into mem
		loc = ac[pop].is_variant().compute()
		pop_ac = ac[pop].compress(loc, axis=0).compute()
		pop_pos = pos.compress(loc, axis=0)
		print('computing divergence...', loc.sum())

		vals, windows, counts = allel.windowed_statistic(
			pop_pos, pop_ac.to_frequencies(), compute_divergence, windows=eqa)
		dxy_by_window[chrom][pop] = vals / window_size

	for pop in species.keys():
		print('processing', pop)

		# Faster if we drop non variant loci first, and load into mem
		loc = ac_species[pop].is_variant().compute()
		pop_ac = ac_species[pop].compress(loc, axis=0).compute()
		pop_pos = pos.compress(loc, axis=0)
		print('computing divergence...', loc.sum())

		vals, windows, counts = allel.windowed_statistic(
			pop_pos, pop_ac.to_frequencies(), compute_divergence, windows=eqa)
		dxy_by_window[chrom][pop] = vals / window_size
	# And save that table to file
	dxy_by_window[chrom].to_csv('refdiff_phase2_' + chrom + '_table.csv', sep = '\t')

