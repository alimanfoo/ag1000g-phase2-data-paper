import sys
sys.path.insert(0, '/home/eric/vector-ops/agam-report-base/src/python')
from util import *
import zcache
import veff
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ag1k_dir = '/kwiat/vector/ag1000g/release'
from ag1k import phase2_ar1
phase2_ar1.init(os.path.join(ag1k_dir, 'phase2.AR1'))

# Get the haplotype meta data
hap_meta_auto = pd.read_csv('/kwiat/vector/ag1000g/release/phase2.AR1/haplotypes/haplotypes.autosomes.meta.txt', sep = '\t', index_col = 1)
hap_meta_X = pd.read_csv('/kwiat/vector/ag1000g/release/phase2.AR1/haplotypes/haplotypes.X.meta.txt', sep = '\t', index_col = 1)
hap_meta = {'2L': hap_meta_auto, 
            '2R': hap_meta_auto, 
            '3L': hap_meta_auto, 
            '3R': hap_meta_auto, 
            'X' : hap_meta_X}

# Here are the chromosomes
chroms = ['2L', '2R', '3L', '3R', 'X']
# Here is where we will save the boxplot figure
boxplot_filename = 'refdiff_phase2_boxplot.eps'
# This is where we will store the output 
dxy_by_pop = dict()
hamming = pd.DataFrame(columns = chroms)
NEED TO MAKE SURE THE SAME SAMPLES ARE OBTAINED FOR ALL THE CHROMS
accessible_sites = pd.Series(index = chroms)

# Get the haplotypes
haplotypes = phase2_ar1.callset_phased
# Get the accessibility
accessible = phase2_ar1.accessibility

for chrom in chroms:
	# First, for each haplotype, calculate the hamming distance between the haplotype and the reference
	# genome, which is just the sum of all the alt alleles. 
	ham = np.count_nonzero(haplotypes[chrom]['calldata/genotype'], 0)
	# If running on a computer with low available RAM, it will struggle with the above line of code. The
	# following three lines can be used instead, but runs much, much more slowly.
#	ham = np.empty((0,2))
#	for i in range(haplotypes[chrom]['calldata/genotype'].shape[1]):
#		ham = np.concatenate((ham, np.expand_dims(np.count_nonzero(haplotypes[chrom]['calldata/genotype'][:,i,:], 0), 0)), 0)
	# Then divide by the number of accessible positions to get the dxy. 
	acc = np.count_nonzero(accessible[chrom]['is_accessible'])
	# Add these distance and accessible site values to the storing dictionaries
	hamming[chrom] = ham.flatten()
	accessible_sites[chrom] = acc
	# Calculate dxy for this chromosome by dividing the hamming distance by the number of accessible sites. 
	dxy = ham.flatten()/acc
	# Now, for each population, calculate the mean and standard error of the dxy
	dxy_by_pop[chrom] = pd.DataFrame([], columns = ['mean', 'std'])
	for p in phase2_ar1.pop_ids:
		these_dxy = dxy[hap_meta[chrom]['population'] == p]
		dxy_by_pop[chrom].loc[p, :] = [np.mean(these_dxy), np.std(these_dxy)]

# Calculate the genome-wide dxy for each sample
all_the_ham = np.sum(hamming, 1)
all_the_sites = np.sum(accessible_sites)
all_the_dxy = pd.DataFrame(all_the_ham / all_the_sites, columns = ['dxy'])
		
# Plot the boxplot	
def dxy_boxplot(data, group, fig, ax, pal):
	grouped_dxy = data.groupby(list(group))
	bx = ax.boxplot([x[1]['dxy'] for x in list(grouped_dxy)], notch = True, bootstrap = 1000, whis = [5, 95], showfliers = False, 
					patch_artist= True, medianprops = dict(linestyle='-', color='k'),
					whiskerprops = dict(linestyle = '-', linewidth = 1.5, color = 'k'),
					boxprops = dict(linestyle = '-', linewidth = 1.5, color = 'k'))
	ax.set_xticklabels(list(grouped_dxy.groups.keys()), rotation = 45, ha = 'right', fontsize = 8) 
	ax.set_xlabel('Population')
	ax.set_ylabel('Dxy')
	ax.grid(axis = 'y')
	for patch, color in zip(bx['boxes'], palette):
		patch.set_facecolor(color)
	fig.tight_layout()

palette = sns.color_palette('husl',n_colors = len(phase2_ar1.pop_ids))
fig, ax = plt.subplots(figsize=(7, 4))
dxy_boxplot(all_the_dxy, hap_meta[chrom]['population'], fig, ax, palette)
plt.savefig(output_filename, format = 'eps', dpi = 1000)

# Now we calculate the dxy in windows across the different chromosomes and plot them against their
# genomic position

# Let's just do X first
# Calculate dxy in windows across the chromosome. 
win_size = 10000
chrom_size = len(phase2_ar1.genome_agamp3[chrom])
# We use the lower bound of the fraction for the number of windows, which means that the final, incomplete
# window of the chromosome won't get used.
win_num = chrom_size // win_size
dxy_by_window = dict(gam = dict(), col = dict())
for chrom in chroms:
	dxy_by_window['gam'][chrom] = [0]*win_num
	dxy_by_window['col'][chrom] = [0]*win_num
	# This will be quicker if we actually load the calldata fully into memory first, but that
	# won't work on the laptop (or even possibly the desktop)
	for i in range(win_num):
		these_positions = SOMETHING WITH POS AND SLICE
		ham = np.count_nonzero(haplotypes[chrom]['calldata/genotype'][i:(i+win_size),:,:], 0).flatten()
		acc = np.count_nonzero(accessible[chrom]['is_accessible'][i:(i+win_size)])
		dxy = ham/acc
		# Now get the mean dxy by species. 
		dxy_by_window['gam'][chrom][i] = 
		dxy_by_window['col'][chrom][i] = 


