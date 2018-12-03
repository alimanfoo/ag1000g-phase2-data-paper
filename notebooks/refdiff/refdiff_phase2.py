import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
import seaborn as sns
sns.set_context('paper')
sns.set_style('white')
sns.set_style('ticks')
rcParams = plt.rcParams
rcParams['font.family'] = 'sans-serif'
_font_size = 7
rcParams['font.size'] = _font_size
rcParams['figure.titlesize'] = _font_size
rcParams['axes.titlesize'] = _font_size
rcParams['axes.labelsize'] = _font_size
rcParams['xtick.labelsize'] = _font_size
rcParams['ytick.labelsize'] = _font_size
rcParams['legend.fontsize'] = _font_size
_line_width = .5
rcParams['axes.linewidth'] = _line_width
rcParams['lines.linewidth'] = _line_width
rcParams['patch.linewidth'] = _line_width
rcParams['ytick.direction'] = 'out'
rcParams['xtick.direction'] = 'out'
rcParams['savefig.jpeg_quality'] = 100
rcParams['figure.dpi'] = 120
rcParams['lines.markeredgewidth'] = _line_width
rcParams['figure.figsize'] = (4.85, 3)
#rcParams['text.usetex'] = True
#rcParams['text.latex.unicode'] = True
rcParams['text.usetex'] = False
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.tt'] = 'Courier New'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.sf'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['svg.fonttype'] = 'none'
#%config InlineBackend.figure_formats = {'retina', 'png'}

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

# Set up the chromatin data for the genome
_data_chromatin = b"""CHX     chro    X       20009764        24393108
CH2R    chro    2R      58984778        61545105
CH2L    chro    2L      1       2431617
PEU2L   chro    2L      2487770 5042389
IH2L    chro    2L      5078962 5788875
IH3R    chro    3R      38988757        41860198
CH3R    chro    3R      52161877        53200684
CH3L    chro    3L      1       1815119
PEU3L   chro    3L      1896830 4235209
IH3L    chro    3L      4264713 5031692
"""
tbl_chromatin = (
    etl
    .fromtext(etl.MemorySource(_data_chromatin))
    .split('lines', '\s+', ['name', 'type', 'chrom', 'start', 'stop'])
    .convert(('start', 'stop'), int)
    .cutout('type')
)

def fig_linear_genome(plotf, genome, chromosomes=None, fig=None, 
					  bottom=0, height=1, width_factor=1.08, chrom_pad=0.035, 
					  clip_patch_kwargs=None, **kwargs):
	if chromosomes is None:
		chromosomes = ['2R', '2L', '3R', '3L', 'X']
	genome_size = sum(len(genome[chrom]) for chrom in chromosomes)

	from matplotlib.path import Path

	if fig is None:
		fig = plt.figure(figsize=(8, 1))

	left = 0

	if clip_patch_kwargs is None:
		clip_patch_kwargs = dict()
	clip_patch_kwargs.setdefault('edgecolor', 'k')
	clip_patch_kwargs.setdefault('facecolor', 'none')
	clip_patch_kwargs.setdefault('lw', 1)

	axs = dict()
	for chrom in chromosomes:

		# calculate width needed for this chrom
		width = len(genome[chrom]) / (genome_size * width_factor)

		# create axes
		ax = fig.add_axes([left, bottom, width, height])
		#ax.set_axis_bgcolor((1, 1, 1, 0));
		axs[chrom] = ax

		# construct clip path
		if chrom in {'2R', '3R'}:
			verts = [(0.01, 0.02), (0.9, 0.02), (1.01, 0.3), (1.01, 0.7), (0.9, .98), (0.01, .98), (0.01, 0.02)]
			codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
		elif chrom == "X":
			verts = [(0.01, 0.02), (0.9, 0.02), (0.99, 0.3), (0.99, 0.7), (0.9, .98), (0.01, .98), (0.01, 0.02)]
			codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
		else:
			verts = [(0.1, 0.02), (.99, 0.02), (.99, .98), (.1, .98), (-0.01, .7), (-0.01, .3), (0.1, 0.02)]
			codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
		path = Path(verts, codes)
		clip_patch = mpl.patches.PathPatch(path, transform=ax.transAxes, **clip_patch_kwargs)

		# do the plotting
		plotf(chrom=chrom, ax=ax, clip_patch=clip_patch, genome=genome, **kwargs)

		# increment left coordinate
		left += len(genome[chrom]) / (genome_size * width_factor)
		if chrom in {'2L', '3L'}:
			left += chrom_pad

	return axs


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

# Identify samples from the wild populations
populations = phase2_ar1.pop_ids
hap_wild_auto = np.where(np.isin(hap_meta_auto['population'], populations))[0]
hap_wild_X = np.where(np.isin(hap_meta_X['population'], populations))[0]
hap_wild = {'2L': hap_wild_auto, 
			'2R': hap_wild_auto, 
			'3L': hap_wild_auto, 
			'3R': hap_wild_auto, 
			'X' : hap_wild_X}


# Here are the chromosomes
chroms = ['2L', '2R', '3L', '3R', 'X']
autosomes = chroms[:4]
# Here is where we will save the figures
boxplot_filename = 'refdiff_phase2_boxplot.jpg'
chromplot_filename = 'refdiff_phase2_chromplot.jpg'
# This is where we will store the output 
hamming = pd.DataFrame(columns = autosomes)
accessible_sites = pd.Series(index = autosomes)
dxy_by_pop = dict()
dxy_by_window = dict(gam = dict(), col = dict())

# Get the haplotypes
haplotypes = phase2_ar1.callset_phased
# Get the accessibility
accessible = phase2_ar1.accessibility

win_size = 500000
for chrom in chroms:
	# I think this will be quicker if we actually load the calldata fully into memory first.
	these_haplotypes = haplotypes[chrom]['calldata/genotype'].value
	this_accessibility = accessible[chrom]['is_accessible'].value
	# First, calculate the dxy across each chromosome
	# For each haplotype, calculate the hamming distance between the haplotype and the reference
	# genome, which is just the sum of all the alt alleles. 
	ham = np.count_nonzero(these_haplotypes, 0)
	# Then count the number of accessible positions to get the dxy. 
	acc = np.count_nonzero(this_accessibility)
	# Calculate dxy for this chromosome by dividing the hamming distance by the number of accessible sites. 
	dxy = ham.flatten()/acc
	# Now, for each population, calculate the mean and standard error of the dxy
	dxy_by_pop[chrom] = pd.DataFrame([], columns = ['mean', 'std'])
	for p in populations:
		these_dxy = dxy[hap_meta[chrom]['population'] == p]
		dxy_by_pop[chrom].loc[p, :] = [np.mean(these_dxy), np.std(these_dxy)]
	# If this is an autosome, add the hamming and accessible site values to the storing objects, keeping only 
	# the hamming values for the wild populations
	if chrom in autosomes:
		hamming[chrom] = ham.flatten()[hap_wild_auto]
		accessible_sites[chrom] = acc

	# Second, calculate dxy by window
	pos = allel.SortedIndex(haplotypes[chrom]['variants']['POS'])
	chrom_size = len(phase2_ar1.genome_agamp3[chrom])
	# We use the lower bound of the fraction for the number of windows, which means that the final, incomplete
	# window of the chromosome won't get used.
	win_num = chrom_size // win_size
	# Set up the lists that will store the output
	dxy_by_window['gam'][chrom] = [0]*win_num
	dxy_by_window['col'][chrom] = [0]*win_num
	for i in range(win_num):
		# Find the indices of the variants present in this section of the genome. If there aren't any, we 
		# return 0 (so if there are no accessible sites, we default to a dxy of 0).
		try:
			loc = pos.locate_range(i*win_size + 1, (i+1)*win_size)
		except KeyError:
			dxy_by_window['gam'][chrom][i] = (0, i*win_size + (win_size/2))
			dxy_by_window['col'][chrom][i] = (0, i*win_size + (win_size/2))
		else:
			acc = np.count_nonzero(this_accessibility[(i*win_size):((i+1)*win_size)])
			# Calculate the dxy
			ham = np.count_nonzero(these_haplotypes[loc,:,:], 0).flatten()
			dxy = ham/acc
			# Now get the mean dxy by species, and record the mean position of each window
			dxy_by_window['gam'][chrom][i] = (np.mean(dxy[hap_meta[chrom]['m_s'] == 'S']), i*win_size + (win_size/2))
			dxy_by_window['col'][chrom][i] = (np.mean(dxy[hap_meta[chrom]['m_s'] == 'M']), i*win_size + (win_size/2))


# delete the haplotype object to clear memory
del these_haplotypes

# Calculate the genome-wide dxy for each sample
all_the_ham = np.sum(hamming, 1)
all_the_sites = np.sum(accessible_sites)
all_the_dxy = pd.DataFrame(all_the_ham / all_the_sites, columns = ['dxy'])
all_the_dxy['population'] = list(hap_meta_auto['population'][hap_wild_auto])
all_the_dxy.to_csv('all_the_dxy.csv', sep = '\t')
		
# Plot the boxplot	
def dxy_boxplot(data, group, fig, ax, pal, fn=None, save_dpi=200):
	grouped_dxy = data.groupby(list(group))
	bx = ax.boxplot([x[1] for x in list(grouped_dxy)], notch = True, bootstrap = 1000, whis = [5, 95], showfliers = False, 
					patch_artist= True, medianprops = dict(linestyle='-', color='k'),
					whiskerprops = dict(linestyle = '-', linewidth = 1.5, color = 'k'),
					boxprops = dict(linestyle = '-', linewidth = 1.5, color = 'k'))
	ax.set_xticklabels(list(grouped_dxy.groups.keys()), rotation = 45, ha = 'right', fontsize = 8) 
	ax.set_xlabel('Population')
	ax.set_ylabel('Autosomal Dxy')
	ax.grid(axis = 'y')
	for patch, color in zip(bx['boxes'], palette):
		patch.set_facecolor(color)
	fig.tight_layout()
	if fn:
		if re.search('\.eps', fn):
			fig.savefig(fn, format='eps', dpi=save_dpi, bbox_inches='tight')
		else:
			fig.savefig(fn, jpeg_quality=100, dpi=save_dpi, bbox_inches='tight')

palette = sns.color_palette('husl',n_colors = len(populations))
fig, ax = plt.subplots(figsize=(7, 4))
dxy_boxplot(all_the_dxy['dxy'], all_the_dxy['population'], fig, ax, palette, boxplot_filename)


# Now plot the dxy along the chromsomomes
def plot_dxy(chrom, ax, dxy, ymax=0.02, offset=5, legend_frameon=False, **kwargs):
	if chrom == '2R':
		sns.despine(ax=ax, bottom=True, offset=offset)
		ax.set_yticks([0, .01, .02], minor=False)
		ax.set_yticks(np.arange(0, ymax, .002), minor=True)
		ax.set_yticklabels([0, .01, .05])
		ax.set_ylabel('Dxy', ha='left', va='top', rotation=0)
		ax.yaxis.set_label_coords(-0.25, 1.3)
	else:
		sns.despine(ax=ax, left=True, bottom=True, offset=offset)
		ax.set_yticks([])
		ax.set_title(chrom)

	# Set the values of the x axis as being the mid-point of each window
	ygam,xgam = zip(*dxy['gam'][chrom])
	ycol,xcol = zip(*dxy['col'][chrom])
	
	h1 = ax.plot(xgam, ygam, color='red', label = 'An. gambiae')
	h2 = ax.plot(xcol, ycol, color='blue', label = 'An. coluzzii')

	if chrom == 'X':
		handles = []
		for h in [h1, h2][::-1]:
			handles.append(plt.plot([0, 1], [0, 0], color = h[0].get_c(), lw=2, label=h[0].get_label())[0])
		ax.legend(handles=handles, bbox_to_anchor=(1, 1.5), loc='upper left', frameon=legend_frameon, fancybox=False, 
				  title='Species', labelspacing=.1);

	ax.set_ylim(0, ymax)
	ax.set_xlim(0, max(xgam))
	ax.set_xticks([])
	ax.set_title(chrom)

def plot_heterochromatin(chrom, ax, clip_patch, tbl_chr, genome, colors=None, legend_frameon=False, **kwargs):
    if colors is None:
        colors = {'PE': 'w', 'IH': 'k', 'CH': 'k'}

    recs = tbl_chr.eq('chrom', chrom).records()
    color = [colors[rec.name[:2]] for rec in recs]
    xranges = [(rec.start, rec.stop-rec.start) for rec in recs]
    a = ax.broken_barh(xranges, yrange=(0, 1), color=color)
    ax.set_xlim(0, len(genome[chrom]))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()
    
    a.set_clip_path(clip_patch) 
    ax.add_patch(clip_patch)

    if chrom == 'X':
        a = plt.Rectangle((0, 0.2), 1, 1, ec='k', fc='w', label='euchromatin')
        b = plt.Rectangle((0, 0.2), 1, 1, ec='k', fc='k', label='heterochromatin')
        handles = [a, b]
        ax.legend(handles=handles, bbox_to_anchor=(1, 2), loc='center left', frameon=legend_frameon, 
                  fancybox=True, labelspacing=.1, title='Chromatin state');


def fig_assemble(fw=4.3, fn=None, dpi=150, save_dpi=200, legend_frameon=False, fh=0.8):
	figsize = fw, fh
	fig = plt.figure(figsize=(fw, fh), dpi=dpi)
	fig_linear_genome(plot_dxy, phase2_ar1.genome_agamp3, chromomsomes = chroms, fig=fig, bottom=0.1, height=.7, legend_frameon=legend_frameon, dxy = dxy_by_window)
	fig_linear_genome(plot_heterochromatin, phase2_ar1.genome_agamp3, fig=fig, bottom=0, height=.07, legend_frameon=legend_frameon, clip_patch_kwargs=dict(edgecolor='k', lw=.5), tbl_chr = tbl_chromatin)

	if fn:
		if re.search('\.eps', fn):
			fig.savefig(fn, format = 'eps', dpi=save_pdi, bbox_inches='tight')
		else:
			fig.savefig(fn, jpeg_quality=100, dpi=save_dpi, bbox_inches='tight')


fig_assemble(fn = chromplot_filename)
