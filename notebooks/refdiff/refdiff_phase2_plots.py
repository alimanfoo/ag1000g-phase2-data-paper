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
_data_chromatin = b"""CHX	 chro	X	   20009764		24393108
CH2R	chro	2R	  58984778		61545105
CH2L	chro	2L	  1	   2431617
PEU2L   chro	2L	  2487770 5042389
IH2L	chro	2L	  5078962 5788875
IH3R	chro	3R	  38988757		41860198
CH3R	chro	3R	  52161877		53200684
CH3L	chro	3L	  1	   1815119
PEU3L   chro	3L	  1896830 4235209
IH3L	chro	3L	  4264713 5031692
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


#ag1k_dir = '/kwiat/vector/ag1000g/release'
# On Liverpool server:
ag1k_dir = '/home/elucas/data'
from ag1k import phase2_ar1
phase2_ar1.init(os.path.join(ag1k_dir, 'phase2.AR1'))
chroms = list(phase2_ar1.callset_phased.keys())
colours = phase2_ar1.pop_colors
colours['An. gambiae'] = 'red'
colours['An. coluzzii'] = 'blue'

subpops = phase2_ar1.df_samples.groupby('population').indices
subpops_ix = {k: list(v) for k, v in subpops.items()}
species = phase2_ar1.df_samples.groupby('m_s').indices
species['An. gambiae'] = species.pop('S')
species['An. coluzzii'] = species.pop('M')
species['Unknown'] = species.pop('M/S')
species_ix = {k: list(v) for k, v in species.items()}

window_size = 100000
dxy_by_window = dict()

for chrom in chroms:
	dxy_by_window[chrom] = pd.read_csv('refdiff_phase2_' + chrom + '_table.csv', sep = '\t')

# Now plot the dxy along the chromsomomes
def plot_dxy(chrom, ax, dxy, pops=None, offset=5, legend_frameon=False, **kwargs):
	if chrom == '2R':
		sns.despine(ax=ax, bottom=True, offset=offset)
		ax.set_yticks([0, 0.005, 0.01, 0.015], minor=False)
		ax.set_yticks(np.arange(0, 0.015, 0.001), minor=True)
		ax.set_yticklabels([0, 0.005, 0.01, 0.015])
		ax.set_ylabel('Dxy', ha='left', va='top', rotation=0)
		ax.yaxis.set_label_coords(-0.25, 1.3)
	else:
		sns.despine(ax=ax, left=True, bottom=True, offset=offset)
		ax.set_yticks([])
		ax.set_title(chrom)

	# If pops is None, plot all of the populations
	if pops == None:
		pops = list(dxy['2L'].columns)

	# Get the x axis 
	x = dxy[chrom].index
	# Create the handles for each population
	h = []
	for p in pops:
		h.append(ax.plot(x, dxy[chrom][p], color=colours[p], label = p))

	# Draw the legend
	if chrom == 'X':
		handles = []
		for this_h in h[::-1]:
			handles.append(plt.plot([0, 1], [0, 0], color = this_h[0].get_c(), lw=2, label=this_h[0].get_label())[0])
		ax.legend(handles=handles, bbox_to_anchor=(1, 1.5), loc='upper left', frameon=legend_frameon, fancybox=False, 
				  title='Species', labelspacing=.1);

	ax.set_ylim(0, 0.015)
	ax.set_xlim(0, max(x))
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


def fig_assemble(populations=None, fw=4.3, fn=None, dpi=150, save_dpi=200, legend_frameon=False, fh=0.8):
	figsize = fw, fh
	fig = plt.figure(figsize=(fw, fh), dpi=dpi)
	fig_linear_genome(plot_dxy, phase2_ar1.genome_agamp3, chromomsomes = chroms, fig=fig, bottom=0.1, height=.7, legend_frameon=legend_frameon, dxy = dxy_by_window, pops = populations)
	fig_linear_genome(plot_heterochromatin, phase2_ar1.genome_agamp3, fig=fig, bottom=0, height=.07, legend_frameon=legend_frameon, clip_patch_kwargs=dict(edgecolor='k', lw=.5), tbl_chr = tbl_chromatin)

	if fn:
		if re.search('\.eps', fn):
			fig.savefig(fn, format = 'eps', dpi=save_pdi, bbox_inches='tight')
		else:
			fig.savefig(fn, jpeg_quality=100, dpi=save_dpi, bbox_inches='tight')


fig_assemble(['BFcol', 'AOcol', 'GHcol'], fn = 'refdiff_phase2_col.jpg')
fig_assemble(['CMgam', 'UGgam', 'BFgam'], fn = 'refdiff_phase2_gam.jpg')
fig_assemble(['GM', 'GW', 'KE'], fn = 'refdiff_phase2_unk.jpg')
fig_assemble(['An. gambiae', 'An. coluzzii'], fn = 'refdiff_phase2_species.jpg')

