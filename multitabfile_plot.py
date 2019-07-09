#!/usr/bin/python2.7

# Author: Patrick Denis Browne
# e-mail: pdbr@plen.ku.dk
# Licensed under the terms of the MIT license
# see: https://github.com/padbr/gcbias/LICENSE for details

Usage="""\
A one off to plot to look at log transformed coverage data between
multiple organisms/projects.

The format of the input files is a tab-delimited one, with percent
GC in the first column, and coverage of windows with that GC in the
following columns. All numbers in this must be integers. Lines in the
input file may begin with a '#', in which case they will be ignored.
Comments cannot be placed in the end of a line without raising, most
likely, a ValueError.

For each file, a separate, contrasting colour will be used to plot.
All coverages will be divided by the average of the 49% GC window.
For now, the only way to change this behaviour is to change the 49 in
the line below (`norm_pct = 49') to another suitable integer.
The coverages, after normalisation using the average of the 49% window,
will then be log transformed and averages and standard deviations for
plotting will be calculated.
All input files will be plotted on one space, and a binomial will be
fitted to all of the combined data.

If any of the genomes to be plotted has less than `min_points'
(default = 5) windows with the GC content specified in
`norm_pct' (default = 49) an error will be raised. The only way to
solve this is to change the `norm_pct' variable to a GC content for
which all genomes to be plotted have sufficient windows, or to lower
the value of `norm_pct', the latter option being less desireable.
For many projects, it is likely that this limitation will make it
impossible to plot all genomes together - a limitation for which I
have no intended solution.

usage: tabfile_plot.py file1.tab file2.tab file3.tab ... (max 26 files)
All input file namess must end with '.tab' (case sensitive)
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def tabbed_coverage_parse(tabfile, norm_pct, min_points, log_transform):
    lines = [
        line.rstrip('\n').split('\t')
        for line in open(tabfile, 'r').readlines()
        if not line.startswith('#')
    ]
    lines = [line for line in lines if len(line) > min_points]
    xs = [int(line[0]) for line in lines] # GC percent values
    covs = [[int(cov) for cov in covs[1:]] for covs in lines]
    if log_transform == 'yes':
        covs = [[int(cov)+1 for cov in covs[1:]] for covs in lines]
    weights = [len(cov) for cov in covs]
    norm_index = xs.index(norm_pct)
    norm_factor = np.average(covs[norm_index])
    covs = [[float(num)/norm_factor for num in cov] for cov in covs]
    if log_transform == 'yes':
        covs = [[math.log(num,10) for num in cov] for cov in covs]
    ys = [np.average(cov) for cov in covs]
    stdevs = [np.std(cov) for cov in covs]
    return(xs, ys, stdevs, weights)
    
    

########################################################################
# Define Key Constants
norm_pct = 49
min_points = 5
# List of twenty-six contrasting colours
rgb_colours = [
    (240,163,255), (0,117,220), (153,63,0), (76,0,92), (25,25,25), (0,92,49),
    (43,206,72), (255,204,153), (128,128,128), (148,255,181), (143,124,0),
    (157,204,0), (194,0,136), (0,51,128), (255,164,5), (255,168,187),
    (66,102,0), (255,0,16), (94,241,242), (0,153,143), (224,255,102),
    (116,10,255), (153,0,0), (255,255,128), (255,255,0), (255,80,5),
]
rgb_colours = [tuple(float(col[i])/255 for i in range(3))
    for col in rgb_colours]
log_transform = 'yes'
weighted_binomial = 'yes'
decplcs = 6 # Demical places for rounding in the fitted binomial equation
########################################################################

infiles = [str(File) for File in sys.argv if str(File).endswith('.tab')]
xs = []
ys = []
stdevs = []
colours = []
weights = []
fig = plt.figure(num=1)
ax = plt.subplot(111)
for i in range(len(infiles)):
    file_xs, file_ys, file_stdevs, file_weights \
    = tabbed_coverage_parse(infiles[i], norm_pct, min_points, log_transform)
    file_color = rgb_colours[i]
    ax.plot(file_xs, file_ys,marker='o', markerfacecolor=file_color,
            linestyle='', label=infiles[i].split('_miseq_')[0])
    for j in range(len(file_xs)):
        xs.append(file_xs[j])
        ys.append(file_ys[j])
        stdevs.append(file_stdevs[j])
        colours.append(file_color)
        weights.append(file_weights[j])

for i in range(len(xs)):
    ax.errorbar(xs[i], ys[i], stdevs[i], ecolor='k')
if weighted_binomial == 'yes':
    fitted_binomial = np.polyfit(xs, ys, deg=2, w=weights)
else:
    fitted_binomial = np.polyfit(xs, ys, deg=2)
eq_print = ('trendline: '
            + str(round(fitted_binomial[0],decplcs))
            + 'x^2 + '
            + str(round(fitted_binomial[1],decplcs))
            + 'x + '
            + str(round(fitted_binomial[2],decplcs))
           )
ylims = ax.get_ylim()
xlims = ax.get_xlim()
fit_xs = []
fit_ys = []
for i in np.arange(xlims[0],xlims[1],0.1):
    j = fitted_binomial[0]*(i**2) + fitted_binomial[1]*i + fitted_binomial[2]
    if j >= ylims[0] and j <= ylims[1]:
        fit_xs.append(i)
        fit_ys.append(j)
ax.plot(fit_xs, fit_ys, linestyle='-', linewidth=2, color='k',
        marker='', label=eq_print)
legend = ax.legend(loc='lower right', shadow=True, fontsize='x-large')
ax.set_xlabel('GC content(%)', size=14)
ax.set_ylabel('Average normalized log transformed coverage', size=18)
ax.set_title('MiSeq', size=20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
plt.show()
