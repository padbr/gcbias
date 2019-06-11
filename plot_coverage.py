#!/usr/bin/python

# Author: Patrick Denis Browne
# e-mail: pdbr@plen.ku.dk
# Licensed under the terms of the MIT license
# see: https://github.com/padbr/gcbias/LICENSE for details

import sys
import re
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st


usage = """\
plot_coverage.py -i infile.seqteck -o outfile.svg -norm_pct <int> \\
-min_vals <int> -yaxis_min <float> -yaxis_max<float> -xaxis_min \\
<float> -xaxis_max <float> -title <string> log_scale y_errs avg_gc

Required Arguments:
  -i          input file

Optional Arguments:
  -o          output file: default is to append '.svg' to infile name
  -norm_pct   integer: default is to normalize by total read abundance
  -min_vals   integer: minimum number of windows in any give GC bracket
                       defaults to 1
  -yaxis_min  float: default to matplotlib decision
  -yaxis_max  float: default to matplotlib decision
  -xaxis_min  float: default to matplotlib decision
  -xaxis_max  float: default to matplotlib decision
  -title      title to put on plot: if '# Title = <string>' tag is in
              infile, that will be the default. Otherwise, no title
              will be included
  log_scale   plot y-axis on a log scale: default to straight scale
  y_errs      show y-error bars (+- 1 stdev)
  avg_gc      include a vertical line showing average gc content -
              requires '# Background GC = <float>' tag in infile

The input file format consists of 3 types of lines:
1. A comment: begins with '#' and has no '=' sign
2. A tag: begins with '#' and has and '=' sign
3. A list of tab-delimited integers
Lines of type 1. are ignored by this script.
Lines of type 2. are used to store metadata:
    'Background GC = <some float>' and 'Title = <some string>' are
    explicitly used by this script - to add a vertical red line
    indicating average GC content and to add a title, respectively.
    Other metadata won't be used by this script, but can be useful to
    record things such as window size and step size. Since these meta-
    data are not explicitly used, they could be type 1 also.
Lines of type 3. are a list of 'GC' precentages of windows in the first
    column, followed by total read depth on a per-window basis for
    every window of the genome/contigs within the value of the first
    column. The GC percentages in the first column must be unique.
    Repeated GC values will be plotted as separate points otherwise.
"""

try:
    infile = sys.argv[sys.argv.index('-i')+1]
except ValueError:
    print "Error: An input file was not correctly specified"
    print(usage)
    quit()

try:
    outfile = sys.argv[sys.argv.index('-o')+1]
    if not outfile.endswith('.svg'):
        print "Changing outfile name from %s to %s" \
            % (outfile, outfile+'.svg')
        outfile = outfile + '.svg'
except ValueError:
    outfile = infile + '.svg'

try:
    norm_pct = int(sys.argv[sys.argv.index('-norm_pct')+1])
except ValueError:
    norm_pct = None

try:
    min_vals = int(sys.argv[sys.argv.index('-min_vals')+1])
except ValueError:
    min_vals = 1

tags = {}
Title = None
pct_values = []
raw_coverages = []
for line in open(infile, 'r').readlines():
    line = line.rstrip('\n')
    if line.startswith('#'):
        if line.count('=') == 1:
            k = line.lstrip('# ').split('=')[0].rstrip().upper()
            v = line.split('=')[1].strip()
            tags[k] = v
    else:
        cols = [int(col) for col in line.split('\t')]
        if len(cols) > min_vals:
            pct_values.append(cols.pop(0))
            raw_coverages.append(cols)

if 'avg_gc' in sys.argv:
    avg_gc = float(tags['Background GC'.upper()])

if 'title'.upper() in tags.keys():
    Title = tags['TITLE']
if '-title' in sys.argv:
    Title = str(sys.argv[sys.argv.index('-title')+1])

num_windows = 0
total_coverage = 0
for pct_win in raw_coverages:
    num_windows += len(pct_win)
    total_coverage += sum(pct_win)
raw_norm_coverage = float(total_coverage)/num_windows

if norm_pct:
    assert norm_pct in pct_values, \
        'Cannot normalise to %i %%GC because it is not in %s' \
            % (norm_pct,infile)
    norm_index = pct_values.index(norm_pct)
    norm_coverage = float(sum(raw_coverages[norm_index])) \
        / len(raw_coverages[norm_index])
else:
    norm_coverage = float(raw_norm_coverage)

normalised_coverages = []
for pct_win in raw_coverages:
    normalised_coverage = []
    for i in range(len(pct_win)):
        normalised_coverage.append(float(pct_win[i])/norm_coverage)
    normalised_coverages.append(normalised_coverage)

averaged_coverages = [np.mean(norm_pcts) for norm_pcts in normalised_coverages]
standard_deviations = [np.std(norm_pcts) for norm_pcts in normalised_coverages]
num_data_pts = [len(norm_pcts) for norm_pcts in normalised_coverages]
max_num_data_pts = max(num_data_pts)

print "Percent GC\t%s" % ('\t'.join([str(pct) for pct in pct_values]))
print "Averaged coverage\t%s" % ('\t'.join([str(avg_cvg)
    for avg_cvg in averaged_coverages]))
print "Standard deviation\t%s" % ('\t'.join([str(stddev)
    for stddev in standard_deviations]))
print "Number of points\t%s" % ('\t'.join([str(i) for i in num_data_pts]))

fig = plt.figure(num=1)
cmap = plt.get_cmap('Blues')
ax1 = fig.add_subplot(111)
fig.canvas.draw()
plt.xlabel('Percent GC')
plt.ylabel('Relative Coverage')

cols = []
for num_data_pt in num_data_pts:
    cols.append(cmap(math.log(float(num_data_pt))
        / math.log(float(max_num_data_pts)))[0:3])

for i in range(len(averaged_coverages)):
    if '-yaxis_min' in sys.argv:
        if averaged_coverages[i] \
            < float(sys.argv[sys.argv.index('-yaxis_min')+1]):
            continue
    if '-yaxis_max' in sys.argv:
        if averaged_coverages[i] \
            > float(sys.argv[sys.argv.index('-yaxis_max')+1]):
            continue
    if len(raw_coverages[i]) < 5:
        continue
    ax1.plot([pct_values[i]], [averaged_coverages[i]],
        color=cols[i], marker='o', linestyle='')
    if 'y_errs' in sys.argv:
        ax1.errorbar(pct_values[i], averaged_coverages[i],
            yerr=standard_deviations[i], ecolor='k')

if '-xaxis_min' in sys.argv:
    min_xplot = float(sys.argv[sys.argv.index('-xaxis_min')+1])
else:
    min_xplot = min(pct_values)
if '-xaxis_max' in sys.argv:
    max_xplot = float(sys.argv[sys.argv.index('-xaxis_max')+1])
else:
    max_xplot = max(pct_values)
ax1.plot([min_xplot,max_xplot], [1,1], 'b--')

if '-yaxis_min' in sys.argv:
    plt.ylim(ymin=float(sys.argv[sys.argv.index('-yaxis_min')+1]))
if '-yaxis_max' in sys.argv:
    plt.ylim(ymax=float(sys.argv[sys.argv.index('-yaxis_max')+1]))

if '-xaxis_min' in sys.argv:
    plt.xlim(xmin=float(sys.argv[sys.argv.index('-xaxis_min')+1]))
if '-xaxis_max' in sys.argv:
    plt.xlim(xmax=float(sys.argv[sys.argv.index('-xaxis_max')+1]))

if 'avg_gc' in sys.argv:
    ax1.plot([avg_gc,avg_gc], [ax1.get_ylim()[0],ax1.get_ylim()[1]], 'g-')

if '-norm_pct' in sys.argv:
    ax1.plot([ax1.get_xlim()[0],ax1.get_xlim()[1]],
        [float(raw_norm_coverage)/norm_coverage,
         float(raw_norm_coverage)/norm_coverage
        ], 'r--')

if 'log_scale' in sys.argv:
    ax1.set_yscale("log", nonposy='clip')

if Title:
    plt.title(Title)

plt.savefig(outfile)
