#!/usr/bin/python2.7

# Author: Patrick Denis Browne
# e-mail: pdbr@plen.ku.dk
# Licensed under the terms of the MIT license
# see: https://github.com/padbr/gcbias/LICENSE for details

import os
import math
import copy
import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class depthTab:
    def __init__(self,dFile,minWndws=5):
        """
        The format of the input files is a tab-delimited one, with
        percent GC in the first column, and coverage of windows with
        that GC in the following columns. All numbers in this must be
        integers. Lines in the input file may begin with a '#', in which
        case they will be ignored. Comments cannot be placed in the end
        of a line without raising, most likely, a ValueError.
        
        Four comments are mandotory in the input file. These must be in
        the format of `# Key = Value' where the keys are: (i) Background
        GC, (ii), Window size, (iii), Step size, and (iv) Title.
        The value for 'Background GC' must be a number and may
        optionally contain a decimal point. The values for 'Window size'
        and 'Step size' must be integers. The value for 'Title' can be
        any string - It may be used as a title for a chart later.
        
        The minWndws variable means the number of integer values that
        must follow the releveant GC content in order to be included in
        the output.
        """
        comments = [line.rstrip('\n').lstrip('# ').split('=') 
                    for line in open(dFile,'r').readlines()
                    if line.startswith('#')
                    and '=' in line and len(line.split('=')) == 2]
        self.metadata = {}
        for i in range(len(comments)):
            comments[i][0] = comments[i][0].strip()
            comments[i][1] = comments[i][1].strip()
            self.metadata[comments[i][0]] = comments[i][1]
        for k in ['Background GC', 'Window size', 'Step size', 'Title']:
            assert k in self.metadata.keys(), \
            "There is no value for %s in %s" % (k, dFile)
        self.metadata['Background GC'] = float(self.metadata['Background GC'])
        self.metadata['Window size'] = int(self.metadata['Window size'])
        self.metadata['Step size'] = int(self.metadata['Step size'])
        # 'Title' is not in self.metadata despite it being a mandatory input!
        self.gcCoverage = {}
        for line in open(dFile,'r').readlines():
            line = line.rstrip('\n')
            if not line.startswith('#'):
                cols = [int(col) for col in line.split('\t')]
                gcpct = cols.pop(0)
                if len(cols) >= minWndws:
                    self.gcCoverage[gcpct] = cols
        self.gcCoverageKeys = [pctgc for pctgc in self.gcCoverage.keys()]
        self.gcCoverageKeys.sort()
    
    def __len__(self):
        return(len(self.gcCoverageKeys))
    
    def compare_within(self):
        self.within_comparison = {}
        for i in range(len(self.gcCoverageKeys)):
            for j in range(len(self.gcCoverageKeys)):
                if not i == j:
                    pcti = self.gcCoverageKeys[i]
                    pctj = self.gcCoverageKeys[j]
                    avgi = float(sum(self.gcCoverage[pcti])) \
                           / len(self.gcCoverage[pcti])
                    avgj = float(sum(self.gcCoverage[pctj])) \
                           / len(self.gcCoverage[pctj])
                    ratioij = avgi / avgj
                    self.within_comparison[(pcti, pctj)] = ratioij



def metaGCstandard_curve_from_tabFiles(
        minWndws=5, outfile='within_comparisons.tsv', suppressPlot=False):
    """
    Makes a 3D-standard curve of GC bias in a metagenome.
    It is assumed that all *.tab files in the current directory each
    represent a contig, with certain metadata in comments ('Background 
    GC', 'Window size', 'Step size', 'Title'). Without comments it
    expects tab-delimited integers. For each line, it expects a GC value
    (rounded to nearest integer) and all following columns are coverage
    values for each genomic window with that GC-content. The relative
    coverage of windows with differing GC-contents within each contig
    yields a profile of GC-bias in the dataset, which is written to 
    'within_comparison.tsv', presented a 3D plot and exported in the
    form of a dictionary where the keys are 2-element tuples, the first
    element being the Numerator GC content, and the second element being
    the Denominator GC content, and the values being the log-transformed
    (base 10) coverage ratios)
    """
    tFiles = [File for File in os.listdir(os.getcwd())
              if File.endswith('.tab')]
    depthObjs = [depthTab(tFile, minWndws=minWndws) for tFile in tFiles]
    dropIndices = []
    for i in range(len(depthObjs)):
        if len(depthObjs[i]) == 0:
            dropIndices.append(i)
    dropIndices.reverse()
    for dropIndex in dropIndices:
        a = tFiles.pop(dropIndex)
        a = depthObjs.pop(dropIndex)
    within_comparisons = {}
    for i in range(len(depthObjs)):
        depthObjs[i].compare_within()
        for k,v in depthObjs[i].within_comparison.items():
            within_comparisons[k] = within_comparisons.get(k,[]) + [v]
    outh = open(outfile, 'w')
    header = '\t'.join(['Numerator GC', 'Denominator GC', 'Coverage Ratio',
                        'log(Coverage Ratio, 10)\n'])
    outh.write(header)
    xs = []
    ys = []
    zs = []
    for k,v in within_comparisons.items():
        avgRatio = float(sum(v)) / len(v)
        logRatio = math.log(avgRatio, 10)
        output = "%i\t%i\t%s\t%s\n" % (k[0], k[1],
                                       str(avgRatio), str(logRatio))
        outh.write(output)
        xs.append(k[0])
        ys.append(k[1])
        zs.append(logRatio)
    outh.close()
    if not suppressPlot:
        fig = plt.figure(num=1)
        ax = fig.add_subplot(111, projection='3d')
        ax.plot([min(xs),max(xs)], [min(xs),max(xs)], [0,0], c='r', marker='.')
        ax.scatter(xs, ys, zs, c='g', marker='o')
        ax.set_xlabel('Numerator GC%')
        ax.set_ylabel('Denominator GC%')
        ax.set_zlabel('log(Coverage Ratio)')
        plt.show()
    return(within_comparisons)

def metaGCcurve_pickle(minWndws=5, outfile='metaGCcurve.p'):
    '''
    Makes a 3D-plot of GC bias in a metagenome and stores the figure in 
    a pickle object, for later recall and viewing with matplotlib.
    It is assumed that all *.tab files in the current directory each
    represent a contig, with certain metadata in comments ('Background 
    GC', 'Window size', 'Step size', 'Title'). Without comments it
    expects tab-delimited integers. For each line, it expects a GC value
    (rounded to nearest integer) and all following columns are coverage
    values for each genomic window with that GC-content. The relative
    coverage of windows with differing GC-contents within each contig
    yields a profile of GC-bias in the dataset, which is written to 
    'within_comparison.tsv', presented a 3D plot and exported in the
    form of a dictionary where the keys are 2-element tuples, the first
    element being the value of the Numerator GC content, and the second
    element being the value of the Denominator GC content, and the
    values being the log-transformed (base 10) coverage ratios)
    '''
    tFiles = [File for File in os.listdir(os.getcwd()) 
              if File.endswith('.tab')]
    depthObjs = [depthTab(tFile, minWndws=minWndws) for tFile in tFiles]
    dropIndices = []
    for i in range(len(depthObjs)):
        if len(depthObjs[i]) == 0:
            dropIndices.append(i)
    dropIndices.reverse()
    for dropIndex in dropIndices:
        a = tFiles.pop(dropIndex)
        a = depthObjs.pop(dropIndex)
    within_comparisons = {}
    for i in range(len(depthObjs)):
        depthObjs[i].compare_within()
        for k,v in depthObjs[i].within_comparison.items():
            within_comparisons[k] = within_comparisons.get(k,[]) + [v]
    outh = open(outfile, 'wb')
    xs = []
    ys = []
    zs = []
    for k,v in within_comparisons.items():
        avgRatio = float(sum(v)) / len(v)
        logRatio = math.log(avgRatio, 10)
        xs.append(k[0])
        ys.append(k[1])
        zs.append(logRatio)
    coords = (xs, ys, zs)
    pickle.dump(coords, outh, 2)

if __name__ == '__main__':
    # get an interactive viaualisation of the GC-bias and a file called
    # `within_comparisons.tsv' to see the results in a spreadsheet
    wc = metaGCstandard_curve_from_tabFiles()
    # make a pickled file for python to quickly load the 3-d data later
    # by default this makes a file called metaGCcurve.p
    metaGCcurve_pickle()
    print("See 'within_comparisons.tsv' for results in spreadsheet format")
    print("The results are also encoded in 'metaGCcurve.p' for easy importing")
    
