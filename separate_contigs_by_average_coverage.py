#!/usr/bin/python2.7

# Author: Patrick Denis Browne
# e-mail: pdbr@plen.ku.dk
# Licensed under the terms of the MIT license
# see: https://github.com/padbr/gcbias/LICENSE for details

import sys
import math
import gzip
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio import SeqIO

usage = """separate_contigs_by_average_coverage.py -d aln.depth -ref
contigs.fasta -thresh <float> -norm_out <file> -outliers <file>
-revise

Assesses coverages of contigs following genome assembly and provides a
means of separating normally covered from over covered contigs. In the
plot, contigs with sufficient length and normal coverage will be
represented in green. Contigs with sufficient length and abberant
coverage will be red. Contigs with insufficient length will be grey.
Because GC-content can be a major source of coverage bias, the plots are
presented in 3D (contig length VS coverage VS percent GC content).

This is intended to be a supervised method for distinguishing abberant
contigs in an assembly. Unresolved repeats, plasmids, and phages are
examples of elements that may result in abberantly covered contigs. It
is essential that the "supervised" element of this script is respected.
Otherwise, undercovered abberants may also cause valid contigs to be
erroneously removed.



Required Arguments:
  -d            <file>: file with coverage depths in the same format as
                'samtools depth -a' output
  -ref          <file>: fasta file with the reference
                genome(s)/contig(s) specified in the first column of -d
  -thresh       <float>: The modified Z-score to use as a threshold.
                Observations with a modified Z-score (based on median
                absolute deviation) greater than this value will be
                classified as outliers.

Optional Arguments:
  -min_size     <int>: Default=10,000. Minimum size of contigs to include
                in the analysis. Normal and abberant coverage will be
                assessed after the removal of small contigs.
  -norm_out     <file>: A file to write contigs with normal coverage.
                Fasta format.
  -outliers     <file>: A file to write contigs with abnormal coverage.
                Fasta format. Does not include small contigs.
  -revise       Including this flag will give an interactive prompt to
                revise the threshold during runtime. This will allow for
                a visual examination of the contigs (coloured red) that
                will be treated as over the threshold.
  

If neither '-norm_out' or '-outliers' are specified in the arguments,
this script will effectively be used for data exploration - examining
the sizes and coverages of contigs.

Nucleotide cases are ignored in the reference genome.
The first column of such a depth file must be an exact match to a
sequence identifier in the reference fasta.
The genomic coordinates in the depth file must be 1-based.
"""


########################################################################
#                         Define Functions                             #
########################################################################


def NGS_depth(annotation, seq_format,depth):
    """
    This is basically intended to act as an extension to BioPython and
    parse coverage depth per nucleotide into a biopython-like object.
    
    annotation: reference sequence annotation
                the sequence identifiers must exactly match those in
                depth
    seq_format: 'fasta' or 'genbank' or any other biopython format
                string
    depth:      a file produced from the output of:
                `samtools depth -a file.bam'
    
    returns: a list of biopython-like records
             each record contains the usual biopython annotations.
             additionally each record contains 'per_letter_depth'
             annotations, which track the sequencing depth at each
             position in the record
    """
    records = [record for record in SeqIO.parse(annotation, seq_format)]
    record_indices = {}
    for i in range(len(records)):
        records[i].annotations['per_letter_depth'] \
            = [0 for nucl in str(records[i].seq)]
        record_indices[str(records[i].id)] = i
    
    if depth.endswith('.gz'):
        depth_h = gzip.open(depth, 'r')
    else:
        depth_h = open(depth, 'r')
    for line in depth_h.readlines():
        line = line.rstrip('\n')
        cols = line.split('\t')
        assert len(cols) == 3, \
            "Exactly three columns were not found in line:\n%s" % line
        rec_id = str(cols[0])
        z_pos = int(cols[1]) - 1
        coverage = int(cols[2])
        rec_index = record_indices[rec_id]
        records[rec_index].annotations['per_letter_depth'][z_pos] = coverage
    return(records)


def is_outlier(avg_coverages, thresh=5.0):
    """
    Returns a boolean array with True if points are outliers and
    False otherwise.
    
    Parameters
    avg_coverages: A one (or multi) dimensional list.
    thresh:        A modified z-score used as a threshold. Higher
                   values mean lower stringency.
    
    Taken almost verbatim from:
    github.com/joferkington/oost_paper_code/blob/master/utilities.py
    """
    points = np.array(avg_coverages)
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points-median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745*diff/med_abs_deviation
    return modified_z_score > thresh


def filter_ngs_coverage_outliers(ngs_depth_objs, threshold):
    average_coverages = [float(sum(ngs_depth_obj.annotations['per_letter_depth'
        ]))/len(ngs_depth_obj.seq) for ngs_depth_obj in ngs_depth_objs]
    outliers = is_outlier(average_coverages, thresh=threshold)
    norm_coverage_objs = [ngs_depth_objs[i] for i in range(len(ngs_depth_objs))
        if outliers[i] == False]
    high_coverage_objs = [ngs_depth_objs[i] for i in range(len(ngs_depth_objs))
        if outliers[i] == True]
    return(norm_coverage_objs, high_coverage_objs, outliers)


def pct_gc(ngs_depth_objs):
    gcs = []
    for rec in ngs_depth_objs:
        seq = str(rec.seq).lower()
        atgc_count = 0
        gc_count = 0
        for nucl in seq:
            if nucl in ['a', 't', 'g', 'c']:
                atgc_count += 1
            if nucl in ['g', 'c']:
                gc_count += 1
        if atgc_count == 0:
            gcs.append(float(100))
        else:
            gc_pct = float(gc_count*100) / atgc_count
            gcs.append(gc_pct)
    return(gcs)


if __name__ == '__main__':
    ########################################################################
    #                   Parse Command Line Options                         #
    ########################################################################
    try:
        depth_f = sys.argv[sys.argv.index('-d')+1]
        ref_f = sys.argv[sys.argv.index('-ref')+1]
    except ValueError:
        print("\nERROR: Specify input files with: -d and -ref\n\n")
        print(usage)
        quit()
    try:
        thresh = float(sys.argv[sys.argv.index('-thresh')+1])
    except ValueError:
        print("\nERROR: Sepcify a threshold: e.g: -thresh = 10.0\n\n")
        print(usage)
        quit()
    if '-norm_out' in sys.argv:
        norm_f = sys.argv[sys.argv.index('-norm_out')+1]
    if '-outliers' in sys.argv:
        outliers_f = sys.argv[sys.argv.index('-outliers')+1]
    if '-revise' in sys.argv:
        revise = 'yes'
    else:
        revise = 'no'
    if '-min_size' in sys.argv:
        min_size = int(sys.argv[sys.argv.index('-min_size')+1])
    else:
        min_size = 10000
    
    
    ########################################################################
    #             Parse Raw Data and Filter Abberant Contigs               #
    ########################################################################
    ngs_records = NGS_depth(ref_f, 'fasta', depth_f)
    raw_records = [rec for rec in ngs_records if len(rec.seq) >= min_size]
    short_records = [rec for rec in ngs_records if len(rec.seq) < min_size]
    gc_pcts = pct_gc(raw_records)
    short_gcs = pct_gc(short_records)
    short_xs = [len(rec.seq) for rec in short_records]
    short_zs = [float(sum(rec.annotations['per_letter_depth'])) / len(rec.seq)
        for rec in short_records]
    print "Contig\tLength\tPercent GC\tAverage Coverage"
    for i in range(len(raw_records)):
        recid = str(raw_records[i].id)
        reclength = str(len(raw_records[i].seq))
        gc_content = "%.2f" % gc_pcts[i]
        avg_coverage = float(sum(raw_records[i].annotations['per_letter_depth']
            )) / len(raw_records[i].seq)
        avg_coverage = "%.2f" % avg_coverage
        print "%s\t%s\t%s\t%s" % (recid, reclength, gc_content, avg_coverage)
        
    if revise == 'no':
        norm_cov_contigs, high_cov_contigs, outliers \
            = filter_ngs_coverage_outliers(raw_records, thresh)
        norm_xs = [len(record.seq) for record in norm_cov_contigs]
        norm_zs = [float(sum(record.annotations['per_letter_depth']))
            / len(record.seq) for record in norm_cov_contigs]
        high_xs = [len(record.seq) for record in high_cov_contigs]
        high_zs = [float(sum(record.annotations['per_letter_depth']))
            / len(record.seq) for record in high_cov_contigs]
        norm_ys = [gc_pcts[i] for i in range(len(raw_records))
            if outliers[i] == False]
        high_ys = [gc_pcts[i] for i in range(len(raw_records))
            if outliers[i] == True]
        fig = plt.figure(figsize=(9,6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(short_xs, short_gcs, short_zs, marker='o', color='grey')
        ax.scatter(norm_xs, norm_ys, norm_zs, marker='o', color='green')
        ax.scatter(high_xs, high_ys, high_zs, marker='o', color='red')
        ax.set_xlabel('Contig length (bp)')
        ax.set_zlabel('Average Coverage')
        ax.set_ylabel('GC content (percent)')
        plt.show()
    
    while revise == 'yes':
        norm_cov_contigs, high_cov_contigs, outliers \
            = filter_ngs_coverage_outliers(raw_records, thresh)
        norm_xs = [len(record.seq) for record in norm_cov_contigs]
        norm_zs = [float(sum(record.annotations['per_letter_depth']))
            / len(record.seq) for record in norm_cov_contigs]
        high_xs = [len(record.seq) for record in high_cov_contigs]
        high_zs = [float(sum(record.annotations['per_letter_depth']))
            / len(record.seq) for record in high_cov_contigs]
        norm_ys = [gc_pcts[i] for i in range(len(raw_records))
            if outliers[i] == False]
        high_ys = [gc_pcts[i] for i in range(len(raw_records))
            if outliers[i] == True]
        fig = plt.figure(figsize=(9,6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(short_xs, short_gcs, short_zs, marker='o', color='grey')
        ax.scatter(norm_xs, norm_ys, norm_zs, marker='o', color='green')
        ax.scatter(high_xs, high_ys, high_zs, marker='o', color='red')
        ax.set_xlabel('Contig length (bp)')
        ax.set_zlabel('Average Coverage')
        ax.set_ylabel('GC content (percent)')
        plt.show()
        plt.clf()
        fig = ''
        ax = ''
        revise_update = raw_input("Revise threshold? y/n:").lower()
        if revise_update.startswith('y'):
            thresh = float(raw_input("Input new threshold value:"))
        else:
            revise = 'no'
    
    if '-norm_out' in sys.argv:
        written = SeqIO.write(norm_cov_contigs, norm_f, 'fasta')
        print "Wrote %i records to %s" % (written, norm_f)
    if '-outliers' in sys.argv:
        written = SeqIO.write(high_cov_contigs, outliers_f, 'fasta')
        print "Wrote %i records to %s" % (written, outliers_f)
