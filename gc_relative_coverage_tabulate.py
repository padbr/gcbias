#!/usr/bin/python2.7

# Author: Patrick Denis Browne
# e-mail: pdbr@plen.ku.dk
# Licensed under the terms of the MIT license
# see: https://github.com/padbr/gcbias/LICENSE for details

import sys
from Bio import SeqIO


def gc_pct(seq):
    seq = seq.upper()
    total = 0
    GC = 0
    for nucl in seq:
        if nucl in ['A', 'T', 'G', 'C']:
            total += 1
        if nucl in ['G', 'C']:
            GC += 1
    if total == 0:
        return(None)
    else:
        pct_gc = float(100*GC) / float(total)
        return(pct_gc)


usage = """Usage:
gc_relatvie_coverage_tabulate.py -d aln.depth -ref reference.fasta \\
-window <int> -step <int> -title > outfile

Required Arguments:
  -d            file with coverage depths in the same format as
                'samtools depth -a' output
  -ref          fasta file with the reference genome(s)/contig(s)
                specified in the first column of -d

Optional Arguments:
  -window       integer: default is 500: Window size in which GC
                contents and normalized coverage depths are calculated
  -step         integer: default is half the window size (rounded down):
                Amount by which the window slides along a genome/contig
  -title        string: title to put over the plot. The default is to
                have no title
  
Output is written to STDOUT

Any non ATGC nucleotides are ignored while calculating GC percent.
Nucleotide cases are ignored in the reference genome.
This script is written to parse depth files produced by
'samtools depth -a something.bam'.
The first column of such a depth file must be an exact match to a
sequence identifier in the reference fasta.
The genomic coordinates in the depth file must be 1-based.

GC content is assessed in a sliding window.
There are defaults set for window size and step size which can be
changed. The window defaults to 500. The step size defaults to half the
window size - rounded down. The window sized would perhaps ideally be
around the average insert size of the sequenced fragment (assuming
Illumina HiSeq/MiSeq/NextSeq). The step size would normally not be
larger than the window size.
"""


# Parsing the command line options
try:
    depth_f = sys.argv[sys.argv.index('-d')+1]
    ref_f = sys.argv[sys.argv.index('-ref')+1]
except ValueError:
    print("ERROR: You must use both '-d' and '-ref' options!\n\n")
    print(usage)
    quit()
try:
    win_size = int(sys.argv[sys.argv.index('-window')+1])
except ValueError:
    win_size = 500
try:
    stp_size = int(sys.argv[sys.argv.index('-step')+1])
except ValueError:
    stp_size = int(win_size)/2


records_indices = {}
records = [record for record in SeqIO.parse(ref_f, 'fasta')]
for i in range(len(records)):
    records[i].annotations['per_letter_annotations'] \
        = [0 for nucl in str(records[i].seq)]
    records_indices[str(records[i].id)] = i

bkg_gc = gc_pct(''.join([str(rec.seq) for rec in records]))
totalSeqLength = len(''.join([str(rec.seq) for rec in records]))

for line in open(depth_f, 'r').readlines():
    line = line.rstrip()
    cols = line.split('\t')
    assert len(cols) == 3, "Exactly 3 columns were not found in:\n%s" % line
    rec_id = cols[0]
    z_pos = int(cols[1]) - 1
    coverage = int(cols[2])
    records[records_indices[rec_id]].annotations[
        'per_letter_annotations'][z_pos] = coverage

gc_coverages = {}
for record in records:
    i = 0
    while i + win_size <= len(str(record.seq)):
        subseq = str(record.seq)[i:i+win_size]
        valid_nucls = 0
        for nucl in subseq:
            if nucl.upper() in ['A','T','G','C']:
                valid_nucls += 1
        if valid_nucls > 0:
            pct_gc = int(round(gc_pct(subseq),0))
            raw_coverage = sum(record.annotations[
                'per_letter_annotations'][i:i+win_size])
            if not pct_gc in gc_coverages:
                gc_coverages[pct_gc] = [raw_coverage]
            else:
                gc_coverages[pct_gc].append(raw_coverage)
        i += stp_size

print '# Background GC = %s' % str(round(bkg_gc,1))
print '# Window size = %i' % (win_size)
print '# Step size = %i' % (stp_size)
print '# Total Sequence Length = %i' % (totalSeqLength)
if '-title' in sys.argv:
    title = sys.argv[sys.argv.index('-title')+1]
    print '# Title = %s' % (title)

for k in sorted(gc_coverages.keys()):
    v = '\t'.join([str(i) for i in gc_coverages[k]])
    print '%i\t%s' % (k,v)
