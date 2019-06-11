#!/usr/bin/python2.7

# Author: Patrick Denis Browne
# e-mail: pdbr@plen.ku.dk
# Licensed under the terms of the MIT license
# see: https://github.com/padbr/gcbias/LICENSE for details

import os
from sys import argv
from Bio import SeqIO

Usage = """
extract_long_contigs.py -i <infile.fa> [options]
Takes all sequence records (contigs) of sufficient length from an input
fasta file and writes them each, one by one to individual output files.

REQUIRED:
-i <infile.fasta>   Input file of sequences. This must be in fasta format.

OPTIONAL:
-o <outdir>         Output directory. Use only alphanumeric characters
                    and underscores to be safe.
                    Default is the current directory.

-l <length>         The minimum size of fasta sequence to be kept. The
                    option specified must be an integer. Default = 10000

-s <.fa>            The file extension to use. The default is '.fa'.
                    Other sane values might be '.fasta', '.ffn', ... etc

-y                  Assume yes for the prompt about whether or not it's
                    okay to pollute the output directory with a lot of
                    files.
"""

for arg in argv:
    if arg in ['-h', '--h', '-help', '--help']:
        print(Usage)
        quit()

for arg in argv:
    if arg.startswith('-'):
        if not arg.lstrip('-') in ['i', 'o', 'l', 's', 'y']:
            print "Error: '%s' was not a recognised argument" % arg
            print(Usage)
            quit()

if not '-i' in argv:
    print """\
Error: A fasta format input file must be specified with '-i <input.fa>'.
For more help, use option '--help' to see a full help message.
"""
    quit()
infile = argv[argv.index('-i')+1]

try:
    outdir = argv[argv.index('-o')+1]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
except ValueError:
    outdir = os.getcwd().split('/')[-1]

try:
    minlen = int(argv[argv.index('-l')+1])
except ValueError:
    minlen = 10000

try:
    fnsuffix = argv[argv.index('-s')+1]
except ValueError:
    fnsuffix = '.fa'

recs = [rec for rec in SeqIO.parse(infile,'fasta')
        if len(rec) >= minlen]

if len(recs) == 0:
    print "No records in %s were longer than %i.\nNothing to do." \
        % (infile, minlen)
    quit()

if '-y' in argv:
    pass
else:
    print "%i fasta files will be written to %s.\nIs this ok? (Y/N)" \
        % (len(recs), outdir)
    if not raw_input().lower().startswith('y'):
        quit()

if '-o' in argv:
    os.chdir(outdir)
for rec in recs:
    outname = str(rec.id) + fnsuffix
    _ = SeqIO.write([rec],outname,'fasta')
