#!/bin/bash

# Author: Patrick Denis Browne
# e-mail: pdbr@plen.ku.dk
# Licensed under the terms of the MIT license
# see: https://github.com/padbr/gcbias/LICENSE for details

# It is assumed that you have sequencing reads which have been cleaned
# for adapter sequences. You may have single reads, paired reads, or a
# mixture of single and paired reads.
# 
# The second major assumption (prerequisite) is that these reads have
# been assembled, using a genome assembler, or that an assembly is
# available. The assembly can be very good, basically representing a 
# closed genome, or it could be in hundreds of contigs.
#
# It is aslo assumed that you have bwa mem and the samtools package
# installed and callable on your system.
# Read the manual pages and online help for bwa and samtools for more
# information about these tools if you need it.

help=

if [ "$1" == "-h" ]; then
echo "
Options:
    -g  Reference sequence in fasta format
    -f  Forward reads in fastq format
    -r  Reverse reads in fastq format
    -s  Unpaired reads in fastq format
    -t  Number of threads to use (default = 2)
    -h  Show this message

There are three options for reads:
(i):   Use paired (forward AND reverse) reads without unpaired reads
(ii):  Use unpaired reads without paired reads
(iii): Use paired reads and unpaired reads
"
exit
fi

paired='False'
unpaired='False'
nt='2'
showhelp='NO'
while getopts ":g:f:r:s:t:" opt; do
case $opt in

g)
ref=${OPTARG}
;;

f)
forward=${OPTARG}
paired='True'
;;

r)
reverse=${OPTARG}
;;

s)single=${OPTARG}
unpaired='True'
;;

t)
nt=${OPTARG}
;;

\?)
echo "Invalid option: ${OPTARG}"
exit
esac
done

if [ ${paired} == 'False' ] && [ ${unpaired} == 'False' ]; then
echo "
Some reads must be specified. Run:
$(basename "$0") -h
for help"
exit
fi

# Run this the first time you need to map reads to this assembly only
indref=${ref}".bwt"
if ! [ -f $indref ]; then
echo "Indexing Reference"
bwa index ${ref}
fi

# Map the merged reads to the assembly
if [ ${unpaired} == 'True' ]; then
bwa mem -t ${nt} ${ref} ${single} | samtools view -@ ${nt} -b - | \
samtools sort -@ ${nt} - -o merged.bam --reference ${ref} -T temp
fi

# Map the non-merged reads
if [ ${paired} == 'True' ]; then
bwa mem -t ${nt} ${ref} ${forward} ${reverse} | \
samtools view -@ ${nt} -b - | \
samtools sort -@ ${nt} - -o pairs.bam --reference ${ref} -T temp
fi

# Combine merged and non-merged read mappings, if necessary
if [ ${paired} == 'True' ] && [ ${unpaired} == 'True' ]; then
echo "Running in mixed mode"
samtools merge aln.bam pairs.bam merged.bam
elif [ ${paired} == 'True' ] && [ ${unpaired} == 'False' ]; then
echo "Running in paired only mode"
mv pairs.bam aln.bam
else
echo "Running in single read mode"
mv merged.bam aln.bam
fi

samtools view -@ ${nt} -b -F `echo "0x904"` aln.bam > aln.904.bam
samtools sort -@ ${nt} -o aln.904.sort.bam aln.904.bam --reference ${ref} \
-T temp
samtools index aln.904.sort.bam

# Calculate the depth at every position of every contig in ref.fa
samtools depth -a aln.904.sort.bam > aln.904.depth
