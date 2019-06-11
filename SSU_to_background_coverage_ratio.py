#!/usr/bin/python2.7

# Author: Patrick Denis Browne
# e-mail: pdbr@plen.ku.dk
# Licensed under the terms of the MIT license
# see: https://github.com/padbr/gcbias/LICENSE for details

########################################################################
# Change 'infile' and 'foreground_loci' variables below as necessary

infile = 'genome.depth' # See README for the format expected in this file
foreground_loci = [(416816,421540,'complement'), (667521,672245,'complement'),
    (932359,937403,'complement'), (1049652,1054375,'complement'),
    (1104569,1109292,'complement'), (1302472,1307208,''), (1308738,1313475,''),
    (1318283,1323007,''), (1427519,1432411,'')]
########################################################################

def is_foreground(target, foreground_loci):
    for locus in foreground_loci:
        if target >= min([int(locus[0]), int(locus[1])]) \
                and target <= max([int(locus[0]), int(locus[1])]):
            return True
    return False


foreground_length = 0
foreground_coverage = 0
background_length = 0
background_coverage = 0

for line in open(infile, 'r').readlines():
    line = line.rstrip()
    cols = line.split('\t')
    coord = int(cols[1])
    coverage = int(cols[2])
    if is_foreground(coord, foreground_loci):
        foreground_length += 1
        foreground_coverage += coverage
    else:
        background_length += 1
        background_coverage += coverage

print "The foreground had a length of %i nt and a total coverage of %i" \
    % (foreground_length, foreground_coverage)
print "The background had a length of %i nt and a total coverage of %i" \
    % (background_length, background_coverage)

foreground_specific_coverage \
    = float(foreground_coverage) / float(foreground_length)
background_specific_coverage \
    = float(background_coverage) / float(background_length)
fg_bg_ratio \
    = float(foreground_specific_coverage) / float(background_specific_coverage)
print "The ratio of foreground to background coverage is %s" \
    % (str(fg_bg_ratio))
