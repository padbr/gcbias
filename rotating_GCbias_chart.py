#!/usr/bin/python2.7

# Author: Patrick Denis Browne
# e-mail: pdbr@plen.ku.dk
# Licensed under the terms of the MIT license
# see: https://github.com/padbr/gcbias/LICENSE for details

import os
import sys
import math
import subprocess
import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

###################### DEFINE SOME CONSTANTS ###########################
xlims = [14.25, 85.0]
ylims = [14.25, 85.0]
zlims = [-0.75, 0.75] # This could be too restrictive in some cases
camera_angles = [(-90,0), (-1,0), (-1,5), (-90,5),
                 (-90,10), (-1,10), (-1,15), (-90,15)]
# The camera angles are given as (azimuth, elevation) in degrees.
########################################################################

Usage = """rotating_chart.py -i <infile.p> -p <prefix>

infile.p    A pickle dump file with x, y and z values to plot as lists.
prefix      A prefix to name all of the image files with.
"""

# Grab command line arguments
assert '-i' in sys.argv, \
"An infile must be specified with '-i <infile>'\n%s" % (Usage)
assert '-p' in sys.argv, \
"A prefix must be specified with '-p <prefix>'\n%s" % (Usage)
infile = open(sys.argv[sys.argv.index('-i')+1],'rb')
prefix = sys.argv[sys.argv.index('-p')+1]

# Make the basic figure
xs, ys, zs = pickle.load(infile)
fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot([min(xs),max(xs)], [min(xs),max(xs)], [0,0], c='r', marker='.')
ax.scatter(xs, ys, zs, c='g', marker='o')
ax.set_xlabel('Numerator GC%')
ax.set_ylabel('Denominator GC%')
ax.set_zlabel('log(Coverage Ratio)')
ax.set_xlim(xlims)
ax.set_ylim(ylims)
ax.set_zlim(zlims)

# Rotate the figure around and produce images
suffix = 0
ax.view_init(azim=camera_angles[0][0], elev=camera_angles[0][1])
plt.savefig("%s_%s.png" % (prefix, str(suffix).zfill(4)))
for i in range(len(camera_angles)-1):
    start_angles = camera_angles[i]
    stop_angles = camera_angles[i+1]
    assert start_angles[0] == stop_angles[0] \
           or start_angles[1] == stop_angles[1]
    if start_angles[0] != stop_angles[0]:
        if start_angles[0] < stop_angles[0]:
            for j in xrange(start_angles[0], stop_angles[0], 1):
                suffix += 1
                ax.view_init(azim=j, elev=start_angles[1])
                plt.savefig("%s_%s.png" % (prefix, str(suffix).zfill(4)))
        if start_angles[0] > stop_angles[0]:
            js = [k for k in xrange(stop_angles[0], start_angles[0], 1)]
            js.reverse()
            for j in js:
                suffix += 1
                ax.view_init(azim=j, elev=start_angles[1])
                plt.savefig("%s_%s.png" % (prefix, str(suffix).zfill(4)))
    else:
        if start_angles[1] < stop_angles[1]:
            for j in xrange(start_angles[1], stop_angles[1], 1):
                suffix += 1
                ax.view_init(azim=start_angles[0], elev=j)
                plt.savefig("%s_%s.png" % (prefix, str(suffix).zfill(4)))
        if start_angles[1] > stop_angles[1]:
            js = [k for k in xrange(stop_angles[1], start_angles[1], 1)]
            js.reverse()
            for j in js:
                suffix += 1
                ax.view_init(azim=start_angles[0], elev=j)
                plt.savefig("%s_%s.png" % (prefix, str(suffix).zfill(4)))

command = "ffmpeg -f image2 -r 24 -i %s_%s04d.png -c:v libx264 \
-pix_fmt yuv420p %s.mp4" % (prefix,'%',prefix)
print "Done. Just tidying up a bit now."
subprocess.check_call(command, shell=True)
command = "tar -capf %s_pngs.tar.gz %s_*.png" % (prefix, prefix)
subprocess.check_call(command, shell=True)
command = "rm %s_*.png" % (prefix)
subprocess.check_call(command, shell=True)
