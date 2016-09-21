#!/usr/bin/env python
"""
Usage:
    postprocess_tasks_MPI.py input.dat output.dat

Post-processes a thread info file (created by a run of SWIFT with the -y flag)
into a human readable format. That is convert task numbers into names and
tic/toc values into offset milliseconds from the start (in their rank) and the
number of milliseconds taken by the task.

This file is part of SWIFT.

Copyright (C) 2016 Peter W. Draper (p.w.draper@durham.ac.uk)
All Rights Reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import pylab as pl
import sys

#  Tasks and subtypes. Indexed as in tasks.h, must be kept up to date.
TASKTYPES = ["none", "sort", "self", "pair", "sub_self", "sub_pair", "init",
             "ghost", "extra_ghost", "kick", "kick_fixdt", "send", "recv",
             "grav_gather_m", "grav_fft", "grav_mm", "grav_up",
             "grav_external", "cooling", "count"]
SUBTYPES = ["none", "density", "gradient", "force", "grav", "tend", "count"]

#  Show docs if help is requested.
if len( sys.argv ) == 2 and ( sys.argv[1][0:2] == "-h" or sys.argv[1][0:3] == "--h" ):
    from pydoc import help
    help( "__main__" )
    sys.exit( 0 )

#  Handle command-line.
if len( sys.argv ) != 3:
    print "Usage: ", sys.argv[0], "input.dat output.dat"
    sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]

#  Read input.
data = pl.loadtxt( infile )

#  First tic, need one for each rank.
start = {}
nranks = int(max(data[:,0])) + 1
for rank in range(nranks):
    rdata = data[data[:,0] == rank]
    start[rank] = float(min(rdata[:,5]))
    print "Time begins at: ", start[rank]

#  CPU clock.
CPU_CLOCK = float(data[0,:][-1])
print "CPU frequency:", CPU_CLOCK

#  Process.
with open(outfile,'w') as outf: 
    outf.write("# rank rid itype subitype pair tic toc ci-count cj-count " +
               "ci-gcount cj-count flags start end dt tasktype subtype\n");
    for line in data:
        rank = int(line[0])
        task = int(line[2])
        subtask = int(line[3])
        tic = float(line[5])
        toc = float(line[6])

        for word in line:
            outf.write(" " + str(int(word)))
        
        outf.write(" " + str((tic-start[rank]) / CPU_CLOCK * 1000.0))
        outf.write(" " + str((toc-start[rank]) / CPU_CLOCK * 1000.0))
        outf.write(" " + str((toc - tic) / CPU_CLOCK * 1000.0))
        outf.write(" " + TASKTYPES[task])
        outf.write(" " + SUBTYPES[subtask])
        outf.write("\n")

sys.exit(0)
