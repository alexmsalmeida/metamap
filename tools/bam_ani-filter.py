#!/usr/bin/env python

import pysam
import sys

if len(sys.argv) == 1:
    print("usage: script.py in.bam id_thresh[%] cov_thresh[%]")
    sys.exit(1)

samfile = pysam.AlignmentFile(sys.argv[1], "rb")
outfile = pysam.AlignmentFile("-", "w", template=samfile)
for read in samfile.fetch():
    try:
        ani = (read.query_alignment_length-read.get_tag("NM"))/float(read.query_alignment_length)*100
        cov = read.query_alignment_length/float(read.query_length)*100
        if ani >= float(sys.argv[2]) and cov >= float(sys.argv[3]):
            outfile.write(read)
    except:
        continue
