#!/usr/bin/env python

import os
import sys
import math
import numpy

if len(sys.argv) < 3:
    print("usage: script.py bwa_depth.tab bwa_depth-pos.tab mode[contigs/complete] suffix")
    sys.exit()

name = os.path.basename(sys.argv[1]).split(sys.argv[4])[0]
spec_stats = {}

# recover genome length and total counts per species (for average depth)
with open(sys.argv[1], "r") as f:
    for line in f:
        if line[0] != "*": # exclude unmapped line
            cols = line.strip("\n").split("\t")
            if sys.argv[3] == "contigs":
                species = "_".join(cols[0].split("_")[:-1]) # species name
            elif sys.argv[3] == "complete":
                species = cols[0]
            clength = int(cols[1]) # contig length
            counts = int(cols[2]) # read counts
            if species not in spec_stats.keys(): # first time species is found
                spec_stats[species] = [clength, counts, []] # dict of length, counts and covPos
            else:
                spec_stats[species][0] += clength
                spec_stats[species][1] += counts

# recover covered positions (for calc mean depth, coverage and evenness)
with open(sys.argv[2], "r") as f:
    for line in f:
        if line[0] != "*": # exclude unmapped line
            cols = line.strip("\n").split("\t")
            if sys.argv[3] == "contigs":
                species = "_".join(cols[0].split("_")[:-1])
            elif sys.argv[3] == "complete":
                species = cols[0]
            depth = int(cols[2]) # depth of covered position
            spec_stats[species][2].append(depth)

print(f"Genome\t{name}_Length\t{name}_Counts\t{name}_MeanDepth\t{name}_Coverage\t{name}_ExpCoverage\t{name}_CoeffVar")

# combine stats and print per species
for species in spec_stats.keys():
    length = spec_stats[species][0]
    counts = spec_stats[species][1]
    covPos = spec_stats[species][2]
    covBases = len(covPos)
    coverage = float(covBases)/length*100
    if covBases < 1:
        covSd = 0
    else:
        covSd = numpy.std(covPos)
    meanDepth = float(sum(covPos))/max(length,1)
    if float(meanDepth) > 0:
        expCov = (1.00 - numpy.exp(-0.883*meanDepth))*100
        meanDepth_covered = float(sum(covPos))/max(len(covPos),1)
        cV = covSd/meanDepth_covered
    else:
        expCov = 0
        cV = 0
    print(f"{species}\t{length}\t{counts}\t{meanDepth:.2f}\t{coverage:.2f}\t{expCov:.2f}\t{cV:.2f}")