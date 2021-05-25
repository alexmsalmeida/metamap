#!/usr/bin/env python

import glob
import sys
import os
import argparse

def getCounts(folder, ftype):
    counts = []
    file_count = 0
    for f in glob.glob(os.path.join(folder, "*/*"+ftype+".tab")):
        file_count += 1
        print("%i\t%s" % (file_count, f))
        with open(f, "r") as bwa_data:
            linen = 0
            for line in bwa_data:
                linen += 1
                cols = line.strip("\n").split("\t")
                if file_count == 1:
                    counts.append([cols[0]])
                if linen == 1:
                    name = cols[2].replace("_Counts", "")
                    counts[linen-1].append(name)
                else:
                    count = int(cols[2])
                    counts[linen-1].append(str(count))
    return counts

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse BWA results to extract counts')
    parser.add_argument('-i', dest='in_folder', help='Input folder (files should be in subdirectories)', required=True)
    parser.add_argument('-f', dest='ftype', help='\'unique\' or \'total\' counts', required=True)
    parser.add_argument('-o', dest='out_file', help='Output CSV file', required=True)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
        counts = getCounts(args.in_folder, args.ftype)
        with open(args.out_file, "w") as f_out:
            for data in counts:
                f_out.write(",".join(data)+"\n")
