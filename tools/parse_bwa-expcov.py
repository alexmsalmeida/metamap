#!/usr/bin/env python

import glob
import sys
import os
import argparse

def getCov(folder):
    covs = []
    file_count = 0
    for f in glob.glob(os.path.join(folder, "*/*total.tab")):
        file_count += 1
        print("%i\t%s" % (file_count, f))
        with open(f, "r") as bwa_data:
            linen = 0
            for line in bwa_data:
                linen += 1
                cols = line.strip("\n").split("\t")
                if file_count == 1:
                    covs.append([cols[0]])
                if linen == 1:
                    name = cols[5].replace("_ExpCoverage", "")
                    covs[linen-1].append(name)
                else:
                    cov = float(cols[5])
                    covs[linen-1].append(str(cov))
    return covs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse BWA results to extract coverage')
    parser.add_argument('-i', dest='in_folder', help='Input folder (files should be in subdirectories)', required=True)
    parser.add_argument('-o', dest='out_file', help='Output CSV file', required=True)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
        covs = getCov(args.in_folder)
        with open(args.out_file, "w") as f_out:
            for data in covs:
                f_out.write(",".join(data)+"\n")            
