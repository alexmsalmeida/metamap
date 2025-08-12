#!/usr/bin/env python3

import argparse

def load_name_mapping(name_list_file):
    """Load the MAGs to Contigs mapping from name list file"""
    mag_to_contigs = {}
    with open(name_list_file, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line:  # Skip empty lines
                parts = line.split(maxsplit=1)
                if len(parts) == 2:
                    mag, contig = parts
                    mag_to_contigs[mag] = contig
                else:
                    print(f"Warning: Invalid line format in {name_list_file} (line {line_num}): {line}")
    return mag_to_contigs

def ren_fasta(args):
    """Rename sequences in a FASTA file by replacing MAG substrings with Contigs values and removing '_k141'"""
    mag_to_contigs = load_name_mapping(args.name_list_file)
    
    with open(args.fasta_file, "r") as f, open(args.output_file, "w") as out:
        for line in f:
            if line.startswith(">"):
                original_header = line.strip()
                header_content = original_header[1:]  
                new_header = header_content.split('_k141')[0]

                modified_content = new_header
                for mag, contig in mag_to_contigs.items():
                    if mag in modified_content:  
                        modified_content = modified_content.replace(mag, contig)
                        break  
                out.write(f">{modified_content}\n")
            else:
                out.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rename FASTA headers by replacing MAG substrings with Contig values and removing "_k141"')
    parser.add_argument('-f', dest='fasta_file', help='Input FASTA file', required=True)
    parser.add_argument('-o', dest='output_file', help='Output file', default='rename.fa')
    parser.add_argument('-n', dest='name_list_file', help='Path to file containing MAGs and Contigs mappings', required=True)

    args = parser.parse_args()
    ren_fasta(args)
