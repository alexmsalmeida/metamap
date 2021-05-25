# snakemake workflow for mapping metagenomes against reference database

import os
import glob

# preparing files
configfile: "config.yml"

INPUT_FILE = config["input"]
OUTPUT_DIR = config["output"]
REFNAME = os.path.basename(config["database"]["bwa_db"]).split(".fa")[0]

os.system("chmod -R +x tools")

samp2path = {}
with open(INPUT_FILE) as f:
    for line in f:
        cols = line.strip().split()
        fwd = cols[0]
        rev = cols[1]
        sample = os.path.basename(cols[0]).split("_1.fastq")[0]
        samp2path[sample] = [fwd, rev]

for sample in samp2path.keys():
    if not os.path.exists(OUTPUT_DIR+"/mapping/"+sample+"/logs"):
        os.makedirs(OUTPUT_DIR+"/mapping/"+sample+"/logs")

if not os.path.exists(OUTPUT_DIR+"/summary/logs"):
    os.makedirs(OUTPUT_DIR+"/summary/logs")

# rule that specifies the final expected output files
rule all:
    input:
        expand(OUTPUT_DIR+"/mapping/{sample}/done.txt", sample=samp2path.keys()),
        OUTPUT_DIR+"/summary/bwa_cov-est.csv",
        OUTPUT_DIR+"/summary/bwa_cov-exp.csv",
        OUTPUT_DIR+"/summary/bwa_counts-total.csv",
        OUTPUT_DIR+"/summary/bwa_counts-unique.csv"

# map metagenomes
rule bwa_mem:
    input:
        fwd = lambda wildcards: samp2path[wildcards.sample][0],
        rev = lambda wildcards: samp2path[wildcards.sample][1]
    output:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_raw.bam"
    params:
        bwa_db = config["database"]["bwa_db"],
    conda:
        "envs/metamap.yml"
    shell:
        """
        bwa mem -M -t 8 {params.bwa_db} {input.fwd} {input.rev} | samtools view -@ 7 -F 256 -uS - > {output}
        """

# parse with samtools
rule samtools_sort:
    input:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_raw.bam"
    output:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_raw_sorted.bam"
    params:
        tmp = lambda wildcards: glob.glob(os.path.join(OUTPUT_DIR, "mapping", wildcards.sample, "*bam.tmp*"))
    conda:
        "envs/metamap.yml"
    shell:
        """
        rm -rf {params.tmp}
        samtools sort -@ 7 {input} -o {output}
        """

rule samtools_filter_alns:
    input:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_raw_sorted.bam"
    output:
        idx = OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_raw_sorted.bam.bai",
        bam = OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_total_sorted.bam",
        bai = OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_total_sorted.bam.bai"
    conda:
        "envs/metamap.yml"
    shell:
        """
        samtools index -@ 7 {input}
        tools/bam_ani-filter.py {input} 90 60 | samtools view -u - -o {output.bam}
        samtools index -@ 7 {output.bam}
        """

rule samtools_extract_unique:
    input:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_total_sorted.bam"
    output:
        bam = OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_unique_sorted.bam",
        bai = OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_unique_sorted.bam.bai"
    conda:
        "envs/metamap.yml"
    shell:
        """
        samtools view -@ 7 -q 1 -f 2 -u {input} -o {output.bam}
        samtools index -@ 7 {output.bam}
        """

rule samtools_idx:
    input:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_{type}_sorted.bam"
    output:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_{type}_depth.tab"
    conda:
        "envs/metamap.yml"    
    shell:
        "samtools idxstats {input} > {output}"

rule samtools_depth:
    input:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_{type}_sorted.bam"
    output:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_{type}_depth-pos.tab"
    conda:
        "envs/metamap.yml"    
    shell:
        "samtools depth {input} > {output}"

rule samtools_combine:
    input:
        depth = OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_{type}_depth.tab",
        pos = OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_{type}_depth-pos.tab"
    output:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_{type}.tab"
    params:
        reftype = config["database"]["type"],
        suffix = "_"+REFNAME,
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/parse_bwa-depth.py {input.depth} {input.pos} {params.reftype} {params.suffix} > {output}
        """

# parse final output
rule parse_cov:
    input:
        expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_total.tab", sample=samp2path.keys())
    output:
        OUTPUT_DIR+"/summary/bwa_cov-est.csv"
    params:
        outdir = OUTPUT_DIR+"/mapping"
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/parse_bwa-cov.py -i {params.outdir} -o {output}
        """

rule parse_expcov:
    input:
        expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_total.tab", sample=samp2path.keys())
    output:
        OUTPUT_DIR+"/summary/bwa_cov-exp.csv"
    params:
        outdir = OUTPUT_DIR+"/mapping"
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/parse_bwa-expcov.py -i {params.outdir} -o {output}
        """

rule parse_counts:
    input:
        lambda wildcards: expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_{type}.tab", sample=samp2path.keys(), type=wildcards.type)
    output:
        OUTPUT_DIR+"/summary/bwa_counts-{type}.csv"
    params:
        outdir = OUTPUT_DIR+"/mapping"
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/parse_bwa-counts.py -i {params.outdir} -f {wildcards.type} -o {output}
        """

# clean tmp
rule clean_tmp:
    input:
        lambda wildcards: expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_{type}.tab", sample=wildcards.sample, type=["unique", "total"])
    output:
        touch(OUTPUT_DIR+"/mapping/{sample}/done.txt")
    params:
        raw = lambda wildcards: glob.glob(os.path.join(OUTPUT_DIR, "mapping", wildcards.sample, "*raw*")),
        sort = lambda wildcards: glob.glob(os.path.join(OUTPUT_DIR, "mapping", wildcards.sample, "*sorted*")),
        depth = lambda wildcards: glob.glob(os.path.join(OUTPUT_DIR, "mapping", wildcards.sample, "*depth*"))
    shell:
        "rm -rf {params.raw} {params.sort} {params.depth}"
