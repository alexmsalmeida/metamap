# snakemake workflow for mapping metagenomes against reference database

import os

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
        OUTPUT_DIR+"/summary/bwa_cov-est.csv",
        OUTPUT_DIR+"/summary/bwa_cov-exp.csv",
        OUTPUT_DIR+"/summary/bwa_counts-total.csv",
        OUTPUT_DIR+"/summary/bwa_counts-unique.csv"

# map metagenomes
rule map2ref:
    input:
        fwd = lambda wildcards: samp2path[wildcards.sample][0],
        rev = lambda wildcards: samp2path[wildcards.sample][1]
    output:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_total.tab",
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REFNAME+"_unique.tab"
    params:
        outpref = OUTPUT_DIR+"/mapping/{sample}/{sample}",
        bwa_db = config["database"]["bwa_db"],
        reftype = config["database"]["type"]
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/map2ref.sh -t 8 -i {input.fwd} -n {input.rev} -r {params.bwa_db} -o {params.outpref} -c {params.reftype}
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
