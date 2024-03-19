#!/usr/bin/env snakemake

import oncopipe as op
import os
import pandas

##### SNAKEMAKE COMMAND
# snakemake -np --use-conda -s mod_mitohpc_snakemake.smk all --cores 1 --config unix_group=all_the_things
## had to increase mem_mb very high. double check this.
# snakemake -p --use-conda --scheduler greedy --jobs 5000 --latency-wait 120 --keep-going --default-resources mem_mb=2000 disk_mb=2000  --cluster-sync "srun -p upgrade -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads}" -s last.dnarrange.smk all -np

##### SETUP VARIABLES

##### CONFIG FILES

##### SETUP SAMPLES
SAMPLES = pandas.read_csv('/projects/rmorin_scratch/ONT_scratch/samples.last.tsv', sep = "\t")
CONTROLS = op.filter_samples(SAMPLES, pathology = ["BL", "DLBCL"])
MCLS = op.filter_samples(SAMPLES, pathology = "MCL")

#print(SAMPLES)
rule symlink_fastq:
    input:
        fq = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/data/promethION_fastq/genome/{sample_id}.fastq.gz"
    output:
        fq = "00-inputs/fastq/{sample_id}.fq.gz"
    run:
        op.absolute_symlink(input.fq, output.fq)

rule last_alignment:
    input:
        fq = str(rules.symlink_fastq.output.fq),
    output:
        maf = "01-last/maf/{sample_id}.maf.gz"
    params:
        train = "/projects/rmorin_scratch/ONT_scratch/results/last/01-20985T.train",
        ref_db = directory("/projects/rmorin_scratch/ONT_scratch/results/ref/last-db/mydb")
    conda: "/projects/rmorin_scratch/ONT_scratch/results/last/last.yaml"
    shell:
        op.as_one_line("""
            lastal -P24 --split -p {params.train} {params.ref_db} {input.fq} | gzip > {output.maf}
        """)

rule dnarrange:
    input:
        control = expand(str(rules.last_alignment.output.maf), 
        zip,
        sample_id=CONTROLS['sample_id']),
        maf = str(rules.last_alignment.output.maf)
    output:
        maf = "02-dnarrange/maf/{sample_id}.maf"
    conda: "/projects/rmorin_scratch/ONT_scratch/results/last/dnarrange.yaml"
    shell:
        op.as_one_line("""
            dnarrange {input.maf} : {input.control} > {output.maf}
        """)

#dnarrange stricter can run on existing maf
rule dnarrange_strict:
    input:
        maf = str(rules.dnarrange.output.maf)
    output:
        maf = "02-dnarrange/maf/{sample_id}.strict.maf"
    conda: "/projects/rmorin_scratch/ONT_scratch/results/last/dnarrange.yaml"
    shell:
        op.as_one_line("""
            dnarrange -s3 {input.maf} > {output.maf}
        """)

## merge dnarrange into consensus seq
rule dnarrange_merge:
    input:
        fq = str(rules.symlink_fastq.output.fq),
        maf = str(rules.dnarrange_strict.output.maf)
    output:
        maf = "04-consensus-seq/maf/{sample_id}.merged.fa"
    params:
        train = "01-20985T.train"
    conda: "/projects/rmorin_scratch/ONT_scratch/results/last/dnarrange.yaml"
    shell:
        op.as_one_line("""
            dnarrange-merge {input.fq} \
            {params.train} {input.maf} > {output.maf}
        """)


## realign merged dnarrange to genome with new DB
#lastal -P24 -m20 --split -p ../../ref/01-20985.merged.par ../../ref/last-db_NEAR/mydb ../04-consensus-seq/01-20985T.merged.fa > 01-20985T.merged.maf

## draw pics
#last-multiplot 01-20985T.final.maf 01-20985T-pics

rule all:
    input:
        expand(rules.last_alignment.output.maf,
        zip,
        sample_id=SAMPLES['sample_id']),
        expand(rules.dnarrange.output.maf,
        zip,
        sample_id=MCLS['sample_id'])
