#!/usr/bin/env snakemake

import oncopipe as op
import os
import pandas

##### SNAKEMAKE COMMAND
# snakemake -p --use-conda --scheduler greedy --jobs 5000 --latency-wait 120 --keep-going --default-resources mem_mb=380000 disk_mb=2000  --cluster-sync "srun -p upgrade -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads}" -s last.dnarrange.smk all -np
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
        fa = "04-consensus-seq/maf/{sample_id}.merged.fa"
    params:
        train = "01-20985T.train"
    conda: "/projects/rmorin_scratch/ONT_scratch/results/last/dnarrange.yaml"
    shell:
        op.as_one_line("""
            dnarrange-merge {input.fq} \
            {params.train} {input.maf} > {output.fa}
        """)


## realign merged dnarrange to genome with new DB
rule last_realign:
    input:
        fa = str(rules.dnarrange_merge.output.fa)
    output:
        maf = "05-realign-consensus/maf/{sample_id}.merged.maf"
    params:
        opts = "--split -p",
        threads = "24",
        sensitivity = "100",
        train = "/projects/rmorin_scratch/ONT_scratch/results/ref/01-20985.merged.par",
        ref_db = directory("/projects/rmorin_scratch/ONT_scratch/results/ref/last-db_NEAR/mydb")
    conda: "/projects/rmorin_scratch/ONT_scratch/results/last/dnarrange.yaml"
    shell:
        op.as_one_line("""
            lastal -P{params.threads} -m{params.sensitivity} {params.opts} {params.train} \
            {params.ref_db} \
            {input.fa} > {output.maf}
        """)

## draw pics
#last-multiplot 01-20985T.final.maf 01-20985T-pics
rule last_multiplot:
    input:
        maf = str(rules.last_realign.output.maf)
    output:
        pic_dir = directory("06-realign-pics/{sample_id}")
    conda: "/projects/rmorin_scratch/ONT_scratch/results/last/dnarrange.yaml"
    shell:
        op.as_one_line("""
            last-multiplot {input.maf} {output.pic_dir}
        """)

## edit maf to remove mismap rate
rule remove_mismap:
    input:
        maf = str(rules.last_realign.output.maf)
    output:
        maf = "05-realign-consensus/maf/{sample_id}.no_mismap.maf"
    shell:
        op.as_one_line("""
            sed '1 i ##maf version=1' {input.maf} | sed -r 's/mismap=[0-9\.e-]*//' > {output.maf}
        """)

## parse maf for coordinates
rule parse_maf:
    input:
        alignment_maf = str(rules.remove_mismap.output.maf)
    output:
        alignment_df = "07-tsv-coordinates/tsv/{sample_id}.tsv"
    conda: "/projects/rmorin_scratch/ONT_scratch/results/last/parse_maf.yaml"
    script:
        "/projects/rmorin_scratch/ONT_scratch/results/last/parse_maf.py"

rule all:
    input:
        expand(rules.last_alignment.output.maf,
        zip,
        sample_id=SAMPLES['sample_id']),
        expand(rules.last_multiplot.output.pic_dir,
        zip,
        sample_id=MCLS['sample_id']),
        expand(rules.parse_maf.output.alignment_df,
        zip,
        sample_id=MCLS['sample_id'])
