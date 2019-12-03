#!/usr/bin/env python3

import fileinput
import multiprocessing
import os

def find_basecalled_fastq_files():
    glob_results = snakemake.io.glob_wildcards(
        'basecall_guppy2.3.7/basecalled/{fc}/pass/{fq}.fastq')
    all_files = snakemake.io.expand(
        'basecall_guppy2.3.7/basecalled/{fc}/pass/{fq}.fastq',
        zip,
        fc=glob_results.fc,
        fq=glob_results.fq)
    return [x for x in all_files if os.path.isfile(x)]


fc_list = ['20181102_0046_ASW47-t7-3',
    '20181108_0252_ASW47-B2',
    '20181102_2243_ASW47-t7-3',
    '20181109_0257_ASW47-C1',
    '20181104_2357_ASW47-t7-124',
    '20181112_0016_ASW47-C2',
    '20181107_0140_ASW47-B1']

pigz_container = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
guppy_container = 'shub://TomHarrop/singularity-containers:guppy_2.3.7'
minionqc_container = 'shub://TomHarrop/singularity-containers:minionqc_1.4.1'


rule target:
    input:
        'basecall_guppy2.3.7/minionqc/combinedQC/summary.yaml',
        'basecall_guppy2.3.7/merged/all_pass.fastq.gz'


rule minionqc:
    input:
        expand('basecall_guppy2.3.7/basecalled/{fc}/sequencing_summary.txt',
               fc=fc_list)
    output:
        'basecall_guppy2.3.7/minionqc/combinedQC/summary.yaml'
    params:
        search_dir = 'basecall_guppy2.3.7/basecalled',
        outdir = 'basecall_guppy2.3.7/minionqc'
    threads:
        min(len(fc_list), multiprocessing.cpu_count())
    priority:
        1
    singularity:
        minionqc_container
    log:
        'basecall_guppy2.3.7/logs/minionqc.log'
    shell:
        'MinIONQC.R '
        '--processors={threads} '
        '--input={params.search_dir} '
        '--outputdirectory={params.outdir} '
        '&> {log}'


rule compress:
    input:
        'basecall_guppy2.3.7/merged/all_pass.fastq'
    output:
        'basecall_guppy2.3.7/merged/all_pass.fastq.gz'
    singularity:
        pigz_container
    threads:
        multiprocessing.cpu_count()
    shell:
        'pigz --processes {threads} --best --keep {input}'


rule combine:
    input:
        expand('basecall_guppy2.3.7/basecalled/{fc}/sequencing_summary.txt',
               fc=fc_list)
    output:
        fq = temp('basecall_guppy2.3.7/merged/all_pass.fastq')
    params:
        files = lambda wildcards: find_basecalled_fastq_files()
    run:
        with open(output.fq, 'wt') as f:
            for line in fileinput.input(params.files):
                f.write(line)


rule basecall:
    input:
        'basecall_guppy2.3.7/temp/{fc}'
    output:
        'basecall_guppy2.3.7/basecalled/{fc}/sequencing_summary.txt'
    params:
        outdir = 'basecall_guppy2.3.7/basecalled/{fc}'
    log:
        'basecall_guppy2.3.7/logs/{fc}_guppy.log'
    resources:
        gpu = 1
    singularity:
        guppy_container
    shell:
        'guppy_basecaller '
        '--flowcell FLO-MIN106 '
        '--kit SQK-LSK109 '
        '--input {input} '
        '--save_path {params.outdir} '
        '--qscore_filtering true '
        '--enable_trimming true '
        '--recursive '
        '--device auto '
        '--num_callers 16 '
        '--chunks_per_runner 96 '
        '&> {log}'


rule untar:
    input:
        'data/{fc}.tar.gz'
    output:
        temp(directory('basecall_guppy2.3.7/temp/{fc}'))
    log:
        'basecall_guppy2.3.7/logs/{fc}_untar.log'
    singularity:
        pigz_container
    shell:
        'mkdir -p {output} ; '
        'tar -zx -f {input} -C {output} --strip-components 1 &> {log}'


