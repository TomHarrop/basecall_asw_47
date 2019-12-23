#!/usr/bin/env python3

import fileinput
import multiprocessing
from pathlib import Path


# GLOBALS
outdir = 'basecall_guppy341'
tempdir = Path(outdir, 'temp')
logdir = Path(outdir, 'logs')

fc_list = [
    '20181102_0046_ASW47-t7-3',
    '20181102_2243_ASW47-t7-3',
    '20181104_2357_ASW47-t7-124',
    '20181107_0140_ASW47-B1',
    '20181108_0252_ASW47-B2',
    '20181109_0257_ASW47-C1',
    '20181112_0016_ASW47-C2',
    'sample_10',
    'sample_3',
    'sample_5',
    'sample_7',
    'sample_8']

pool_list = [x for x in fc_list
             if x.startswith('sample')]

asw47_list = [x for x in fc_list
              if x.startswith('201811')]


pigz_container = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
guppy_container = 'shub://TomHarrop/ont-containers:guppy_3.4.1'
minionqc_container = 'shub://TomHarrop/singularity-containers:minionqc_1.4.1'


# rule target:
#     input:
#         'basecall_guppy2.3.7/minionqc/combinedQC/summary.yaml',
#         'basecall_guppy2.3.7/merged/all_pass.fastq.gz'


# rule minionqc:
#     input:
#         expand('basecall_guppy2.3.7/basecalled/{fc}/sequencing_summary.txt',
#                fc=fc_list)
#     output:
#         'basecall_guppy2.3.7/minionqc/combinedQC/summary.yaml'
#     params:
#         search_dir = 'basecall_guppy2.3.7/basecalled',
#         outdir = 'basecall_guppy2.3.7/minionqc'
#     threads:
#         min(len(fc_list), multiprocessing.cpu_count())
#     priority:
#         1
#     singularity:
#         minionqc_container
#     log:
#         'basecall_guppy2.3.7/logs/minionqc.log'
#     shell:Path(outdir, 'basecalled', '{fc}', 'sequencing_summary.txt')
#         'MinIONQC.R '
#         '--processors={threads} '
#         '--input={params.search_dir} '
#         '--outputdirectory={params.outdir} '
#         '&> {log}'


# rule compress:
#     input:
#         'basecall_guppy2.3.7/merged/all_pass.fastq'
#     output:
#         'basecall_guppy2.3.7/merged/all_pass.fastq.gz'
#     singularity:
#         pigz_container
#     threads:
#         multiprocessing.cpu_count()
#     shell:
#         'pigz --processes {threads} --best --keep {input}'


rule target:
    input:
        expand(Path(tempdir, '{group}', 'all_pass.fastq').as_posix(),
               group=['pool', 'asw47'])


def group_bc_output(wildcards):
    my_fc_list = pool_list if wildcards.group == 'pool' else asw47_list
    return expand(Path(outdir,
                       'basecalled',
                       '{fc}',
                       'sequencing_summary.txt').as_posix(),
                  fc=my_fc_list)


def find_basecalled_fastq_files(wildcards):
    my_fc_list = pool_list if wildcards.group == 'pool' else asw47_list
    my_files = []
    for fc in my_fc_list:
        my_path = Path(outdir, 'basecalled', fc, 'pass', '{fq}.fastq')
        glob_results = glob_wildcards(my_path)
        all_files = expand(
            my_path,
            fq=glob_results.fq)
        my_files.append(x for x in all_files if Path(x).is_file())
    return my_files


rule combine:
    input:
        group_bc_output
    output:
        fq = Path(tempdir, '{group}', 'all_pass.fastq')
    params:
        files = lambda wildcards: find_basecalled_fastq_files()
    run:
        with open(output.fq, 'wt') as f:
            for line in fileinput.input(params.files):
                f.write(line)


rule basecall:
    input:
        Path(tempdir, '{fc}')
    output:
        Path(outdir, 'basecalled', '{fc}', 'sequencing_summary.txt')
    params:
        outdir = lambda wildcards:
            Path(outdir, 'basecalled', wildcards.fc)
    log:
        Path(logdir, 'basecall.{fc}.log')
    resources:
        gpu = 1
    singularity:
        guppy_container
    shell:
        'guppy_basecaller '
        '--flowcell FLO-MIN106 '
        '--kit SQK-LSK109 '
        '-i {input} '
        '-s {params.outdir} '
        '--qscore_filtering true '
        '--recursive '
        '--device auto '
        '&> {log}'

rule untar:
    input:
        'data/{fc}.tar.gz'
    output:
        temp(directory(Path(tempdir, '{fc}')))
    log:
        Path(logdir, 'untar.{fc}.log')
    singularity:
        pigz_container
    shell:
        'mkdir -p {output} ; '
        'tar -zx -f {input} -C {output} --strip-components 1 &> {log}'


