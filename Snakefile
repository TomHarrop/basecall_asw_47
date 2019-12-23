#!/usr/bin/env python3

import fileinput
import multiprocessing
from pathlib import Path

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

rule target:
    input:
        # expand(Path(outdir, 'merged', '{group}.fq.gz').as_posix(),
        #        group=['pool', 'asw47']),
        expand(Path(outdir,
                    'minionqc',
                    '{group}').as_posix(),
               group=['pool', 'asw47'])

rule minionqc:
    input:
        Path(tempdir, 'minionqc' '{group}')
    output:
        directory(Path(outdir, 'minionqc', '{group}'))
    threads:
        min(len(fc_list), multiprocessing.cpu_count())
    singularity:
        minionqc_container
    log:
        Path(logdir, 'minionqc.{group}.log')
    shell:
        'sleep 100 ; '
        'ls -alh . ; '
        'pwd ; '
        'MinIONQC.R '
        '--processors={threads} '
        '--input={input} '
        '--outputdirectory={output} '
        # '&> {log}'


rule manual_shadow:
    input:
        group_bc_output
    output:
        temp(directory(Path(tempdir, 'minionqc' '{group}')))
    singularity:
        minionqc_container
    shell:
        'mkdir -p {output} ; '
        'ln -s {input} {output}/ ; '
        'ls -lh {output}'

rule compress:
    input:
        Path(tempdir, '{group}', 'all_pass.fastq')
    output:
        Path(outdir, 'merged', '{group}.fq.gz')
    singularity:
        pigz_container
    threads:
        multiprocessing.cpu_count()
    shell:
        'pigz --processes {threads} --best --keep {input} '
        '; '
        'mv {input}.gz {output}'


rule combine:
    input:
        group_bc_output
    output:
        fq = temp(Path(tempdir, '{group}', 'all_pass.fastq'))
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


