#!/usr/bin/env python3

import fileinput

with open(snakemake.output['fq'], 'wt') as f:
    for line in fileinput.input(snakemake.params['files']):
        f.write(line)
