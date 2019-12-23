#!/usr/bin/env python3

import fileinput

print(snakemake.params)
print(snakemake.output)

quit(1)

with open(snakemake.output['fq'], 'wt') as f:
    for line in fileinput.input(snakemake.params['files']):
        f.write(line)
