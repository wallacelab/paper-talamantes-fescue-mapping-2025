
# This file is a snake make tutorial.
# you must run it with a prefix to tell the pipline how to name files that have {prefix}

# This is how to run this file 
# snakemake --snakefile Rules_Practice.smk --cores 1
import pandas as pd

rule Hello:
    conda: 
        "/home/drt06/applications/miniconda3/envs/Alignment"
    output:
        '{prefix}.txt'  
    shell:
        'echo Hello world > {output}'

#reads output from last time and adds something to it
rule world:
    input: 
        '{prefix}.txt'
    output:
        '{prefix}.second.csv'
    shell:
        '''
        cat {input} > {output}
        echo This is the second line >> {output}
        '''

rule bello:
    input:
        '{prefix}.txt'
    output:
        '{prefix}.bello.tsv'
    run:
        data = pd.read_csv(input[0], header=None)
        data.iloc[0, 0] = 1 
        data.to_csv(output[0], sep = '\t')