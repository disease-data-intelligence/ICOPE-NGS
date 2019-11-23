#! /usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from datetime import datetime
import pdb
import sys
import subprocess
import os
from io import StringIO
# imports from own repo's
from utils_py.version import print_modules, imports
from computerome.somatic_setup import find_pairs
from quality.gene_coverage import calculate_statistics


def get_parser():
    parser = argparse.ArgumentParser(
        description="Run coverage stats for all")
    parser.add_argument('-bam', dest="bam", help="bam-file", type=str)
    parser.add_argument('-panel', '--gene-panel', dest='panel',
                        default='/home/projects/HT2_leukngs/data/references/general/315_genes_of_interest.txt',
                        help="(Default: /home/projects/HT2_leukngs/data/references/general/315_genes_of_interest.txt)")
    parser.add_argument('-bed', '--bedfile', type=str, dest='bed',
                        default="/home/projects/HT2_leukngs/data/references/hg37/test.bed",
                        help="Bed-file with intervals to look at. "
                             "(Default: /home/projects/HT2_leukngs/data/references/hg37/USCS.hg37.canonical.exons.bed)")
    return parser


def get_args(args=None):
    parser = get_parser()
    args = parser.parse_args(args)
    return args


def run_samtools(bam, bed):
    print("# Running samtools for collecting coverage stats")
    assert os.path.exists(bam), "does not exist"
    assert os.path.exists(bed), "does not exist"
    try:
        input = subprocess.Popen(["samtools", "bedcov", bed, bam], stdout=subprocess.PIPE)
    except OSError as e:
        print("# Could not run samtools, try to use: module load samtools/1.9 and/or check path")
    print("# Decoding bedcov output ... ")
    decoding = StringIO(input.communicate()[0].decode('utf-8'))
    data = pd.read_csv(decoding, sep='\t', header=None)
    data.columns = ['chromosome', 'start', 'end', 'gene', 'exon', 'strand', 'coverage']
    data.drop(columns=['strand'], inplace=True)
    return data


def calculate_coverage_stats(data, panel):
    print("# Calculating coverage stats ... ")
    data['mean_cov'] = data['coverage'] / (data['end'] - data['start'])
    data['exons_above20x_frac'] = data['mean_cov'].apply(lambda x: int(x > 20.0))
    coverage_chromosomes = data.groupby('chromosome').mean()
    # filter gene coverage data
    data_interest = data[data['gene'].isin(panel)]
    coverage_genes = data_interest.groupby('gene').mean()
    low_coverage_exons = data_interest[data_interest['mean_cov'] < 20.0]
    return coverage_genes.loc[:, ['mean_cov', 'exons_above20x_frac']], \
           low_coverage_exons.loc[:, ['chromosome', 'gene', 'exon', 'mean_cov']], \
           coverage_chromosomes.loc[:, ['mean_cov', 'exons_above20x_frac']]


def main(bam, bed, panel_file):
    panel = list(pd.read_csv(panel_file).values.flatten())
    print(f"# Processing {bam} ... ")
    coverage_genes, low_cov_exons, coverage_chrom = calculate_coverage_stats(run_samtools(bam, bed), panel)
    sample = os.path.basename(bam).replace('.bam', '')
    path = os.path.dirname(bam)
    output = path + sample
    print(f"# Writing output to in {output}")
    coverage_genes.to_csv(sample + '_genes.tsv', sep='\t')
    coverage_chrom.to_csv(sample + '_chromosomes.tsv', sep='\t')
    low_cov_exons.to_csv(sample + '_low_cov_exons.tsv', sep='\t', index=False)


if __name__ == "__main__":
    start_time = datetime.now()
    args = get_args()
    print("# args:", args)
    global_modules = globals()
    modules = imports(global_modules)
    print_modules(list(modules))
    main(args.bam, args.bed, args.panel)
    end_time = datetime.now()
    print("# Done!")
    print('# Duration: {}'.format(end_time - start_time))
