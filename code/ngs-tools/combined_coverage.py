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
import matplotlib.pyplot as plt
import seaborn as sns


# imports from own repo's
from utils_py.version import print_modules, imports
from computerome.somatic_setup import find_pairs
from quality.gene_coverage import calculate_statistics

# Set up command line parser
def get_parser():
    parser = argparse.ArgumentParser(
        description="Run coverage stats for all")
    parser.add_argument('-in', dest="infile", help="Coverage file from samtools or bam-file", type=str)
    parser.add_argument('-panel', '--gene-panel', dest='panel',
                        default='/home/projects/HT2_leukngs/data/references/general/315_genes_of_interest.txt',
                        help="(Default: /home/projects/HT2_leukngs/data/references/general/315_genes_of_interest.txt)")
    parser.add_argument('-bed', '--bedfile', type=str, dest='bed',
                        default="/home/projects/HT2_leukngs/data/references/hg37/test.bed",
                        help="Bed-file with intervals to look at. "
                             "(Default: /home/projects/HT2_leukngs/data/references/hg37/USCS.hg37.canonical.exons.bed)")
    parser.add_argument('-out', dest='out')
    return parser


def get_args(args=None):
    parser = get_parser()
    args = parser.parse_args(args)
    return args


def run_samtools(bam, bed):
    print("# Running samtools for collecting coverage stats")
    # Bedcov reports the total read base count (i.e. the sum of per base read depths)
    # for each genomic region specified in the supplied BED file.
    assert os.path.exists(bam), "does not exist"
    assert os.path.exists(bed), "does not exist"
    try:
        input = subprocess.Popen(["samtools", "bedcov", bed, bam], stdout=subprocess.PIPE)
    except OSError as e:
        print("# Could not run samtools, try to use: module load samtools/1.9 and/or check path")
    print("# Decoding bedcov output ... ")
    decoding = StringIO(input.communicate()[0].decode('utf-8'))
    data = pd.read_csv(decoding, dtype={0: str}, sep='\t', header=None)
    data.columns = ['chromosome', 'start', 'end', 'gene', 'exon', 'strand', 'coverage']
    data.drop(columns=['strand'], inplace=True)
    return data


def load_data(infile):
    print(f"# Loading data from {infile}")
    data = pd.read_csv(infile, dtype={0: str}, sep='\t', header=None)
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
    print("# Number of unique genes found:", len(data['gene'].unique()))
    print("# Found", len(data_interest), "out of", len(panel))
    return coverage_genes.loc[:, ['mean_cov', 'exons_above20x_frac']], \
           low_coverage_exons.loc[:, ['chromosome', 'gene', 'exon', 'mean_cov']], \
           coverage_chromosomes.loc[:, ['mean_cov', 'exons_above20x_frac']]


def plot_distribution(genes, name, plots=4, sort='value'):
    if sort == 'value':
        data = genes.sort_values(by='mean_cov')
    elif sort == 'alphabet':
        data = genes.sort_values(by='gene')
    else:
        data = genes
    rows = int(np.ceil(len(data)/plots))
    fig, axes = plt.subplots(plots, 1, figsize=(10, plots*3))
    cmap = sns.color_palette("GnBu")
    for i, ax in enumerate(fig.axes):
        sns.barplot(x='gene', ax=ax, y='mean_cov', palette=cmap, data=data[i*rows:(i+1)*rows])
        ax.set_title('Coverage pr. target gene')
        ax.set_ylabel('Mean coverage')
        ax.set_xlabel('Gene name')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, verticalalignment='top', fontweight='light')
    fig.tight_layout()
    plotname = name + '_coverage_pr_gene.png'
    print("# Saving to", plotname)
    plt.savefig(plotname)


def write_excel(output, coverage_genes, coverage_chrom, low_cov_exons):
    print("# Writing tsv files ... ")
    coverage_genes.to_csv(output + '_genes.tsv', sep='\t')
    coverage_chrom.to_csv(output + '_chromosomes.tsv', sep='\t')
    low_cov_exons.to_csv(output + '_low_cov_exons.tsv', sep='\t', index=False)

    with pd.ExcelWriter(output + 'xlsx') as excel_obj:
        print(f"# Writing xlsx file to {output}")
        coverage_genes.to_excel(excel_obj, index_label='Gene coverage', sheet_name='Gene coverage')
        coverage_chrom.to_excel(excel_obj, index_label='Chromosome', sheet_name='Chromosome coverage in exonic regions')
        low_cov_exons.to_excel(excel_obj, sheet_name='Low coverage exons', index=False)


def main(infile, bed, panel_file, file_suffix=None):
    gene_panel = list(pd.read_csv(panel_file).values.flatten())

    if infile.endswith('.bam'):
        suffix = '.bam'
        coverage_genes, low_cov_exons, coverage_chrom = calculate_coverage_stats(run_samtools(infile, bed), gene_panel)
    elif infile.endswith('.bed'):
        suffix = '.bed'
        coverage_genes, low_cov_exons, coverage_chrom = calculate_coverage_stats(load_data(infile), gene_panel)
    else:
        print("# File format not recognized, exiting ... ")
        exit(1)

    if file_suffix:
        sample = os.path.basename(infile).replace(suffix, '') + '_' + file_suffix
    else:
        sample = os.path.basename(infile).replace(suffix, '')

    path = os.path.dirname(infile)
    output = path + sample
    print(f"# Writing output to in {output}")
    write_excel(output, path, sample, coverage_genes, coverage_chrom, low_cov_exons)
    plot_distribution(coverage_genes, output)


if __name__ == "__main__":
    start_time = datetime.now()
    args = get_args()
    print("# args:", args)
    global_modules = globals()
    modules = imports(global_modules)
    print_modules(list(modules))
    main(args.infile, args.bed, args.panel, args.out)
    end_time = datetime.now()
    print("# Done!")
    print('# Duration: {}'.format(end_time - start_time))
