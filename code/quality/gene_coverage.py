#! /usr/bin/env python3

import seaborn as sns
import pandas as pd
import sys
import re
import pdb
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
sys.path.append("/home/projects/HT2_leukngs/apps/github/code/utilities")
import version

def read_input(filename):
    data = pd.read_csv(filename, sep='\t')
    data.columns = ['chromosome', 'start', 'end', 'gene', 'exon', 'strand', 'coverage']
    data.drop(columns=['strand'], inplace=True)
    return data


def calculate_statistics(data, panel):
    data['mean_cov'] = data['coverage'] / (data['end'] - data['start'])
    genes = data.groupby('gene')['mean_cov'].mean().to_frame().reset_index()
    genes = genes[genes['gene'].isin(panel)].reset_index()
    return genes


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
    plotname = name + "_" + sort[0] + "sort_" + 'coverage_pr_gene.png'
    print("# Saving to", plotname)
    plt.savefig(plotname)


if __name__ == '__main__':
    print("# Running gene coverage function")
    global_modules = globals()
    modules = version.imports(global_modules)
    version.print_modules(list(modules))
    filename = sys.argv[1]
    gene_file = sys.argv[2]
    gene_panel = list(pd.read_csv(gene_file).values.flatten())
    print("# plotting, ", filename)
    data = read_input(filename)
    genes = calculate_statistics(data, panel=gene_panel)
    print("# Number of unique genes plotted:", len(genes['gene'].unique()))
    samplename = filename.replace('.bed', '')
    for sorting in ['value', 'alphabet']:
        plot_distribution(genes, name=samplename, plots=4, sort=sorting)



