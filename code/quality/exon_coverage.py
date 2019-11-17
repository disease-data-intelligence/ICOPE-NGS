#! /usr/bin/env python3

import seaborn as sns
import pandas as pd
import sys
import re
import matplotlib.pyplot as plt
import pdb
import numpy as np
sys.path.append("/home/projects/HT2_leukngs/apps/github/code/utilities")
import version 

def read_input(filename):
    data = pd.read_csv(filename, sep='\t')
    data.columns = ['chromosome', 'start', 'end', 'gene', 'exon', 'strand', 'coverage']
    data.drop(columns=['strand'], inplace=True)

    return data


def calculate_statistics(data):
    data['mean_cov'] = data['coverage'] / (data['end'] - data['start'])
    data['above20x'] = data['mean_cov'].apply(lambda x: int(x > 20.0))
    data['above10x'] = data['mean_cov'].apply(lambda x: int(x > 10.0))
    data['below20x'] = data['mean_cov'].apply(lambda x: int(x < 20.0))
    data['below10x'] = data['mean_cov'].apply(lambda x: int(x < 10.0))
    transcript_coverage = data.groupby('gene')

    # get size of each group:
    t_mean = transcript_coverage.mean()
    t_mean['nr_exons'] = transcript_coverage.size()
    t_mean['below10_count'] = transcript_coverage.sum().loc[:, 'below10x']
    t_mean['below20_count'] = transcript_coverage.sum().loc[:, 'below20x']

    return data, t_mean


def aggregate_results(t_mean, gene_panel, name):
    columns = ['mean_cov', 'above20x', 'above10x', 'below20x', 'below10x', 'nr_exons', 
               'below20_count', 'below10_count']
    t_mean['Gene of interest'] = t_mean.index.isin(gene_panel).astype(int)
    columns.append('Gene of interest')
    with pd.ExcelWriter(name + '.summary.xlsx') as writer:
        sheet_name='All canonical transcripts'
        t_mean.loc[:, columns].to_excel(writer, sheet_name)
        excel_columns = [{'header': 'gene'}] + [{'header': header} for header in columns]
        options = {'columns': excel_columns, 'style': 'Table Style Light 11'}
        writer.sheets[sheet_name].add_table(0, 0, t_mean.loc[:, columns].shape[0], t_mean.loc[:, columns].shape[1], options)


if __name__ == '__main__':
    print("# Running exon coverage function")
    global_modules = globals()
    modules = version.imports(global_modules)
    version.print_modules(list(modules))
    filename = sys.argv[1]
    gene_file = sys.argv[2]
    gene_panel = list(pd.read_csv(gene_file).values.flatten())
    sample = filename.replace('.bed', '')
    data = read_input(filename)
    data, t_mean = calculate_statistics(data)
    aggregate_results(t_mean, gene_panel, name=sample)


