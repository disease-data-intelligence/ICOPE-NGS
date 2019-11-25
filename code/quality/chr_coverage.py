#! /usr/bin/env python3

import pandas as pd 
import seaborn as sns 
import sys 
import numpy as np 
import pdb
import matplotlib
import subprocess
import os
from io import StringIO
import argparse

matplotlib.use('Agg')
import matplotlib.pyplot as plt 
from utils_py.version import print_modules, imports
from utils_py.pprinting import print_overwrite


def get_parser():
    parser = argparse.ArgumentParser(
        description="Run coverage stats for all")
    parser.add_argument('-in', dest="infile", help="Coverage file from samtools or bam-file", type=str)
    parser.add_argument('-out', dest='out', help="Name of outfile. Will be placed in current working directory if "
                                                 "there is no path.")
    parser.add_argument('-limit', dest='limit', help="Limit for plotting coverage. Default: 150", default=150)
    return parser


def get_args(args=None):
    parser = get_parser()
    args = parser.parse_args(args)
    return args



def run_genome_cov(bam):
    print("# Input is a bam-file, we have to run genomecov ... ")
    try:
        input = subprocess.Popen(["bedtools", "genomecov", "-ibam", bam, "-max", "150"], stdout=subprocess.PIPE)
    except OSError as e:
        print("# Could not run bedtools or find bam, try to use: module load bedtools/2.28.0 and check path")
    decoding = StringIO(input.communicate()[0].decode('utf-8'))
    cov = pd.read_csv(decoding, sep='\t', header=None)
    cov.columns = ['chr', 'cov', 'obs_bases', 'total', 'frac']
    return cov


def read_coverage_file(covfile):
    cov = pd.read_csv(covfile, sep='\t')
    print("Succesfully read ", filename)
    cov.columns = ['chr', 'cov', 'obs_bases', 'total', 'frac']
    return cov


def plot_collect_coverage(cov, outname, input_upper_limit):
    upper_limit = 50
    plots = cov['chr'].unique()
    # only plot chr 1-22, X, Y and genome
    plots = [x for x in plots if (x[0].isdigit() or x.startswith('X') or
                                  x.startswith('Y') or x.startswith('g'))]
    plots = [x for x in plots if (x[0].isdigit() or x[0] in ['X', 'Y', 'g', 'M'])]
    summary = dict()
    fig, axes = plt.subplots(int(np.ceil(len(plots) / 3)), 3, figsize=(30, len(plots)))
    for frag, ax in zip(plots, fig.axes):
        print_overwrite("# Now plotting region: ", frag)
        covered_fraction = cov[cov['chr'] == frag][0:upper_limit]['frac'].sum()
        skip = 5
        while covered_fraction < 0.90 and upper_limit <= input_upper_limit:
            upper_limit += 1
            covered_fraction = cov[cov['chr'] == frag][0:upper_limit]['frac'].sum()
            # print("# Recomputed upper limit to:", upper_limit, "covered fraction", covered_fraction)
        if upper_limit > 100:
            skip = 10
        sns.barplot(x='cov', ax=ax, y='frac', data=cov[cov['chr'] == frag][0:upper_limit])
        total = cov[cov['chr'] == frag].iloc[0, 3]
        mean_cov = (cov[cov['chr'] == frag].loc[:, 'obs_bases'] * cov[cov['chr'] == frag].loc[:, 'cov']).sum() / total
        summary[frag] = mean_cov
        ax.set_title(
            "Coverage distribution for chromosome / contig " + frag + ". Mean coverage=" + str(round(mean_cov, 3)))
        # set xticks
        ax.set_xticklabels(ax.get_xticklabels(), visible=False)

        for label in ax.xaxis.get_ticklabels()[::skip]:
            label.set_visible(True)
    print("\n# Done plotting ...")
    fig.tight_layout()
    outname = outname + '_coverage_pr_chromosome.png'
    print("# saving coverage pr. chromosome plot to", outname)
    plt.savefig(outname, format='png')
    summary_df = pd.DataFrame(summary, index=[0]).transpose().rename(columns={0: 'Coverage'})
    return summary_df


def main(filename, input_upper_limit, outname=None):
    split = os.path.splittext(filename)
    if outname:
        outname = os.path.join(os.getcwd(), outname)
    else:
        outname = split[0]

    if split[1] == '.bam':
        cov = run_genome_cov(filename)
    elif split[1] == '.cov':
        cov = read_coverage_file(filename)
    print(f"# Saving to {outname}")
    mean_coverage = plot_collect_coverage(cov, outname, input_upper_limit)
    mean_coverage.to_csv(outname + '.tsv', sep='\t')


if __name__ == '__main__':
    print("# Running coverage pr. chromosome function") 
    global_modules = globals()
    modules = imports(global_modules)
    print_modules(list(modules))
    args = get_args()
    print(f"# Input: {args.infile} \t  Upper limit: {args.limit}")
    main(args.infile, args.limit, args.out)
    print("# Done!")






