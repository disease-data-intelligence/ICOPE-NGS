#! /usr/bin/env python3

import pandas as pd 
import seaborn as sns 
import sys 
import numpy as np 
import pdb
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt 
sys.path.append("/home/projects/HT2_leukngs/apps/github/code/utilities")
import version

def main(file):
    cov = pd.read_csv(file, sep='\t')
    cov.columns = ['chr', 'cov', 'obs_bases', 'total', 'frac']
    upper_limit = 50
    plots = cov['chr'].unique()
    # only plot chr 1-22, X, Y and genome
    plots = [x for x in plots if (x[0].isdigit() or x.startswith('X') or
                                  x.startswith('Y') or x.startswith('g'))]

    fig, axes = plt.subplots(int(np.ceil(len(plots) / 3)), 3, figsize=(30, len(plots)))
    for frag, ax in zip(plots, fig.axes):
        covered_fraction = cov[cov['chr'] == frag][0:upper_limit]['frac'].sum()
        skip = 5
        while covered_fraction < 0.90:
            upper_limit += 1
            covered_fraction = cov[cov['chr'] == frag][0:upper_limit]['frac'].sum()
            # print("# Recomputed upper limit to:", upper_limit, "covered fraction", covered_fraction)
        if upper_limit > 150:
            skip = 10
        sns.barplot(x='cov', ax=ax, y='frac', data=cov[cov['chr'] == frag][0:upper_limit])
        total = cov[cov['chr'] == frag].iloc[0, 3]
        mean_cov = (cov[cov['chr'] == frag].loc[:, 'obs_bases'] * cov[cov['chr'] == frag].loc[:, 'cov']).sum() / total
        ax.set_title(
            "Coverage distribution for chromosome / contig " + frag + ". Mean coverage=" + str(round(mean_cov, 3)))
        # set xticks
        ax.set_xticklabels(ax.get_xticklabels(), visible=False)

        for label in ax.xaxis.get_ticklabels()[::skip]:
            label.set_visible(True)

    fig.tight_layout()
    outname = 'coverage_pr_chromosome.png'
    plt.savefig(outname, format='png')


if __name__ == '__main__':
    print("# Running coverage pr. chromosome function") 
    global_modules = globals()
    modules = version.imports(global_modules)
    version.print_modules(list(modules))
    filename = sys.argv[1]
    print("# input:", filename)
    main(filename)






