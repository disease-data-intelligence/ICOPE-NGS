#!/usr/bin/env python3

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
import pdb
import matplotlib 
matplotlib.use('Agg')

sys.path.append("/home/projects/HT2_leukngs/apps/github/code/utilities")
from utilities.version import print_modules, imports


def plot_qual(filename):
    qual = pd.read_csv(filename, sep='\t')
    cols = ['[4]number of SNPs', '[7]number of indels']
    pretty_names = {cols[0]: 'SNPs',
                    cols[1]: 'Indels'}
    # for plotting
    for col in cols:
        fig, ax = plt.subplots()
        ax = qual.plot.bar(x='[3]Quality', y=col, label='Number of '+pretty_names[col], logy=True, color='#007b7f')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(2.5))
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.set_title("Quality score distribution for " + pretty_names[col] + ' ' + filename.split('.')[0].replace('qual',''))
        ax.set_xlabel('Quality score')
        plt.setp(ax.patches, linewidth=0)
        plt.tight_layout()
        plt.savefig(filename.split('.')[0]+pretty_names[col]+'.pdf', dpi=150)


def plot_dp(filename):
    dp = pd.read_csv(filename, sep='\t')
    cols = ['[6]number of sites']   
    pretty_names = {cols[0]: 'Number of sites'}
    # for plotting
    for col in cols:
        fig, ax = plt.subplots()
        ax = dp.plot.bar(x='[3]bin', y=col, label=pretty_names[col], logy=True, color='#ff8f8f')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(2.5))
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.set_xlim(0, 500+5)
        ax.set_title("Depth distribution for " + pretty_names[col] + ' ' + filename.split('.')[0].replace('dp',''))
        ax.set_xlabel('Depth')
        plt.setp(ax.patches, linewidth=0)
        plt.tight_layout()
        plt.savefig(filename.split('.')[0]+pretty_names[col].replace(' ', '_') +'.pdf', dpi=150)


def summarize_sn(filename):
    # for SN
    sn_all = pd.DataFrame()
    for file in filename:
        sn = pd.read_csv(file, sep='\t', index_col='[3]key').drop(columns=['# SN', '[2]id'])
        sn.columns = [file.split('.')[0]]
        sn_all = pd.concat([sn_all, sn], axis=1)
    sn_all.to_excel(filename[0].split('-')[0]+'.xlsx')



if __name__ == '__main__':
    """ Function for plotting either quality distribution or depth distribution for variants in true set. It can also 
    be used for making an excel-file containing information on several VCF-files"""
    function = sys.argv[1].lower()
    filename = sys.argv[2:]  # may be a list of files

    global_modules = globals()
    modules = imports(global_modules)
    print_modules(list(modules))

    print("Input args: \n",
          "function:", function, "\n",
          "filename:", filename)

    if function == 'qual':
        plot_qual(filename[0])
    elif function == 'dp':
        plot_dp(filename[0])
    elif function == 'sn':
        summarize_sn(filename)




