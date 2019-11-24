#!/usr/bin/env python3

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
import pdb
import matplotlib
matplotlib.use('Agg')
import subprocess
import os
from io import StringIO

# import from own repos
from utils_py.version import print_modules, imports


def plot_qual(data, axes):
    cols = ['[4]number of SNPs', '[7]number of indels']
    pretty_names = {cols[0]: 'SNPs',
                    cols[1]: 'Indels'}
    # for plotting
    for col, ax in zip(cols, axes):
        data.plot.bar(x='[3]Quality', y=col, label='Number of '+pretty_names[col], logy=True, color='#007b7f', ax=ax)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(2.5))
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.set_title("Quality score distribution for " + pretty_names[col])
        ax.set_xlabel('Quality score')
        plt.setp(ax.patches, linewidth=0)
        plt.tight_layout()


def plot_dp(data, ax):
    if len(data) == 0:
        # some vcf files do not have this field, data will be empty
        return None
    else:
        cols = ['[6]number of sites']
        pretty_names = {cols[0]: 'Number of sites'}
        # for plotting
        for col in cols:
            data.plot.bar(x='[3]bin', y=col, label=pretty_names[col], logy=True, color='#ff8f8f', ax=ax)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(2.5))
            ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
            ax.set_xlim(0, 500+5)
            ax.set_title("Depth distribution for " + pretty_names[col])
            ax.set_xlabel('Depth')
            plt.setp(ax.patches, linewidth=0)
            plt.tight_layout()


def summarize_sn(data, filename):
    # for SN
    data = data.set_index('[3]key')
    data.drop(columns=['# SN', '[2]id'], inplace=True)
    sheetname = filename.replace('txt', 'xlsx')
    data.to_excel(sheetname)
    print(f"Saved excel to {sheetname}")


def open_vcf_stats(filename, function):
    assert os.path.exists(filename), filename + "does not exist"
    identifier = '^' + function
    input = subprocess.Popen(["grep", identifier, "-B", "1"], stdin=open(filename), stdout=subprocess.PIPE)
    decoding = StringIO(input.communicate()[0].decode('utf-8'))
    data = pd.read_csv(decoding, sep='\t')
    return data


if __name__ == '__main__':
    """ Function for plotting either quality distribution or depth distribution for variants in true set. It can also 
    be used for making an excel-file containing information on several VCF-files"""
    filename = sys.argv[1]
    # to consider: subprocess.Popen("bcftools", "stats", filename, ">" ""$destination/"$opt"vcf_summary.txt)
    global_modules = globals()
    modules = imports(global_modules)
    print_modules(list(modules))
    print("# Gettings stats for", filename)
    fig, axes = plt.subplots(3, figsize=(8,15))
    fields_fun = [summarize_sn, plot_qual, plot_dp]
    fields_args = [filename, axes[1:], axes[0]]
    field_names = ['SN', 'QUAL', 'DP']
    data = map(lambda x: open_vcf_stats(filename,x), field_names)
    for fun, data, args in zip(fields_fun, data, fields_args):
        fun(data, args)

    plot_name = filename.replace('txt', 'pdf')
    plt.savefig(plot_name, dpi=150)
    print(f"Saved plot to {plot_name}")
