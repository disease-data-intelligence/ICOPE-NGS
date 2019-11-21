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


def plot_qual(filename, data):
    cols = ['[4]number of SNPs', '[7]number of indels']
    pretty_names = {cols[0]: 'SNPs',
                    cols[1]: 'Indels'}
    # for plotting
    for col in cols:
        fig, ax = plt.subplots()
        ax = data.plot.bar(x='[3]Quality', y=col, label='Number of '+pretty_names[col], logy=True, color='#007b7f')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(2.5))
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.set_title("Quality score distribution for " + pretty_names[col]
                     + ' ' + filename.split('.')[0].replace('qual', ''))
        ax.set_xlabel('Quality score')
        plt.setp(ax.patches, linewidth=0)
        plt.tight_layout()
        plot_name = filename.replace('txt', '')+pretty_names[col]+'.pdf'
        plt.savefig(plot_name, dpi=150)
        print(f"Saved plot to {plot_name}")


def plot_dp(filename, data):
    cols = ['[6]number of sites']
    pretty_names = {cols[0]: 'Number of sites'}
    # for plotting
    for col in cols:
        fig, ax = plt.subplots()
        ax = data.plot.bar(x='[3]bin', y=col, label=pretty_names[col], logy=True, color='#ff8f8f')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(2.5))
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.set_xlim(0, 500+5)
        ax.set_title("Depth distribution for " + pretty_names[col] + ' ' + filename.split('.')[0].replace('dp',''))
        ax.set_xlabel('Depth')
        plt.setp(ax.patches, linewidth=0)
        plt.tight_layout()
        plot_name = filename.replace('txt', '') + pretty_names[col].replace(' ', '_') + '.pdf'
        plt.savefig(plot_name, dpi=150)
        print(f"Saved plot to {plot_name}")


def summarize_sn(filename, data):
    # for SN
    data.set_index('[3]key')
    data.drop(columns=['# SN', '[2]id'], inplace=True)
    file_name = filename[0].split('-')[0] + '.xlsx'
    data.to_excel(file_name)
    print(f"Saved excel to {filename}")


def open_vcf_stats(filename, function):
    assert os.path.exists(filename), filename + "does not exist"
    identifier = '^' + function
    input = subprocess.Popen(["grep", identifier], stdin=open(filename), stdout=subprocess.PIPE)
    decoding = StringIO(input.communicate()[0].decode('utf-8'))
    data = pd.read_csv(decoding, sep='\t')
    return data


if __name__ == '__main__':
    """ Function for plotting either quality distribution or depth distribution for variants in true set. It can also 
    be used for making an excel-file containing information on several VCF-files"""
    filename = sys.argv[1]
    global_modules = globals()
    modules = imports(global_modules)
    print_modules(list(modules))
    print("# Gettings stats for", filename)
    fields_fun = {'SN': summarize_sn, 'QUAL': plot_qual, 'DP': plot_dp}
    for field, fun in fields_fun.items():
        data = open_vcf_stats(filename, field)
        fun(filename, data)
