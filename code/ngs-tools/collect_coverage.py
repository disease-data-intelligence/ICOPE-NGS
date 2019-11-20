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
sys.path.append("/home/projects/HT2_leukngs/apps/github/code/utilities")
import version
sys.path.append("/home/projects/HT2_leukngs/apps/github/code/computerome")
sys.path.append("/home/projects/HT2_leukngs/apps/github/code/quality")
# imports from own repo's
from utilities.version import print_modules, imports
from computerome.somatic_setup import find_pairs
from quality.gene_coverage import calculate_statistics


def get_parser():
    parser = argparse.ArgumentParser(
        description="Submitting paired tumor and germline analysis for somatic variants. "
                    "Assumes hardcoded Sentieon tumor (PST01) and germline pipeline (PSG01) and only one bam-file for "
                    "each. Pipeline version are the newest per default (this is updated recurrently in this script."
                    "Submission jobs do not request the reserved nodes.")
    parser.add_argument('-samples', dest="samples", nargs='+', help="Sample numbers", type=str)

    # define version of pipelines to use (default should be updated to use the latest version)
    parser.add_argument('-PSG', '--germline-pipeline-number', dest="PSG_version", type=str, default="01",
                        help="Germline pipeline number (Default: 01)")
    parser.add_argument('-PSP', '--paired-pipeline-number', dest="PSP_version", type=str, default="01",
                        help="Paired pipeline number (Default: 01). This defines number of new folders created and"
                             "should match what is linked to from sentieon_paired_frozen.sh")
    parser.add_argument('-PST', '--tumor-pipeline-number', dest="PST_version", type=str, default="02",
                        help="Tumor pipeline number (Default: 02)")
    parser.add_argument('-panel', '--gene-panel', dest='panel',
                        default='/home/projects/HT2_leukngs/data/references/general/300_genes_of_interest.txt',
                        help="(Default: /home/projects/HT2_leukngs/data/references/general/300_genes_of_interest.txt)")
    parser.add_argument('-bed', '--bedfile', type=str, dest='bed',
                        default="/home/projects/HT2_leukngs/data/references/hg37/USCS.hg37.canonical.exons.bed",
                        help="Bed-file with intervals to look at. "
                             "(Default: /home/projects/HT2_leukngs/data/references/hg37/USCS.hg37.canonical.exons.bed)")
    parser.add_argument('-destination', dest='destination')
    return parser


def get_args(args=None):
    parser = get_parser()
    args = parser.parse_args(args)
    args.PSG_version = "PSG" + args.PSG_version
    args.PST_version = "PST" + args.PST_version
    args.PSP_version = "PSP" + args.PSP_version
    if args.destination is None:
        args.destination = os.getcwd()
    return args


def run_samtools(bam, bed):
    assert os.path.exists(bam), "does not exist"
    assert os.path.exists(bed), "does not exist"
    input = subprocess.Popen(["samtools", "bedcov", bed, bam], stdout=subprocess.PIPE)
    decoding = StringIO(input.communicate()[0].decode('utf-8'))
    data = pd.read_csv(decoding, sep='\t', header=None)
    data.columns = ['chromosome', 'start', 'end', 'gene', 'exon', 'strand', 'coverage']
    data.drop(columns=['strand'], inplace=True)
    return data


def calculate_coverage_stats(data, panel):
    data['mean_cov'] = data['coverage'] / (data['end'] - data['start'])
    data['exons_above20x_frac'] = data['mean_cov'].apply(lambda x: int(x > 20.0))
    data = data[data['gene'].isin(panel)]
    coverage = data.groupby('gene').mean()
    return coverage.loc[:, ['mean_cov', 'exons_above20x_frac']]


def process_pair(germline, tumor, bed, destination, panel):
    germline_coverage_genes = calculate_coverage_stats(run_samtools(germline, bed), panel)
    tumor_coverage_genes = calculate_coverage_stats(run_samtools(tumor, bed), panel)
    # sample_name = germline.split('.bam')[0]
    merge = pd.merge(tumor_coverage_genes, germline_coverage_genes,
                     left_index=True, right_index=True, suffixes=('_germline', '_tumor'))
    sample = germline.split('.bam')[0]
    output = destination + sample + '.tsv'
    merge.to_csv(output, sep='\t')


def main(samples, psg, pst, psp, destination, bed, panel_file):
    gene_panel = list(pd.read_csv(panel_file).values.flatten())
    for germline, tumor, _ in find_pairs(samples, psg, pst, psp):
        process_pair(germline, tumor, destination, bed, gene_panel)


if __name__ == "__main__":
    start_time = datetime.now()
    args = get_args()
    print("# args:", args)
    print("# Submitting paired jobs")
    global_modules = globals()
    modules = imports(global_modules)
    print_modules(list(modules))
    main(args.samples, args.PSG_version, args.PST_version, args.PSP_version,
         args.destination, args.bed, args.panel)
    end_time = datetime.now()
    print("# Done!")
    print('# Duration: {}'.format(end_time - start_time))
