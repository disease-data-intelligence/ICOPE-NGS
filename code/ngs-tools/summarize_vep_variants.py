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
#import version

def get_parser():
    parser = argparse.ArgumentParser(
        description="Combines many annotated vcf-files to one overview table. Only keeps few of the annotations and "
                    "few format attributes")
    parser.add_argument('-samples', dest="samples", nargs='+', help="Sample numbers", type=str)
    parser.add_argument('-outfile', dest="outfile", help="Name of output. (Default: variant_collection.tsv)",
                        default='variant_collection.tsv', type=str)
    return parser


def get_args(args=None):
    parser = get_parser()
    args = parser.parse_args(args)
    if not args.outfile.endswith('.tsv'):
        args.outfile += '.tsv'
    return args


def parse_format_rowwise(row):
    fields = row['FORMAT'].split(':')
    for field, normal, tumor in zip(fields, row['FORMAT_NORMAL'].split(':'), row['FORMAT_TUMOR'].split(':')):
        row[field + '_NORMAL'], row[field + '_TUMOR'] = normal, tumor
    return row


def parse_info(data):
    fields = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON",
              "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
              "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "SIFT",
              "PolyPhen"]
    data['INFO'] = data['INFO'].str.replace('CSQ=', '').str.split(',')
    data = data.explode('INFO')
    data[fields] = data['INFO'].str.split('|', expand=True)

    return data


def parse_vcf(sample):
    assert os.path.exists(sample), sample + "does not exist"
    sample_name = sample.split('.vep')[0]

    # read vcf-input
    input = subprocess.Popen(["grep", "-v", "#"], stdin=open(sample), stdout=subprocess.PIPE)
    decoding = StringIO(input.communicate()[0].decode('utf-8'))
    data = pd.read_csv(decoding, sep='\t', header=None)
    data.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
                    'FORMAT', 'FORMAT_NORMAL', 'FORMAT_TUMOR']
    data['Sample'] = sample_name
    data = data.apply(parse_format_rowwise, axis=1)
    data = parse_info(data)
    data.drop(columns=['INFO', 'FORMAT', 'FORMAT_NORMAL', 'FORMAT_TUMOR'], inplace=True)
    return data


def main(samples, outfile):
    all_data = pd.DataFrame()
    for s in samples:
        sample_variants = parse_vcf(s)
        all_data = pd.concat([all_data, sample_variants])
    all_data.to_csv(outfile, sep='\t')



if __name__ == "__main__":
    start_time = datetime.now()
    parsed_args = get_args()
    print("# args:", parsed_args)
    print("# Summarizing variants")
    global_modules = globals()
    #modules = version.imports(global_modules)
    #version.print_modules(list(modules))
    main(parsed_args.samples, parsed_args.outfile)
    end_time = datetime.now()
    print("# Done!")
    print('# Duration: {}'.format(end_time - start_time))
