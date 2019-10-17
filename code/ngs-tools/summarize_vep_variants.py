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


def parse_format_rowwise(row, TNScope=False):
    fields = row['FORMAT'].split(':')
    if TNScope > 10:
        for field, normal, tumor in zip(fields, row['FORMAT_NORMAL'].split(':'), row['FORMAT_TUMOR'].split(':')):
            row[field + '_NORMAL'], row[field + '_TUMOR'] = normal, tumor
    else:
        for field, value in zip(fields, row['FORMAT_NORMAL'].split(':')):
            row[field] = value
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
    sample_name = sample.split('.filter')[0]
    input = subprocess.Popen(["zgrep", "-v", "#"], stdin=open(sample), stdout=subprocess.PIPE)
    decoding = StringIO(input.communicate()[0].decode('utf-8'))
    data = pd.read_csv(decoding, sep='\t', header=None)
    # detect format:
    TNScope = len(data.columns) > 10
    header_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
                   'FORMAT', 'FORMAT_NORMAL']
    drop_cols = ['INFO', 'FORMAT', 'FORMAT_NORMAL']
    if TNScope:
        header_cols.append('FORMAT_TUMOR')
        drop_cols.append('FORMAT_TUMOR')    # for later filtering
    else:
        data.columns = header_cols
    data['Sample'] = sample_name
    data = data.apply(parse_format_rowwise, args=TNScope, axis=1)
    data = parse_info(data)
    data.drop(columns=drop_cols, inplace=True)
    selected_fields = ['Sample', 'CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'ID', 'IMPACT', 'Consequence',
                       'Feature_type', 'SYMBOL', 'Feature', 'Gene', 'SIFT', 'PolyPhen']
    return data


def main(samples, outfile):
    all_data = pd.DataFrame()
    for s in samples:
        sample_variants = parse_vcf(s)
        all_data = pd.concat([sample_variants, all_data], sort=True)
    all_data.index.name = 'Unique variant index'
    all_data.to_csv(outfile, sep='\t')
    selected_fields = ['Sample', 'CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'ID', 'IMPACT', 'Consequence',
                       'Feature_type', 'Feature', 'Gene', 'SIFT', 'PolyPhen']
    print("# Number of annotations before filtering:", len(all_data))
    all_data = all_data[~(all_data['IMPACT'].isin(['LOW']) | all_data['Consequence'].isin(['synonymous_variant']) |
                          all_data['SIFT'].str.startswith('tolerated') | all_data['PolyPhen'].str.startswith('benign'))]
    print("# Number of annotations after filtering:", len(all_data))
    all_data.reindex(columns=selected_fields).to_csv(outfile.replace('.tsv', '_selected_fields.tsv'), sep='\t')



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
