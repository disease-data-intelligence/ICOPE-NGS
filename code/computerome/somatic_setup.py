#! /usr/bin/env python3

import os
import sys
import time
import argparse
import subprocess
sys.path.append("/home/projects/HT2_leukngs/apps/github/code/utilities")
import version
from datetime import datetime


def get_parser():
    parser = argparse.ArgumentParser(
        description="Submitting paired tumor and germline analysis for Somatic variants\n"
                    "Assumes hardcoded Sentieon tumor and germline pipeline vs. 1 and only one bam-file for each.")
    parser.add_argument('-samples', dest="samples", help="Sample numbers", type=list)

    # define version of pipelines to use (default should be updated to use the latest version)
    parser.add_argument('-PSG', dest="PSG_version", type=str, default="01")
    parser.add_argument('-PSP', dest="PSP_version", type=str, default="01")
    parser.add_argument('-PST', dest="PST_version", type=str, default="02")

    return parser


def get_args(args=None):
    parser = get_parser()
    args = parser.parse_args(args)
    args.PSG_version = "PSG" + args.PSG_version
    args.PST_version = "PST" + args.PSG_version
    args.PST_version = "PSP" + args.PSP_version
    return args


def find_pairs(samples, germline_pipeline, tumor_pipeline, paired_pipeline):
    for s in samples:
        print("# Finding files for sample", s)
        tumor, germline, dest = (None, None, None)
        for walk_idx, (r, d, f) in enumerate(os.walk(s, topdown=True)):
            if walk_idx == 0:
                print(len(d))
                assert len(d) == 2, "You are not submitting from right destination or sample is not paired"
            if r.endswith('.'+germline_pipeline):
                name = r.split('/')[-1]
                germline = r + '/' + [x for x in f if (x.endswith('.bam') and x.startswith(name))][0]
                print("# Found germline file",  germline)
            if r.endswith('.'+tumor_pipeline):
                name = r.split('/')[-2:]
                tumor = r + '/' + [x for x in f if (x.endswith('.bam') and x.startswith(name))][0]
                dest = r + '/' + paired_pipeline
                print("# Found tumor file", tumor)
                print("# Defining destination for somatic variants:", dest)
            print(sum(map(lambda x: isinstance(x, str), (tumor, germline, dest))))
            if sum(map(lambda x: isinstance(x, str), (tumor, germline, dest))) == 3:
                yield (tumor, germline, dest)



def submit_pair(tumor, germline, dest):
    # write qsub with paired
    apps='/home/projects/HT2_leukngs/apps/github/code'
    jobname = dest
    submit_string = "{a}/computerome/submit.py '{a}/pipeline/sentieon_paired.sh {g} {t} {d}' " \
                    "--hours 2 -n {j} -np 5 -T -no-nr".format(a=apps, g=germline, t=tumor, d=dest, j=jobname)
    subprocess.run(submit_string, shell=True, check=True)
    time.sleep(1)


def main(samples):
    for pair in find_pairs(samples):
        submit_pair(*pair)


if __name__ == "__main__":
    print("# Submitting paired jobs")
    global_modules = globals()
    modules = version.imports(global_modules)
    version.print_modules(list(modules))
    start_time = datetime.now()
    parsed_args = get_args()
    print("# args:", parsed_args)
    main(parsed_args)
    end_time = datetime.now()
    print("# Done!")
    print('# Duration: {}'.format(end_time - start_time))
