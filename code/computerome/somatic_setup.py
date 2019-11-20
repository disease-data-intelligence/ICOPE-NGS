#! /usr/bin/env python3
import pdb
import os
import sys
import time
import argparse
import subprocess
from utils_py.version import print_modules, imports
from datetime import datetime


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
    return parser


def get_args(args=None):
    parser = get_parser()
    args = parser.parse_args(args)
    args.PSG_version = "PSG" + args.PSG_version
    args.PST_version = "PST" + args.PST_version
    args.PSP_version = "PSP" + args.PSP_version
    return args


def find_pairs(samples, germline_pipeline, tumor_pipeline, paired_pipeline):
    """
    Function for finding paired germline and tumor sample bam-files.
    yields pairs so that you can use it as an iterator.
    :param samples: foldername for walk start, str
    :param germline_pipeline: int
    :param tumor_pipeline: int
    :param paired_pipeline: int
    :return: path for tumor bam, germline bam and destination for somatic variants if used in concert with calling
    the somatic pipeline
    """
    for s in samples:
        print("# Finding files for sample", s)
        tumor, germline, dest = (None, None, None)
        for walk_idx, (r, d, f) in enumerate(os.walk(s, topdown=True)):
            if walk_idx == 0:
                assert len(d) == 2, "You are not submitting from right destination or sample is not paired"
            if r.endswith('.'+germline_pipeline):
                germline_name = r.split('/')[-1]
                germline = r + '/' + [x for x in f if (x.endswith('.bam') and x.startswith(germline_name))][0]
                print("# Found germline file",  germline)
            if r.endswith('.'+tumor_pipeline):
                tumor_name = r.split('/')[-1]
                tumor = r + '/' + [x for x in f if (x.endswith('.bam') and x.startswith(tumor_name))][0]
                tumor_path = '/'.join(r.split('/')[:-1])
                print("# Found tumor file", tumor)
        if sum(map(lambda x: isinstance(x, str), (tumor, germline))) == 2:
            mrd = tumor_name.split('_')[0]
            dest_filename = mrd + '-' + tumor_name.split('_')[1] + '-' + germline_name.split('_')[1] + '-' + paired_pipeline
            dest =  tumor_path + '/' + dest_filename
            yield (tumor, germline, dest)
        else:
            if not tumor:
                print("# Did not find tumor file")
            if not germline:
                print("# Did not find germline file")




def submit_pair(tumor, germline, dest):
    # write qsub with paired
    print("# Defining destination for somatic variants:", dest)
    apps='/home/projects/HT2_leukngs/apps/github/code'
    jobname = dest
    submit_string = "{a}/computerome/submit.py".format(a=apps)
    submit_string += " 'germline={g}\ntumor={t}\ndestination={d}\n".format(g=germline, t=tumor, d=dest)
    submit_string += "{a}/pipeline/sentieon_paired.sh $germline $tumor $destination' " \
                     "--hours 10 -n {j} -np 28 --move-outfiles --tunnel -no-nr".format(a=apps, j=jobname)
    subprocess.run(submit_string, shell=True, check=True)
    time.sleep(1)


def main(samples, psg, pst, psp):
    for pair in find_pairs(samples, psg, pst, psp):
        submit_pair(*pair)


if __name__ == "__main__":
    start_time = datetime.now()
    parsed_args = get_args()
    print("# args:", parsed_args)
    print("# Submitting paired jobs")
    global_modules = globals()
    modules = imports(global_modules)
    print_modules(list(modules))
    main(parsed_args.samples, parsed_args.PSG_version, parsed_args.PST_version, parsed_args.PSP_version)
    end_time = datetime.now()
    print("# Done!")
    print('# Duration: {}'.format(end_time - start_time))
