#!/usr/bin/env python3
# coding=utf-8
import os
import time
import argparse
import re
from pathlib import Path

def get_args(args=None):
    """
    Parse command line arguments to a namespace with values. 
    By default used with args=None to run it from command line but can be run using another script by setting args.
    :param args: string arguments or list of arguments.
    :return: namespace
    """
    # check function arguments
    if isinstance(args, str): args = args.split(' ')

    parser = argparse.ArgumentParser(description="A general submitter script. This should ideally be used for most submissions.")
    parser.add_argument("script", help="A string indicating the code to run.")
    parser.add_argument("-n", "--name", dest="name", default="icope", 
                        help="Name of the submission job. Default is a unique number after 'icope' that doesn't create conflict.")
    parser.add_argument("-no-nr", "--no-numbering", action='store_true', dest="disable_numbering", help="Disable addinga number after the name in filenames to avoid overwriting files.")
    parser.add_argument("-minutes", "--minutes", dest="minutes", type=int, default=0,
                        help="Wall-time minutes. Default is 0 unless no kind of wall-time is provided, then it is 30.")
    parser.add_argument("-hours", "--hours", type=int, default=10, help="Wall-time hours. Default is 0.")
    parser.add_argument("-mem", "--memory", help="Memory in gb. Default is 100.", default=100)
    parser.add_argument("-dir", "--workdir", default=os.getcwd(), 
                        help="The submission files are saved here. Relative filenames in command is relative to this. "
                             "Default is current directory.")
    parser.add_argument("-py2", "--python2", action='store_true', help="A flag set if we want to run the code with python2.")

    parser.add_argument("-w", "--wait", type=int, dest="wait_for",
                        help="Used to wait of a given job to finish successfully. Insert the job number.")
    parser.add_argument("-T", "--tunnel", dest="tunnel", action="store_true",
                        help="Flag set to include opening and closing tunnel for Sentieon pipeline. \
                        Default name are sentieonstart.sh and sentieonstop.sh located in your $HOME")
    parser.add_argument("-R", "--reserve", dest="reserve", action="store_true",
                        help="Use this flag to submit job to our own reserves nodes")
    parser.add_argument("-a", "--array",
                        help="Used to create job array, the job number will be $PBS_ARRAYID. Insert numbers to run eg. 608-631.")
    parser.add_argument("-max_jobs", "--max_jobs", type=int, default=48, 
                        help="Choose maximum number of jobs to run at a time from this call. "
                             "This setting is only relevant for array jobs. "
                             "If more array jobs are started than this number, some of them will wait to run.")
    parser.add_argument("--verbose", action="store_true", help="Add to execute set -xv and more information.")
    parser.add_argument("-np", "--nproc", default=28, dest="nproc", help="Number of processors to use")
    parser.add_argument("-test", "--test", dest="dry_run", action="store_true",
                        help="For testing writing the qsub-script only. Will not submit to Computerome.")
    parser.add_argument("-move", "--move-outfiles", dest="move_outfiles", action="store_true",
                        help="When running the pipeline, enable moving the output files to the sample folder generated "
                             "by the scirpt. The name of the generated folder is inferred from the pipeline version "
                             "and input sample name")
    args = parser.parse_args(args)
    if args.disable_numbering: 
        args.name = get_job_name(args.workdir, args.name, False)
    else: 
        args.name = get_job_name(args.workdir, args.name, True)
    print("# job name set to", args.name)
    print(args.workdir)

    if not args.minutes and not args.hours:
        args.minutes = 30

    if args.array:
        if "$PBS_ARRAYID" not in args.script:
            parser.error("When using -a the script needs to contain $PBS_ARRAYID")

    return args




def path_split_extension(fname):
    """
    Like os.path.splitext but splits before a double extension if one is present.
    It also splits e.g. ".tsv" as ("", ".tsv") so it will be seen as an extension.
    :param fname: string name
    :return: tuple (basename, extension)
    """
    root, extension = os.path.splitext(fname)
    if not extension and root.startswith("."):
        extension = root
        root = ""
    root, extra_extension = os.path.splitext(root)
    # an extension has to be a dot followed by only letters, if there is e.g. underscore it is seen as a basename.
    if re.match('\.[a-zA-Z]+', extra_extension):
        # it is in fact another extension
        return root, extra_extension + extension
    # it is not a real extension
    return root + extra_extension, extension


def number_files(filenames, max_number=99):
    """
    Get filenames for each given that has a number appended. The number increases to avoid file conflict,
    so we don't overwrite existing files.
    :param filenames: single filename or list of filenames to give a number (the same number for all)
    :param max_number: maximum number to give. If this is reached it will be used even if it overwrites.
    :return: filename(s) to use as single string or list of strings
    """
    is_single = isinstance(filenames, str)
    if is_single: filenames = [filenames]

    for i in range(max_number + 1):
        names = [filename_suffix(name, str(i).rjust(len(str(max_number)), '0')) for name in filenames]
        exists = [os.path.isfile(name) for name in names]
        if any(exists): continue
        # non exists if we get to this point
        if is_single: return names[0]
        return names


def filename_suffix(fname, suffix):
    """Add an ending after a '_' to a filename before the file extension.
    If there is a period in the ending, it is assumed this means the extension is chosen to be the part after the period.
    If the ending is only an extension, a change of extension is performed.
    :param fname string
    :param suffix string or number
    :return string"""
    # haven't decided how to deal with reading from stdin exactly.
    # Ideally we would read the name of the file being piped with cat command
    if fname == "/dev/stdin": fname = ''
    fname = fname.rstrip('_')
    suffix = str(suffix).strip('_')
    assert suffix
    root, extension = path_split_extension(fname)
    if root and not suffix.startswith('.'): root += '_'
    if '.' in suffix: return root + suffix
    if extension: return root + suffix + extension
    return root + suffix + '.fa'


def get_job_name(workdir, name, numbering=False):
    """
    Get a job name where we look for name conflict, have possibility of adding unique numbering.
    We also make sure you have not included the icope_ prefix already in the name since it is added later.
    :param workdir: working directory where we save the file
    :param name: initial name
    :param numbering: boolean indicating if we want to append numbering to make name unique.
    :return: possibly corrected name
    """
    
    root = workdir + "/" + name
    # find out what the outfiles are gonna be named to look for file conflict
    outfnames = [root + ext for ext in ['.qsub', '.out', '.err']]
    # append incremental numbers to avoid conflict
    if numbering: outfnames = number_files(outfnames)
    elif any([os.path.isfile(f) for f in outfnames]): print("# submission files will be overwritten")
    # get the name part of one of the resulting filename
    return os.path.splitext(os.path.basename(outfnames[0]))[0]


def get_walltime(hours=0, minutes=0, seconds=0):
    walltime = "%d:%2d:%2d" % (hours, minutes, seconds)
    walltime = walltime.replace(' ', '0')
    return walltime


def configure_outfiles(workdir, script, script_name, move_outfiles=False):
    if move_outfiles:
        pipeline_folder = '/home/projects/HT2_leukngs/apps/github/code/pipeline/'
        pipeline = script.split(' ')[0]
        sample_name = script.split(' ')[1]  # purposely ony take two first element as we can have a third for nproc
        version_suffix = os.readlink(pipeline_folder + pipeline).split('_')[-1].replace('.sh','')
        outbase = workdir + '/' + sample_name + '.' + version_suffix
    else:
        outbase = workdir
    return outbase


def write_tunnel(homedir):
    line = "#PBS -l prologue={dir}/sentieonstart.sh\n" \
           "#PBS -l epilogue={dir}/sentieonstop.sh\n".format(dir=homedir)
    return line

def get_array_PBS(array, n_jobs):
    if array: return "#PBS -t " + array + "%" + str(n_jobs)
    return ""


def write_qsub(name, script, out_base, nproc=1, memory=20, walltime='1:00:00', workdir=None, python=3, wait_for=None,
               extra_PBS="", reserve=False, verbose=False, tunnel=False, move_outfiles=False):
    """
    Write a qsub file that can be run on computerome.
    :param name: name of the job
    :param script: string that is the code to run.
    :param nproc: number of processors to use.
    :param memory: memory in gb to use.
    :param walltime: walltime string indicating walltime to allocate to job.
    :param workdir: the path to tun the code in. Default is using current working directory.
    :param python: version of python to use.
    :param extra_PBS: extra PBS lines to add at the end of the usual ones.
    :return: the filename of the qsub file.
    """

    if not isinstance(memory, str): memory = str(memory) + 'gb'
    icope_account='HT2_leukngs'
    qsub_string = \
        '#!/bin/sh\n' \
        '### Note: No commands may be executed until after the #PBS lines\n' \
        '### Account information\n' \
        '#PBS -W group_list={account} -A {account}\n'.format(account=icope_account)

    if wait_for:
        qsub_string += \
            '#PBS -W depend=afterok:{wait_for}\n'.format(wait_for=wait_for)
    if reserve:
        print("Submitting to nodes reserved for HT2_leukngs")
        qsub_string += \
            '#PBS -l advres=HT2_leukngs.916947\n'

    qsub_string += \
        '### Job name (comment out the next line to get the name of the script used as the job name)\n' \
        '#PBS -N {name}\n' \
        '### Output files (comment out the next 2 lines to get the job name used instead)\n' \
        '#PBS -e {out_base}/{name}.err\n' \
        '#PBS -o {out_base}/{name}.out\n' \
        '### Email: no (n)\n' \
        '#PBS -M n\n' \
        '### Make the job rerunable (y)\n' \
        '#PBS -r y\n' \
        '### Number of nodes\n' \
        '#PBS -l nodes=1:ppn={nproc}:thinnode\n' \
        '#PBS -l walltime={walltime}\n' \
        '#PBS -l mem={memory}\n' \
        '{extra}\n' \
        'echo This is the STDOUT stream from a PBS Torque submission script.\n' \
        '# Go to the directory from where the job was submitted (initial directory is $HOME)\n' \
        'echo Working directory is $PBS_O_WORKDIR\n' \
        'cd $PBS_O_WORKDIR\n' \
        '\n' \
        '# Get number of processors\n' \
        'export NPROC={nproc}\n' \
        'echo "This job has allocated $NPROC nodes"\n' \
        '\n' \
        '# Load user Bash settings:\n' \
        'source /home/projects/HT2_leukngs/apps/github/shared_utils/shared_bash_profile\n' \
        ''.format(extra=extra_PBS, name=name, nproc=nproc, out_base=out_base, memory=memory,
                  walltime=walltime, workdir=workdir)

    if '2' in str(python):
        qsub_string += \
            'module unload anaconda3\n' \
            'module load anaconda2/4.0.0\n'

    qsub_string += \
        '\n' \
        'start=`date +%s` \n' \
        'echo "Now the user defined script is run. After the ---- line, the STDOUT stream from the script is pasted."\n' \
        'echo "Start at `date`"\n' \
        'echo "-----------------------------------------------------------------------------------------------------"\n' 
    
    # print all bash commands run and what they were interpreted as
    if verbose: qsub_string += 'set -vx\n'
        
    qsub_string += '\n' + script + '\n \n'
    if move_outfiles:
        qsub_string += "mv $PBS_O_WORKDIR/$PBS_JOBNAME.qsub {}".format(out_base)
    qsub_string += \
        '\n' \
        'echo "-----------------------------------------------------------------------------------------------------"\n' \
        'echo "End at `date`"\n' \
        'end=`date +%s` \n' \
        'runtime=$((end-start))\n' \
        'echo Runtime: $runtime seconds\n' \
        'sleep 5\n' \
        'exit 0\n'

    qsub_filename = '{workdir}/{name}.qsub'.format(workdir=workdir, name=name)
    with open(qsub_filename, 'w') as outfile: outfile.write(qsub_string)

    return qsub_filename


def submit(fname):
    start_dir = os.getcwd()
    os.chdir(os.path.dirname(fname))
    time.sleep(1)
    os.system('qsub %s' % os.path.abspath(fname))
    time.sleep(1)
    os.chdir(start_dir)


def main(args):
    walltime = get_walltime(args.hours, args.minutes)
    python = 2 if args.python2 else 3

    extra_string = get_array_PBS(args.array, args.max_jobs)
    if args.tunnel:
        home = str(Path.home())
        print("# Configuring tunnel with scripts in", home) 
        extra_string += write_tunnel(home)
    if args.verbose: print("write qsub")

    workdir = os.path.abspath(args.workdir) if args.workdir else os.getcwd()
    out_base = configure_outfiles(workdir, args.script, args.name, args.move_outfiles)

    fname = write_qsub(args.name, args.script, out_base, args.nproc, args.memory, walltime, args.workdir, python,
                       args.wait_for, extra_string, args.reserve, args.verbose, args.tunnel, args.move_outfiles)
    if args.verbose: print("submit qsub")
    if not args.dry_run: 
        submit(fname)
        print("# submitted")


if __name__ == "__main__":
    print("# starting submission...")
    parsed_args = get_args()
    print("# args:", parsed_args)
    main(parsed_args)
    print("# Done!")

