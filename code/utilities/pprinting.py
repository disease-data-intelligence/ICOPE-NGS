#! /usr/bin/env python3

import sys


def print_overwrite(*line, **kwargs):
    """
    Print by overwriting the last print line.
    :param line: string line to print
    :param kwargs: keyword can be sep or delimiter to look for separator character between each element given in line
    :return: 
    """
    if "sep" in kwargs: sep = kwargs["sep"]
    elif "delimiter" in kwargs: sep = kwargs["delimiter"]
    else: sep=" "
    sys.stdout.write("\r" + ''.join(line))
    sys.stdout.flush()
