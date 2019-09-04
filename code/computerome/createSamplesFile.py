#!/usr/bin/env python3
# coding=utf-8

# *********************
# Python script that will take directory names
# and find the corresponding fastqc files 
# They and their paths with get printet into a file
# that can then be used to submit multiple jobs
# NB: Change the path to outputfile to wanted directory 
# ********************** 

import os
import sys

# Calculates number of samples given in commandline
numSamples = len(sys.argv)-2

# Error if germline (G) or tumor (T) is not specified
type = str(sys.argv[1])
if (type != 'G' and type != 'T'):
   print('First argument must be G or T, please try again') 
   sys.exit(1)

# Path to the samples
pathSamples = '/home/projects/HT2_leukngs/data/samples'

# Opens file that will contain the files with corresponding paths
try:
    outFile = open('/home/projects/HT2_leukngs/people/s174333/samplesFile.txt', 'w')
except IOError as error:
    print('Error when trying to open file. Reason: ', str(error))

# Looking for fastqc file for each sample printing it into the outputfile
for i in range(numSamples):
    sample = str(sys.argv[i+2])
    path = pathSamples + '/' + sample
    fastqc = ''
    # r=root, d=directories, f=files
    for r, d, f in os.walk(path):
        for file in f:
            if ('.R1.fq' in file and sample in file and type in file):
               fastqc = os.path.join(r, file)
               print(fastqc)
               print(fastqc, file=outFile)
    if fastqc == '':
        print('Could not find fastqc file for sample:', sample)        

outFile.close()


