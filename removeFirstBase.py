#!/usr/bin/env python
"""removeFirstBase.py

Remove the first base in each read if the is low. Low defined as a Q score <= 10,
['+', '*', ')','(', "'", '&', '%', '$', '#', '"', '!']. Part of the Ribo-Seq pipeline.

Notes
-----
Runs in the current directory and using the *-clean.fastq files in the cutadapt directory.

Parameters
----------
None

Example
-------
    usage:

        removeFirstBase.py
"""
import argparse 
import glob
import itertools
import os
import re
import sys

def removeFirstBase(fastq):
    """removeFirstBase
    Remove reads where the first position quality score is <=10
    """
    qualScores = ['+', '*', ')','(', "'", '&', '%', '$', '#', '"', '!']
    
    name = re.sub('-clean.fastq', '-filt.fastq', fastq)
    with open(fastq, 'r') as f, open(name, 'w') as out:
        for hdr, seq, plus, qual in itertools.zip_longest(*[f]*4):
            firstQual = [*qual][0]
            if not firstQual in qualScores:
                out.write(hdr)
                out.write(seq)
                out.write(plus)
                out.write(qual)

def main():
    cmdparser = argparse.ArgumentParser(description="Remove first base for each read in fastq files if Q <= 10.",
                                        usage='%(prog)s ' ,prog='removeFirstBase.py'  )
    cmdparser.add_argument('-f', '--fastq', action='store', dest='FASTQ', 
                           help='Fastq file to filter by first base.', metavar='')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
    
    # get command line arg
    if cmdResults['FASTQ'] is not None:
        fastq = cmdResults['FASTQ']
    else:
        print('Required fastq file missing.')
        cmdparser.print_help()
        sys.exit()    
    
    # get the fastq files in the cutadapt directory.
    removeFirstBase(fastq)

if __name__ == "__main__":
    main()
