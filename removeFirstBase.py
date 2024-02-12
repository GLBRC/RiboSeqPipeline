#!/usr/bin/env python
"""removeFirstBase.py

Remove the first base in each read if the is low. Low defined as a Q score <= 10,
['+', '*', ')','(', "'", '&', '%', '$', '#', '"', '!']. Part of the Ribo-Seq pipeline.

Notes
-----
Runs in the current directory and globs the file names in the cutadapt directory.

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
import sys

parentDir = os.getcwd() + '/'

def removeFirstBase():
    """removeFirstBase
    Remove reads where the first position quality score is <=10
        """
    qualScores = ['+', '*', ')','(', "'", '&', '%', '$', '#', '"', '!']
    for fastq in glob.glob('cutadapt/*-clean.fastq'):
        name = re.sub('-clean.fastq', '-filt.fastq', fastq)
        with open(fastq, 'r') as f, open(name, 'w') as out:
            for hdr, seq, plus, qual in itertools.zip_longest(*[f]*4):
                firstQual = [*qual][0]
                if not firstQual in qualScores:
                    out.write(hdr)
                    out.write(seq)
                    out.write(plus)
                    out.write(qual)
                else:
                    print(hdr, '  ', qual)    

def main():
    
    cmdparser = argparse.ArgumentParser(description="Remove first base for each read in fastq files if Q <= 10.",
                                        usage='%(prog)s ' ,prog='removeFirstBase.py'  )
    cmdResults = vars(cmdparser.parse_args())
    
    # get the fastq files in the cutadapt directory.
    for fastq in glob.glob('cutadapt/*-clean.fastq'):
        print('name : ', fastq)
   

    


if __name__ == "__main__":
    main()

