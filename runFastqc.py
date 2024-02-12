#!/usr/bin/env python
"""runFastqc 

Run Fastqc on the trimmed fastq files.

Notes
-----
Input trimmed forward fastq file name.

Parameters
----------

f: str
    Name of fastq file to process, only the forward read is analyzed.
                          
"""
import argparse	
import os
import re
import subprocess
import sys

def main():
    """
    Main 
    """   
    cmdparser = argparse.ArgumentParser(description="Run Fastqc w/ rnaseq python environment.",
                                        usage='%(prog)s ' ,prog='runFastqc.py'  )                                  
    cmdparser.add_argument('-f', '--file', action='store', dest='FILE', help='Sequence file name read_1')
    cmdResults = vars(cmdparser.parse_args())

    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)

    # get fastq file names
    if cmdResults['FILE'] is not None:
        infile = cmdResults['FILE']
        print(infile )
    else:
        print('Required fastq file missing')
        sys.exit(1)

    # set up trimmomatic call
    cmd = [ 'fastqc', infile ]
    output = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    result = output[0].decode('utf-8')
    log    = output[1].decode('utf-8')

    # add logging here 
    print(result)
    print('\n\n')
    print(log)
    
if __name__ == "__main__":
    main()
