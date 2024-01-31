#!/usr/bin/env python
"""RiboSeqPipeline.py

Pipeline to align and count start sites for Ribo-Seq data.

Notes
----- 
Put the fastq files in the parent directory within a directory called fastq.

use the gatk4 environment: /home/glbrc.org/mplace/.conda/envs/gatk4

Method
------
    
Parameters
----------
f : str
    A text file with a list of fastq files to process. Only read 1 if paired-end data.

Example
-------
    usage:

        RiboSeqPipeline.py -f fastqfiles.txt
    
Requirements
------------

    1. HTCondor, specifically GLBRC scarcity compute cluster.

References
----------

"""
import argparse 
import itertools
import os
import re
import subprocess
import sys
from interval_tree import IntervalTree

parentDir = os.getcwd() + '/'
referenceDir = '/home/glbrc.org/mplace/scripts/riboSeqPipeline/reference/'

print(parentDir)
print(referenceDir)

def main():
    
    cmdparser = argparse.ArgumentParser(description="Ribo-Seq pipeline, produces alignments and counts.",
                                        usage='%(prog)s -f <fastqFileList.txt>' ,prog='RiboSeqPipeline.py'  )
    cmdparser.add_argument('-f', '--file',    action='store', dest='FILE',    help='Text file, one fastq file name per line, read 1 only if paired-end.', metavar='')
    cmdResults = vars(cmdparser.parse_args())
        
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
    
    if cmdResults['FILE'] is not None:
        fastqFile = cmdResults['FILE']
        fastqLst = []
        with open(fastqFile, 'r') as f:
            for ln in f:
                fastqLst.append(ln.rstrip())
                

    
    
                
        
    


if __name__ == "__main__":
    main()

