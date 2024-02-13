#!/usr/bin/env python
"""removeNonCoding.py

Create new fastq files for reads that DO NOT align to the Non-Coding 
reference genome.

notes
-----
-f, (--require-flags FLAG   ...have all of the FLAGs present)
Do not output alignments with any bits set in FLAG present in the FLAG field, 
only select unaligned reads, sam flag "4", i.e. the reads we want to analyze later.

Runs in the current directory and filters the non-coding alignments in the 
alignNonCodingRNA directory.

Parameters
----------
None

Example
-------
    usage:

        removeNonCoding.py
"""
import argparse 
import glob
import os
import re
import subprocess
import sys

nonCodingOutDir = os.getcwd() + '/alignNonCodingRNA/'

def filter(inSam, outFile):
    """filter
    
    Create new fastq files by filtering for reads that DID NOT ALIGN to the Non-Coding RNA.  
        
    samtools view -f 4 -u -O SAM -o $2 $1
    """
    # Get a list of reads which aligned to Non-Coding RNA and remove them from the *-filt.fastq files
    sampleName = re.sub('.sam', '', os.path.basename(inSam))
    outSam = nonCodingOutDir + sampleName + '-unmapped.sam'
    #cmd = [ 'samtools', 'view', '-f', '4', '-u', '-O', 'SAM', outFile, inSam ]
    #subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()                     
       
    '''
    # Use seqtk to subset the unmapped reads for use with bowtie2
    for fstq in glob.glob(cutadaptOutDir + '*-filt.fastq'):
        print('processing: ', fstq)
        outFile = parentDir + 'alignments/' +  re.sub('-filt.fastq', '.fastq', os.path.basename(fstq))    # create output file name
        nameLst = parentDir + 'alignments/' + re.sub('-filt.fastq', '-names.txt', os.path.basename(fstq)) 
        cmd = [ 'seqtk', 'subseq', fstq, nameLst ]
        # run command and capture output
        output = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()    
        # write results to file
        with open(outFile, 'w') as out:
            out.write(output[0].decode('utf-8'))
        out.close() 
    '''

def main():
    cmdparser = argparse.ArgumentParser(description="Create new fastq files using reads NOT aligned to reference NON-CODING rna.",
                                        usage='%(prog)s  ' ,prog='removeNonCoding.py'  )
    cmdparser.add_argument('-s', '--sam', action='store', dest='SAM',
                           help='Sam file to filter for unaligned reads', metavar='')
    cmdparser.add_argument('-o', '--output', action='store', dest='OUT',
                           help='Output fastq file, where reads are taken from input sam file.', metavar='')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
    
    # gather command line args
    if cmdResults['SAM'] is not None:
        inSam = cmdResults['SAM']
    else:
        print("Missing sam input file.")
        cmdparser.print_help()
        sys.exit(1)
    
    if cmdResults['OUT'] is not None:
        outFile = cmdResults['OUT']
    else:
        print()
        cmdparser.print_help()
        sys.exit(1)
        
    filter(inSam, outFile)
        
        
if __name__ == "__main__":
    main()
