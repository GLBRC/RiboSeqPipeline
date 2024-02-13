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
import sys

# Create new fastq files by filtering for reads that DID NOT ALIGN to the Non-Coding RNA
# First create a new file containing the read names (for reads we want to keep)
# nonCodingOutDir
for unmapped in glob.glob(nonCodingOutDir + '*-unmapped.sam'):
    nameFile = re.sub('-unmapped.sam', '', os.path.basename(unmapped))
    outFile  = parentDir + 'alignments/' + nameFile + '-names.txt'
    with open(unmapped) as f, open(outFile, 'w') as out:
        for line in f:
            name = line.split('\t')[0]
            out.write(f'{name}\n')
    f.close()
    out.close()
        
    
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

def main():
    cmdparser = argparse.ArgumentParser(description="Create new fastq files using reads NOT aligned to reference NON-CODING rna.",
                                        usage='%(prog)s ' ,prog='removeNonCoding.py'  )
    cmdparser.add_argument('-s', '--sam', action='store', dest='SAM', 
                           help='Sam file to filter for unaligned reads', metavar='')
    cmdparser.add_argument('-o', '--output', action='store', dest='OUT', 
                           help='Output fastq file, where reads are taken from input sam file.')
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
        
    print(inSame)
    print(outFile)
        
    
        
    
    
    
    if __name__ == "__main__":
    main()
