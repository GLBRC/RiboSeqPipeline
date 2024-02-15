#!/usr/bin/env python
"""CountRiboStarts.py

Given a sam file produce a table of read start counts.
Part of the Ribo-Seq pipeline.

Notes
-----
Currently only set up to use YPS1009.
   
Parameters
----------
s : str
    Sam alignment file, assumed to be aligned to YPS1009(default) or S288C

Example
-------
    usage:

        CountRiboStarts.py -s sample_X.sam
    

"""
import argparse 
import os
import pickle 
import re
import subprocess
import sys
from interval_tree import IntervalTree  # for efficient searching of gene locations

# Pickled reference genome chromosome sizes
chromosome_info = { 'YPS1009':'/home/glbrc.org/mplace/scripts/riboSeqPipeline/reference/YPS1009_chromSizes.pkl',
              'S288C' : '/home/glbrc.org/mplace/scripts/riboSeqPipeline/reference/S288C_chromSizes.pkl'
              }
# reference genome bed files, contains -72 bases upstream of start and +60 bases from stop. 
genome_info = { 'YPS1009' : '/home/glbrc.org/mplace/scripts/riboSeqPipeline/reference/YPS1009-Genes.bed',
                'S288C' : '/home/glbrc.org/mplace/scripts/riboSeqPipeline/reference/S288C-Genes.bed'  
              }

def chrSize():
    """chrSize

    Create a dictionary of chromosome sizes, required for the interval tree.
    """
    #chromosome sizes, Need these values for the interval tree
    chromSize = {}
    with open(chromosome_info['YPS1009'], 'r') as f:
        for ln in f:
            chrom, length = ln.rstrip().split()
            if chrom not in chromSize:
                chromSize[chrom] = int(length)

    return chromSize

def sortSam(samFile, outSam ):
    """sortSam
    
    Sort input alignment file.    
    """
    # set up samtools sort call, output to bam
    cmd = ['samtools', 'sort', '-O', 'SAM', '-o', outSam, samFile]
    output = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        
    #result = output[0].decode('utf-8')
    #log    = output[1].decode('utf-8') 
    
def countStarts(sortedSam, reference = 'YPS1009'):
    """countStarts
    
    Use interval tree to count alignment start sites.    
    """
    # load the reference genome chromosize file
    sizeDB = open(chromosome_info[reference], 'rb')
    chromSizes = pickle.load(sizeDB)
    sizeDB.close()   
    
    # load the YPS1009 modified bed file
    BED = {}
    with open(genome_info[reference],'r') as f:
        f.readline()
        for ln in f:
            dat = ln.rstrip().split()
            if dat[0] not in BED:
                BED[dat[0]] = {}
                if dat[3] not in BED[dat[0]]:
                    BED[dat[0]][dat[3]] = {'start': int(dat[1]), 'end': int(dat[2]), 'strand':dat[4]}
            else:
                if dat[3] not in BED[dat[0]]:
                    BED[dat[0]][dat[3]] = {'start': int(dat[1]), 'end': int(dat[2]), 'strand':dat[4]}
    
    # set first reference genome chromosome
    if reference == 'YPS1009':
        prevChrom = 'chrI'
    else:
        prevChrom = 'ref|NC_001133|'  # for S288C R64-1-1
    
    # construct interval tree 
    chromBeds = {}
    with open(genome_info[reference], 'r') as f:
        f.readline()                        # skip header
        features = []
        for ln in f:                         
            dat = ln.split()                #  ['ref|NC_001133|', '263', '709', 'YAL069W', '+']
            if dat[0] == prevChrom:
                row = [int(dat[1]), int(dat[2]), dat[3]]
                features.append(row)
            elif dat[0] != prevChrom:
                chromBeds[prevChrom] = IntervalTree(features,1, int(chromSizes[prevChrom])) # write most recent
                if dat[0] not in chromBeds:
                    chromBeds[dat[0]] = None
                features = []
                row = [int(dat[1]), int(dat[2]), dat[3]]
                features.append(row)
                prevChrom = dat[0]
        
        chromBeds[prevChrom] = IntervalTree(features,1, int(chromSizes[prevChrom])) # write most recent
        
    outName = re.sub('-sorted.sam', f'-RiboSeq-all-counts-{reference}_v1.txt', sortedSam)

    print('processing ', sortedSam, outName)
    
    # set up results dictionary for the counting of start positions
    results = {}
    for chromName in chromBeds.keys():
        if chromName not in results:
            results[chromName] = {}   
    
    # open sorted sam file 
    # we are interested in columns  2, 3, 4 (FLAG, RNAME, POS)
    # we use the NM tag to decide if there is an exact match
    # NM Edit distance to the reference, NM:i:0 (exact match), NM:i:1 is one mismatch.
    # Default is to allow 1 mismatch.
    with open(sortedSam, 'r') as f:
        for line in f:
            # skip header information
            if not line.startswith('@'):
                read = line.split('\t')
                if read[1] == '4':             # these reads are unmapped
                    continue
                elif read[1] == '16':          # mapped in reverse direction
                    chrom    = read[2]
                    readStart = int(read[3])
                    tags     = read[11:]       # section contains tags i.e. NM:i:0
                    for tag in tags:           # find NM tag
                        if tag.startswith('NM:'):
                            myTag = tag.split('NM:i:')[1]   
                            if myTag == '0' or myTag == '1':
                                cigar   = re.split('[A-Z]', read[5])  # split cigar string, removing Characters keeping numbers
                                cigar =  [i for i in cigar if i]      # remove any blank slots in list
                                seqLen = 0
                                for ix in cigar:                      # sum cigar lengths 
                                    seqLen += int(ix)  
                                readStart += seqLen       # minus strand, get sequence length to use the 5' end as readStart position
                                readStart -= 1            # subtract one from minus strand to get the read in the right frame 
                                gene = chromBeds[chrom].find_range([readStart, readStart])   # which gene does this read start in
                                if gene is not None:                                       
                                    for g in gene:                                           # this is a list, make sure we have the proper strand 
                                        if BED[chrom][g]['strand'] == '-':
                                            if g not in results[chrom]:
                                                results[chrom][g] = {'start' : BED[chrom][g]['start'], 'end': BED[chrom][g]['end'],
                                                                'strand': BED[chrom][g]['strand'], 'pos' : {readStart : 1 }}
                                            else:
                                                if readStart in results[chrom][g]['pos']:
                                                    results[chrom][g]['pos'][readStart] += 1 
                                                else:
                                                    results[chrom][g]['pos'][readStart] = 1 
                elif read[1] == '0':           # mapped in forward direction
                    chrom    = read[2]
                    readStart = int(read[3])
                    tags     = read[11:]       # section contains NM tag
                    for tag in tags:           # find NM tag
                        if tag.startswith('NM:'):
                            myTag = tag.split('NM:i:')[1]   
                            if myTag == '0' or myTag == '1':
                                gene = chromBeds[chrom].find_range([readStart, readStart])   # which gene does this read start in
                                if gene is not None:                                       
                                    for g in gene:                                         # this is a list, make sure we have the proper strand 
                                        if BED[chrom][g]['strand'] == '+':
                                            if g not in results[chrom]:
                                                results[chrom][g] = {'start' : BED[chrom][g]['start'], 'end': BED[chrom][g]['end'],
                                                                'strand': BED[chrom][g]['strand'], 'pos' : {readStart : 1 }}
                                            else:
                                                if readStart in results[chrom][g]['pos']:
                                                    results[chrom][g]['pos'][readStart] += 1 
                                                else:
                                                    results[chrom][g]['pos'][readStart] = 1   
    f.close()    

    ### Write table of start position counts
    with open(outName, 'w') as out:
        for c in results.keys():      # c = chromosome
            for gene in results[c].keys():
                sortedPos = list(results[c][gene]['pos'].keys())
                startPos = results[c][gene]['start']
                endPos   = results[c][gene]['end']
                strand   = results[c][gene]['strand']
                sortedPos.sort()
                outCodon = []         
                outPos = []           # genome position
                outCnt = []
            
                if strand == "-":
                    # create the codon positions starting at -72
                    for p in range(results[c][gene]['start'], results[c][gene]['end'] + 1):
                        outPos.append(str(p))
                        if p in results[c][gene]['pos']:
                            outCnt.append(str(results[c][gene]['pos'][p]))
                        else:
                            outCnt.append('0')
                    
                    outCnt.reverse()
                    countLine = f'{gene} ' + ' '.join(outCnt)   + '\n'
                    out.write(countLine)
                elif strand == "+":
                    codonPos = -72  
                    for p in range(results[c][gene]['start'], results[c][gene]['end']): 
                        outCodon.append(str(codonPos))
                        codonPos +=1
                        outPos.append(str(p))
                        if p in results[c][gene]['pos']:
                            outCnt.append(str(results[c][gene]['pos'][p]))
                        else:
                            outCnt.append('0')
                        
                    countLine = f'{gene} ' + ' '.join(outCnt)   + '\n'
                    out.write(countLine)

def main():
    
    cmdparser = argparse.ArgumentParser(description="Count Ribo-Seq start sites using reference genome aligned sam file.",
                                        usage='%(prog)s -s <alignment.sam>' ,prog='CountRiboStars.py'  )
    cmdparser.add_argument('-s', '--sam', action='store', dest='SAM', help='Sample reads aligned to reference genome, SAM format.', metavar='')
    cmdparser.add_argument('-o', '--out', action='store', dest='OUT', help='Sorted sam file output name.', metavar='')
    cmdResults = vars(cmdparser.parse_args())
        
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
    # process command line arguments
    if cmdResults['SAM'] is not None:
        samFile = cmdResults['SAM']
    else:
        print()
        cmdparser.print_help()
        sys.exit(1)
        
    if cmdResults['OUT'] is not None:
        outSam = cmdResults['OUT']
    else:
        print()
        cmdparser.print_help()
        sys.exit(1)
        
    sortSam(samFile, outSam)
    countStarts(outSam)
    
    

if __name__ == "__main__":
    main()