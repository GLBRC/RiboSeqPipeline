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
chromosome_info = { 'YPS1009':'/mnt/bigdata/linuxhome/mplace/data/reference/YPS1009/YPS1009_chromSizes.pkl',
              'S288C' : '/mnt/bigdata/linuxhome/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/individualChroms/S288C_chromSizes.pkl'
              }
# Pickled reference genome bed files, contains -72 bases upstream of start and +60 bases from stop. 
genome_info = { 'YPS1009' : '/home/glbrc.org/mplace/scripts/riboSeqPipeline/reference/YPS1009_GeneBed.pkl',
                'S288C' : '/home/glbrc.org/mplace/scripts/riboSeqPipeline/reference/S288C_GeneBed.pkl'  
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
        
    result = output[0].decode('utf-8')
    log    = output[1].decode('utf-8')

    # add logging here 
    print(result)
    #print('\n\n')
    print(log)

    
def countStarts(reference = 'YPS1009'):
    """countStarts
    
    Use interval tree to count alignment start sites.    
    """
    gffDB = open(genome_info[reference], 'rb')    
    geneBed = pickle.load(gffDB)
    for keys in geneBed:
        print(keys, '=>', geneBed[keys])
    gffDB.close()
    
    
    
'''
for aln in glob.glob('alignments/S288C/*-sorted.sam'):    
    #aln = 'alignments/S288C/ANEU-LOG-30-sorted.sam'
    outName = re.sub('-sorted.sam', '-RiboSeq-all-counts-S288C_v1.txt', aln)

    print('processing ', aln, outName)
    
    # set up results dictionary for the counting of start positions
    results = {}
    for chromName in chromBeds.keys():
        if chromName not in results:
            results[chromName] = {}   
    
    # open sam file , MAKE SURE THE SAM FILE IS SORTED (samtools sort)
    # we are interested in columns  2, 3, 4 (FLAG, RNAME, POS)
    # we use the NM tag to decide if there is an exact match
    # NM Edit distance to the reference, NM:i:0 (exact match), NM:i:1 is one mismatch.
    with open(aln, 'r') as f, open('gene_name.txt', 'w') as geneOut:
        for line in f:
            # skip header information
            if not line.startswith('@'):
                read = line.split('\t')
                if read[1] == '4':             # these reads are unmapped
                    continue
                elif read[1] == '16':           # mapped in reverse direction
                    chrom    = read[2]
                    readStart = int(read[3])
                    tags     = read[11:]       # section contains tags i.e. NM:i:0
                    for tag in tags:           # find NM tag
                        if tag.startswith('NM:'):
                            myTag = tag.split('NM:i:')[1]   
                            if myTag == '0' or myTag == '1':
                                cigar   = re.split('[A-Z]', read[5])   # split cigar string, removing Characters keeping numbers
                                cigar =  [i for i in cigar if i]      # remove any blank slots in list
                                seqLen = 0
                                for ix in cigar:                      # sum cigar lengths 
                                    seqLen += int(ix)  
                                readStart += seqLen       # minus strand, get sequence length to use the 5' end as readStart position
                                readStart -= 1            # subtract one from minus strand to get the read in the right frame 
                                gene = chromBeds[chrom].find_range([readStart, readStart])   # which gene does this read start in
                                if gene is not None:                                       
                                    for g in gene:                                         # this is a list, make sure we have the proper strand 
                                        if S288CBED[chrom][g]['strand'] == '-':
                                            #testout = f'{chrom}  {readStart}  {myTag}   {str(cigar)}  {str(seqLen)}  {gene} \n'
                                            #geneOut.write(testout)
                                            #out.write(f'{chrom},  {readStart}, {tag}, {myTag}, {g}\n')
                                            #print( 'start ' , S288CBED[chrom][g]['start'], 'end', S288CBED[chrom][g]['end'],
                                            #                    'strand', S288CBED[chrom][g]['strand'], 'pos', readStart, 'cnts', 1 )
                                            if g not in results[chrom]:
                                                results[chrom][g] = {'start' : S288CBED[chrom][g]['start'], 'end': S288CBED[chrom][g]['end'],
                                                                'strand': S288CBED[chrom][g]['strand'], 'pos' : {readStart : 1 }}
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
                                        if S288CBED[chrom][g]['strand'] == '+':
                                            #out.write(f'{chrom},  {readStart}, {tag}, {myTag}, {g}\n')
                                            #print( 'start ' , S288CBED[chrom][g]['start'], 'end', S288CBED[chrom][g]['end'],
                                            #                    'strand', S288CBED[chrom][g]['strand'], 'pos', readStart, 'cnts', 1 )
                                            if g not in results[chrom]:
                                                results[chrom][g] = {'start' : S288CBED[chrom][g]['start'], 'end': S288CBED[chrom][g]['end'],
                                                                'strand': S288CBED[chrom][g]['strand'], 'pos' : {readStart : 1 }}
                                            else:
                                                if readStart in results[chrom][g]['pos']:
                                                    results[chrom][g]['pos'][readStart] += 1 
                                                else:
                                                    results[chrom][g]['pos'][readStart] = 1   
                
    f.close()    
'''
'''  
    #results['ref|NC_001134|']['YBL039C']
    ### Write table of start position counts
    with open(outName, 'w') as out:
        #with open('MINUS-STRAND-TEST.txt', 'w') as out:
        for c in results.keys():
            for gene in results[c].keys():
                sortedPos = list(results[c][gene]['pos'].keys())
                startPos = results[c][gene]['start']
                endPos   = results[c][gene]['end']
                strand   = results[c][gene]['strand']
                sortedPos.sort()
                #print(c, gene, results[c][gene]['start'], results[c][gene]['end'], results[c][gene]['strand'])
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

                    #for i, tmp in enumerate(outPos, start=-72):
                    #    outCodon.append(str(i))
                    
                    #out.write(f'{c}\n')
                    #outPos.reverse()
                    #posLine   = f'Genome_Position ' + ' '.join(outPos) + '\n'
                    #out.write(posLine)
                    
                    #codonLine = f'{results[c][gene]["strand"]} ' + ' '.join(outCodon) + '\n'
                    #out.write(codonLine)
                    
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
                        
                    #out.write(f'{c}\n')
                    #posLine   = f'Genome_Position ' + ' '.join(outPos) + '\n'
                    #out.write(posLine)
                    #codonLine = f'{results[c][gene]["strand"]} ' + ' '.join(outCodon) + '\n'
                    #out.write(codonLine)
                    countLine = f'{gene} ' + ' '.join(outCnt)   + '\n'
                    out.write(countLine)
'''
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
    countStarts()
    
    

if __name__ == "__main__":
    main()