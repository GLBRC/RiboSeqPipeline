#!/usr/bin/env python
"""CreateGeneBed.py

Pipeline to align and count start sites for Ribo-Seq data.

Notes
----- 
Put the fastq files in the parent directory within a parent directory called fastq.

use the riboSeq environment: /home/glbrc.org/mplace/.conda/envs/riboSeq 

 RiboSeq processing pipeline

Steps

Run Fastqc to check the read length distribution. 

Run cutadapt using the parameters provided by Ezrabio.  
"-j 8 -g "^GGG" -a "A{10}" -n 2 -m 15 --max-n=0.1 --discard-casava -o output.fastq.gz input.fastq.gz"

Remove reads where the first position quality score is <=10

Align reads with Bowtie to non-coding RNA, 
https://downloas.yeastgenome.org/sequence/S288C_reference/rna/archive/rna_coding_R64-1-1_20110203.fasta.gz 
reads that align will be discarded. 

Allow 1 mismatch in bowtie alignment.

Alignment files (SAM/BAM) provide forward/reverse alignment information is provided by the Sam Flag.    0 is read aligned in the Fwd direction
    4 is unaligned
    16 is aligned in the Rvs direction

Align the remaining reads with Bowtie to YPS1009 and S288C reference genomes.Run samtools mpileup to generate base counts for all genes in YPS1009 and S288C.

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
import pickle 
import re
import sys


def createInfoDictionary():
    """createInfoDictionary
    """
    Create information dictionaries for S288C
    gff = '/mnt/bigdata/linuxhome/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff'
    geneLookUp = {}    # dictionary of dictionaries key = chrom : {key = the -72 position of gene value = gene name }
    gffDict = {}       # dictionary of dictionaries key = chrom : {key = {start : 0, end : 0, strand : '+'}
    
    with open(gff, 'r') as g:
        for line in g:
            if line.startswith('#'):   # skip comment rows
                continue
            else:
                dat = line.split('\t')
                chrom = dat[0]
                if chrom not in S288C_chromSizes:
                    continue
                if chrom not in geneLookUp:
                    geneLookUp[chrom] = {}
                # add chrom to gffDict
                if chrom not in gffDict:
                    gffDict[chrom] = {}
                # identify genes
                if dat[2] == 'gene':
                    geneName = re.sub('ID=', '', dat[8].split(';')[0])
                    # identify strand
                    if dat[6] == '+':          # POSIIVE STRAND
                        if int(dat[3]) < 72:          
                            start = 0
                        else:
                            start = int(dat[3]) - 72
                        if start not in geneLookUp:
                            geneLookUp[chrom][start] = geneName
                        else:
                            print('Duplicate start ', geneName, geneLookUp[chrom][start])

                        if not int(dat[4]) + 60 > S288C_chromSizes[chrom]:
                            end   = int(dat[4]) + 60
                        else:
                            end = S288C_chromSizes[chrom]

                        # add final gene info for positive strand 
                        if geneName not in gffDict[chrom]:
                            gffDict[chrom][geneName] = {'strand' : '+', 'start': start, 'end': end }
                    else:                      # NEGATIVE STRAND
                        if  int(dat[3]) - 60 < 60:
                            start = 0
                        else:
                            start = int(dat[3]) - 60
                        if start not in geneLookUp:
                            geneLookUp[chrom][start] = geneName
                        else:
                            print('Duplicate start minus strand ', geneName, geneLookUp[chrom][start]) 

                        if not int(dat[4]) + 72 > S288C_chromSizes[chrom]:
                            end = int(dat[4]) + 72
                        else:
                            end = S288C_chromSizes[chrom]
                        # add final gene info for positive strand 
                        if geneName not in gffDict[chrom]:
                            gffDict[chrom][geneName] = {'strand' : '-', 'start': start, 'end': end }

    # write bed file for use with Samtools  -l flag, start and end positions have been adjust -72 from start and + 60 to end
    # track name=yeast_genes description="All S.cerevisiae Gene locations"
    with open('S288C-Genes.bed', 'w') as out:
        out.write('track name=yeast_genes description="All S.cerevisiae Gene locations"\n')
        for chrom in gffDict.keys():
            for gene, dat in gffDict[chrom].items():
                out.write(f'{chrom} {gffDict[chrom][gene]["start"]} {gffDict[chrom][gene]["end"]} {gene} {gffDict[chrom][gene]["strand"]}\n')

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
               
    # create output files for genome alignments
    if not os.path.exists(parentDir + 'alignments'):
        os.mkdir(parentDir + 'alignments')
        os.mkdir(parentDir + 'alignments/S288C')
        os.mkdir(parentDir + 'alignments/YPS1009') 

    
    
                
        
    


if __name__ == "__main__":
    main()

