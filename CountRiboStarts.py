#!/usr/bin/env python
"""RiboSeqPipeline.py
"""

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