#!/usr/bin/env python
"""RiboSeqPipeline.py

Pipeline to align and count start sites for Ribo-Seq data. Pipeline utilizes HTCondor Dagman 
to manage jobs.

Notes
----- 
Put the fastq files in the parent directory within a directory called fastq.
If data is paired-end only Read 1 is required.
 
    ex:   MyProject
            fastq
                data_A_R1.fastq
                data_B_R1.fastq
                data_C_R1.fastq

use the riboSeq environment: /home/glbrc.org/mplace/.conda/envs/riboSeq 

RiboSeq processing pipeline Steps:

    1) Fastqc to check the read length distribution. 

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
import subprocess
import sys
from pydagman.dagfile import Dagfile
from pydagman.job import Job
from interval_tree import IntervalTree


parentDir = os.getcwd() + '/'
resourceDir = '/home/glbrc.org/mplace/scripts/riboSeqPipeline/reference/'
# bowtie2 reference for genomes
S288C_RefGenome = '/home/glbrc.org/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/s.cerevisiae-R64-1-1' 
YPS1009_RefGenome = '/home/glbrc.org/mplace/data/reference/YPS1009/YPS1009'

def runFastqc():
    """runFastqc
    
    write the fastqc condor submit and shell script files.
    """
    with open('fastqc.jtf', 'w') as submit:
                submit.write( "Universe                 = vanilla\n" )
                submit.write( "Executable               = runFastqc.sh\n")
                submit.write( "Arguments                = $(fastqFile)\n")
                submit.write( "Error                    = fastqc.submit.err\n")
                submit.write( "Log                      = fastqc.submit.log\n")  
                submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
                submit.write( "Queue\n" )
    submit.close()  

    # write shell script to run fastqc
    with open('runFastqc.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write("fastqc $1\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runFastqc.sh', 0o0777)

def runCutadapt():
    """runCutadapt
    
    Run cutadapt on fastq files using the parameters provided by Ezra Bio.
    -j 8 
    -g "^GGG" 
    -a "A{10}" 
    -n 2 
    -m 15 
    --max-n=0.1 
    --discard-casava 
    -o output.fastq.gz input.fastq.gz" .
    """
    # setup input file for cutadapt
    cutadaptOutDir = parentDir + 'cutadapt/'
    if os.path.exists(cutadaptOutDir):
        print("Directory exists.")
    else:
        os.mkdir(cutadaptOutDir)

    with open('inputFastq.txt', 'r') as f, open('cutadaptInput.txt', 'w') as out:
        for fstq in f:
            fstqName = re.sub('.fastq', '-clean.fastq', os.path.basename(fstq.rstrip()))
            fstqOutName = cutadaptOutDir + fstqName
            out.write(f'{fstq.rstrip()} {fstqOutName}\n')
    f.close()
    out.close()
    
        # write the cutadapt condor submit file
    with open('cutadapt.submit', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = runCutAdapt.sh\n")
        submit.write( "Arguments                = $(fastqFile) $(outFastq)\n")
        submit.write( "Error                    = cutadapt.submit.err\n")
        submit.write( "Log                      = cutadapt.submit.log\n")  
        submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
        submit.write( "Queue fastqFile, outFastq from cutadaptInput.txt\n" )
    submit.close()

    # write shell script to run cutadapt
    with open('runCutAdapt.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write("cutadapt -j 8 -g ^GGG -a A{10} -n 2 -m 15 --max-n=0.1 --discard-casava -o $2 $1\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runCutAdapt.sh', 0o0777)

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

def alignNonCoding():
    """alignNonCoding
    
    """
    # Step 4 ) Align reads with Bowtie2 to non-coding RNA, 
    # https://downloads.yeastgenome.org/sequence/S288C_reference/rna/archive/rna_coding_R64-1-1_20110203.fasta.gz 
    # reads that align will be discarded.  Allow 1 mismatch in bowtie2 alignment.
    # bowtie2 -p 8 --phred33 -N 1 -x $REFERENCE -U $file -S $out.sam 
    # -p number of threads
    # -N Sets the number of mismatches
    # -x The basename of the index for the reference genome
    # -U file to align (unpaired)
    # -S File to write SAM alignments to

    # setup input file for bowtie2 alignment to non-coding RNA
    nonCodingOutDir = parentDir + 'alignNonCodingRNA/'
    if os.path.exists(nonCodingOutDir):
        print("Directory exists.")
    else:
            os.mkdir(nonCodingOutDir)

    # get a list of cutadapt cleaned & filtered fastq files for alignment
    with open('alignmentInput.txt', 'w') as out:
        for cleanfstq in glob.glob(cutadaptOutDir + '*-filt.fastq'):
            samFile = re.sub('cutadapt', 'alignNonCodingRNA', re.sub('-filt.fastq', '.sam', cleanfstq))        
            out.write(cleanfstq + ' ' + samFile + '\n')
    out.close()
    
        # write the bowtie2 condor submit file
    with open('ncbowtie2.submit', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = runncBowtie2.sh\n")
        submit.write( "Arguments                = $(fastqFile) $(sam)\n")
        submit.write( "Error                    = ncbowtie2.submit.err\n")
        submit.write( "Log                      = ncbowtie2.submit.log\n")  
        submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
        submit.write( "Queue fastqFile, sam from alignmentInput.txt\n" )
    submit.close()

    # write shell script to run cutadapt
    with open('runncBowtie2.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write("bowtie2 -p 8 --phred33 -N 1 -x /mnt/bigdata/linuxhome/mplace/scripts/riboseq/reference/rna_coding_R64-1-1 -U $1 -S $2\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runncBowtie2.sh', 0o0777)

def removeNonCoding():
    """removeNonCoding
    -F --exclude-flags FLAG ,     Do not output alignments with any bits set in FLAG present in the FLAG field.
    in this case we exclude the reads which did not align to the Non-coding RNA i.e. the reads we want to analyze later.
    """
    # get a list of reads which aligned to Non-Coding RNA and remove them from the clean.fastq files
    #samtools view -F 4 -u  -O SAM -o mapped.sam TestSample1.sam
    with open('filterSamInput.txt', 'w') as out:
        for sam in glob.glob(nonCodingOutDir + '*.sam'):
            sampleName = re.sub('.sam', '', os.path.basename(sam))
            outSam = nonCodingOutDir + sampleName + '-unmapped.sam'
            out.write(sam + ' ' + outSam + '\n')
    out.close()
    
    # write the samtools filter UNMAPPED (reads which did not align to Non-Coding RNA) reads condor submit file
    # these reads will be aligned to S288C and YPS1009
    with open('filter.submit', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = runfilter.sh\n")
        submit.write( "Arguments                = $(sam) $(ncsam)\n")
        submit.write( "Error                    = filter.submit.err\n")
        submit.write( "Log                      = filter.submit.log\n")  
        submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
        submit.write( "Queue sam, ncsam from filterSamInput.txt\n" )
    submit.close()

    # write shell script to run 
    with open('runfilter.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write("samtools view -f 4 -u -O SAM -o $2 $1\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runfilter.sh', 0o0777)

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

def alignBowtie():
    """alignBowtie
    
    Align reads to reference genomes S288C and YPS1009
    align reads to S288C using bowtie2 
    write the bowtie2 condor submit file
    """    
    # Create input file for alignment to the S288C reference genome.
    outputPath = parentDir + 'alignments/S288C/'
    with open('refAlignmentS288C_input.txt', 'w') as out:
        for inFastq in glob.glob(parentDir + 'alignments/*.fastq'):
            sampleName = re.sub('.fastq', '.sam', os.path.basename(inFastq))
            out.write(inFastq + ' ' + outputPath + sampleName + ' /home/glbrc.org/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/s.cerevisiae-R64-1-1' + '\n')
    out.close

    #S288C
    with open('s288cbowtie2.submit', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = runBowtie2.sh\n")
        submit.write( "Arguments                = $(fastqFile) $(sam) $(ref)\n")
        submit.write( "Error                    = s288c_bowtie2.submit.err\n")
        submit.write( "Log                      = s288c_bowtie2.submit.log\n")  
        submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
        submit.write( "Queue fastqFile, sam, ref from refAlignmentS288C_input.txt\n" )
    submit.close()
    #YPS1009
    with open('yps1009bowtie2.submit', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = runBowtie2.sh\n")
        submit.write( "Arguments                = $(fastqFile) $(sam) $(ref)\n")
        submit.write( "Error                    = yps1009_bowtie2.submit.err\n")
        submit.write( "Log                      = yps1009_bowtie2.submit.log\n")  
        submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
        submit.write( "Queue fastqFile, sam, ref from refAlignmentYPS1009_input.txt\n" )
    submit.close()

    # write shell script to run bowtie2
    with open('runBowtie2.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write(f"bowtie2 -p 8 --phred33 -N 1 -x $3 -U $1 -S $2\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runBowtie2.sh', 0o0777)
    
def sortAlignment():
    """sortAlignment
    
    # samtools sort ANEU-LOG-30_S2_R1_001.sam > ANEU-LOG-30_S2_R1_001-sorted.sam
    # Create input file for samtools sort 
    """
    inputFilePath = parentDir + 'alignments/YPS1009/'
    print(inputFilePath)
    with open('sort-YPS1009_input.txt', 'w') as out:
        for sam in glob.glob(inputFilePath + '*.sam'):
            sampleName = re.sub('.sam', '.bam', os.path.basename(sam))
            out.write(sam + ' ' + inputFilePath + sampleName + ' \n')
    out.close

    # Sort the alignment output files S288C and YPS1009
    with open('sort.submit', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = runsort.sh\n")
        submit.write( "Arguments                = $(sam) $(bam)\n")
        submit.write( "Error                    = sort.submit.err\n")
        submit.write( "Log                      = sort.submit.log\n")  
        submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
        #submit.write( "Queue sam, bam from  sort-S288C_input.txt\n" )
        submit.write( "Queue sam, bam from sort-YPS1009_input.txt\n" )
    submit.close()

    # write shell script to run sorting
    with open('runsort.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write("samtools sort $1  > $2\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runsort.sh', 0o0777)

def prepIntervalTree():
    """prepIntervalTree
    
    """        
    prevChrom = 'ref|NC_001133|'
    chromBeds = {}
    with open('S288C-Genes.bed', 'r') as f:
        f.readline()                        # skip header
        features = []
        for ln in f:
            dat = ln.split()                #  ref|NC_001133| 263 709 YAL069W +

            if dat[0] == prevChrom:
                row = [int(dat[1]), int(dat[2]), dat[3]]
                features.append(row)
            elif dat[0] != prevChrom:
                chromBeds[prevChrom] = IntervalTree(features,1, int(S288C_chromSizes[prevChrom])) # write most recent
                if dat[0] not in chromBeds:
                    chromBeds[dat[0]] = None
                features = []
                row = [int(dat[1]), int(dat[2]), dat[3]]
                features.append(row)
                prevChrom = dat[0]
        
        chromBeds[prevChrom] = IntervalTree(features,1, int(S288C_chromSizes[prevChrom])) # write most recent

def countStartPos():
    pass

def cleanUp():
    pass

def writeReadMe():
    pass

def writeCitations():
    
    
    pass

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

    # create Dagfile object, utilize HTCondor Dagman to manage pipeline.
    mydag = Dagfile()
    num = 1

    # Set up Dagman jobs for each input fastq file
    for fsa in fastqLst:
        # 1st step run fastqc
        fastqcJob = Job('fastqc.jtf', 'job' + str(num))   
        fastqcJob.pre_skip(1)
        fastqcJob.add_var('fastq', 'fastq/' + fsa)
        mydag.add_job(fastqcJob)
        num += 1

    
    
    

    mydag.save('MasterDagman.dsf')      # write the dag submit file
    
    # write the required condor files
    runFastqc()
                
        
    


if __name__ == "__main__":
    main()

