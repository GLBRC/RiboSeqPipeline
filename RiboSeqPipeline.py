#!/usr/bin/env python
"""RiboSeqPipeline.py

Pipeline to align and count start sites for Ribo-Seq data. 
Pipeline utilizes HTCondor DAGman to manage jobs.

Notes
----- 
Put the fastq files in the parent directory within a directory called fastq.
If data is paired-end only Read 1 is required, do not include Read 2 files.
 
    ex:   MyProject
            fastq
                data_A_R1.fastq
                data_B_R1.fastq
                data_C_R1.fastq

Alignment files (SAM/BAM) provide forward/reverse alignment information is provided by the Sam Flag.   
    0 is read aligned in the Fwd direction
    4 is read unaligned
    16 is read aligned in the Rvs direction

Must run on scarcity-submit.glbrc.org.
use the riboSeq environment: /home/glbrc.org/mplace/.conda/envs/riboSeq 

Method
------
RiboSeq processing pipeline Steps:

    1) Fastqc to check the read length distribution. 

    2) Run cutadapt using the parameters provided by Ezrabio.  
    
    "-j 8 -g "^GGG" -a "A{10}" -n 2 -m 15 --max-n=0.1 --discard-casava -o output.fastq.gz input.fastq.gz"

    3) Remove reads where the first position quality score is <=10

    4) Align reads with Bowtie to non-coding RNA, reads that align will be discarded. 
        https://downloas.yeastgenome.org/sequence/S288C_reference/rna/archive/rna_coding_R64-1-1_20110203.fasta.gz 
        
    5) Align remaining reads to reference genome allowing 1 mismatch in the bowtie2 alignment.
    
Parameters
----------
f : str
    A text file with a list of fastq files to process. Only the forward read is 
    required if you have paired-end data.

Example
-------
    usage:

        RiboSeqPipeline.py -f fastqfiles.txt
    
Requirements
------------

    1. HTCondor, specifically GLBRC scarcity compute cluster.

"""
import argparse 
import os
import re
import subprocess
import sys
from pydagman.dagfile import Dagfile
from pydagman.job import Job
from interval_tree import IntervalTree

parentDir = os.getcwd() + '/'
resourceDir = '/home/glbrc.org/mplace/scripts/riboSeqPipeline/reference/'
# bowtie2 reference for genomes
refGenomes = { 'S288C' : '/home/glbrc.org/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/s.cerevisiae-R64-1-1' ,
              'YPS1009' : '/home/glbrc.org/mplace/data/reference/YPS1009/YPS1009' }

def runFastqc():
    """runFastqc
    
    write the fastqc condor submit and shell script files.
    """
    with open('fastqc.jtf', 'w') as submit:
                submit.write( "Universe                 = vanilla\n" )
                submit.write( "Executable               = runFastqc.sh\n")
                submit.write( "Arguments                = $(fastqFile)\n")
                submit.write( "Error                    = log/fastqc.$(job).err\n")
                submit.write( "Log                      = log/fastqc.$(job).log\n")  
                submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
                submit.write( "Queue\n" )
    submit.close()  

    # write shell script to run fastqc
    with open('runFastqc.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write("/home/glbrc.org/mplace/scripts/riboSeqPipeline/runFastqc.py -f $1\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runFastqc.sh', 0o0777)

def runCutAdapt():
    """runCutadapt
    
    Run cutadapt on fastq files using the parameters provided by EzraBio.
    -j 8                - number of cores
    -g "^GGG"           - Sequence of an adapter ligated to the 5' end
    -a "A{10}"          - Sequence of an adapter ligated to the 3' end
    -n 2                - Remove up to "x" adapters from each read, 2 in our case
    -m 15               - Discard reads shorter than 15 bp
    --max-n=0.1         - Discard reads with more than COUNT 'N' bases, 0.1 is % of read length
    --discard-casava    - Discard reads that did not pass CASAVA filtering
    -o output.fastq.gz input.fastq.gz" - output name
    """    
    # write the cutadapt condor submit file
    with open('cutadapt.jtf', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = runCutAdapt.sh\n")
        submit.write( "Arguments                = $(infastq) $(outfastq)\n")
        submit.write( "Error                    = log/cutadapt-$(job).submit.err\n")
        submit.write( "Log                      = log/cutadapt-$(job).submit.log\n")  
        submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
        submit.write( "Queue\n" )
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

def runRemoveFirst():
    """runRemoveFirst
    
    write the condor submit and shell script files to run removeFirstBase.py.
    """
    with open('rmFirst.jtf', 'w') as submit:
                submit.write( "Universe                 = vanilla\n" )
                submit.write( "Executable               = runRmFirst.sh\n")
                submit.write( "Arguments                = $(cleanFastq)\n")
                submit.write( "Error                    = log/rmFirst.$(job).err\n")
                submit.write( "Log                      = log/rmFirst.$(job).log\n")  
                submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
                submit.write( "Queue\n" )
    submit.close()  

    # write shell script to run fastqc
    with open('runRmFirst.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write("/home/glbrc.org/mplace/scripts/riboSeqPipeline/removeFirstBase.py -f $1\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runRmFirst.sh', 0o0777)

def alignNonCoding():
    """alignNonCoding
    
    Write the shell script & condor submit file to align reads with Bowtie2 to 
    non-coding RNA. Non-coding RNA reference fasta is from SGD genome version
    R64-1-1, reads which align will be discarded.  Allow 1 mismatch 
    in bowtie2 alignment. 
     
    Command used:
     
    bowtie2 -p 8 --phred33 -N 1 -x $REFERENCE -U $file -S $out.sam 
    
    -p number of threads
    -N Sets the number of mismatches
    -x The basename of the index for the reference genome
    -U file to align (unpaired)
    -S File to write SAM alignments
    """
    # write the bowtie2 condor submit file
    with open('ncbowtie2.jtf', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = runncBowtie2.sh\n")
        submit.write( "Arguments                = $(ncFastq) $(sam) $(unaligned)\n")
        submit.write( "Error                    = log/ncbowtie2.$(job).err\n")
        submit.write( "Log                      = log/ncbowtie2.$(job).log\n")  
        submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
        submit.write( "Queue\n" )
    submit.close()

    # write shell script to run bowtie2
    with open('runncBowtie2.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write("bowtie2 -p 8 --phred33 -N 1 -x /mnt/bigdata/linuxhome/mplace/scripts/riboSeqPipeline/reference/rna_coding_R64-1-1 -U $1 -S $2 --un $3\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runncBowtie2.sh', 0o0777)

def alignBowtie():
    """alignBowtie
    
    Write the bowtie2 shell script & condor submit file for reference
    genome alignment.
    """    
    # generic submit file
    with open('alignment.jtf', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = runAlignment.sh\n")
        submit.write( "Arguments                = $(fastqFile) $(sam) $(ref)\n")
        submit.write( "Error                    = log/alignment.$(job).err\n")
        submit.write( "Log                      = log/alignment.$(job).log\n")  
        submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
        submit.write( "Queue\n" )
    submit.close()

    # write shell script to run bowtie2
    with open('runAlignment.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write(f"bowtie2 -p 8 --phred33 -N 1 -x $3 -U $1 -S $2\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runAlignment.sh', 0o0777)
    
def countStarts():
    """countStarts
    
    Write the Count alignment start positions shell script & condor submit file. 
    """
    with open('counts.jtf', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = runCounts.sh\n")
        submit.write( "Arguments                = $(sam) $(out)\n")
        submit.write( "Error                    = log/counts.$(job).err\n")
        submit.write( "Log                      = log/counts.$(job).log\n")  
        submit.write( "Requirements             = OpSysandVer == \"CentOS7\"\n")
        submit.write( "Queue\n" )
    submit.close()

    # write shell script to run CountRiboStarts.py
    with open('runCounts.sh', 'w') as out:
        out.write("#!/bin/bash\n")
        out.write("source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh\n")
        out.write("unset PYTHONPATH\n")  
        out.write("conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq\n")
        out.write(f"/home/glbrc.org/mplace/scripts/riboSeqPipeline/CountRiboStarts.py -s $1 -o $2\n")
        out.write("conda deactivate")
    out.close()

    os.chmod('runCounts.sh', 0o0777)
    
def cleanUp():
    pass

def writeReadMe():
    pass

def writeCitations():
    pass

def main():
    
    cmdparser = argparse.ArgumentParser(description="Ribo-Seq pipeline, produces alignments and counts.",
                                        usage='%(prog)s -f <fastqFileList.txt>' ,prog='RiboSeqPipeline.py'  )
    cmdparser.add_argument('-f', '--file',    action='store', dest='FILE',
                           help='Text file, one fastq file name per line, read 1 only if paired-end.', metavar='')
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
                
    # setup output dirs for cutadapt
    cutadaptOutDir = parentDir + 'cutadapt/'
    if not os.path.exists(cutadaptOutDir):
        os.mkdir(cutadaptOutDir)
               
    logDir = parentDir + 'log/'
    if not os.path.exists(logDir):
        os.mkdir(logDir)     

    # setup input file for bowtie2 alignment to non-coding RNA
    nonCodingOutDir = parentDir + 'alignNonCodingRNA/'
    if not os.path.exists(nonCodingOutDir):
        os.mkdir(nonCodingOutDir)

    # create output dirs for genome alignments
    if not os.path.exists(parentDir + 'alignments'):
        os.mkdir(parentDir + 'alignments')

    # create Dagfile object, utilize HTCondor Dagman to manage pipeline.
    mydag = Dagfile()
    num = 1

    # Set up Dagman jobs for each input fastq file
    for fsa in fastqLst:
        # 1st step run fastqc
        fastqcJob = Job('fastqc.jtf', 'job' + str(num))   
        fastqcJob.pre_skip(1)
        fastqcJob.add_var('fastqFile', 'fastq/' + fsa)
        fastqcJob.add_var('job', 'job' +  str(num))
        mydag.add_job(fastqcJob)
        num += 1
        
        # 2nd step run cutadapt
        cutadaptOutName = cutadaptOutDir + re.sub('.fastq', '-clean.fastq', os.path.basename(fsa.rstrip()))
        cutadaptJob = Job('cutadapt.jtf', 'job' + str(num))
        cutadaptJob.pre_skip(1)
        cutadaptJob.add_var('infastq', 'fastq/' + fsa)
        cutadaptJob.add_var('outfastq', cutadaptOutName)
        cutadaptJob.add_parent(fastqcJob)
        mydag.add_job(cutadaptJob)
        num += 1    
        
        # 3rd step run removeFirstBase
        rmJob = Job('rmFirst.jtf', 'job' + str(num))
        rmJob.pre_skip(1)
        rmJob.add_var('cleanFastq', cutadaptOutName)
        rmJob.add_parent(cutadaptJob)
        mydag.add_job(rmJob)
        num += 1 
        
        # 4th step run align reads to Non-coding RNA, outputs reads aligned to 
        # Non-coding RNA and a fastq file for reads that did not align.
        ncInputFastq = cutadaptOutDir + re.sub('.fastq', '-filt.fastq', fsa)
        ncSamOut     = nonCodingOutDir + re.sub('.fastq', '.sam', fsa)
        unalignedFastq = nonCodingOutDir + re.sub('.fastq', '-unaligned.fastq', fsa )
        alignNCJob = Job('ncbowtie2.jtf', 'job' + str(num))
        alignNCJob.pre_skip(1)
        alignNCJob.add_var('ncFastq', ncInputFastq)
        alignNCJob.add_var('sam', ncSamOut)
        alignNCJob.add_var('unaligned', unalignedFastq)
        alignNCJob.add_var('job', str(num))
        alignNCJob.add_parent(rmJob)
        mydag.add_job(alignNCJob)
        num += 1 
        
        # 5th step align unaligned reads from previous step to reference genome.
        alignedSam = parentDir + 'alignments/' + re.sub('.fastq', '.sam', fsa)
        alignJob   = Job('alignment.jtf', 'job' + str(num))
        alignJob.pre_skip(1)
        alignJob.add_var('fastqFile',unalignedFastq )
        alignJob.add_var('sam', alignedSam) 
        alignJob.add_var('ref', refGenomes['YPS1009'])
        alignJob.add_var('job', str(num))
        alignJob.add_parent(alignNCJob)
        mydag.add_job(alignJob)
        num+=1     
        
        # 6th step count alignment start sites
        sortedSam = re.sub('.sam', '-sorted.sam' , alignedSam)
        countJob = Job('counts.jtf', 'job' + str(num))       
        countJob.pre_skip(1)
        countJob.add_var('job', str(num))
        countJob.add_var('sam', alignedSam)
        countJob.add_var('out', sortedSam)
        countJob.add_parent(alignJob)
        mydag.add_job(countJob)
        num += 1
            
    mydag.save('MasterDagman.dsf')      # write the dag submit file
    
    # write the required condor files
    runFastqc()
    runCutAdapt()
    runRemoveFirst()
    alignNonCoding()
    alignBowtie() 
    countStarts()         

if __name__ == "__main__":
    main()

