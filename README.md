# RiboSeqPipeline

Generate tables of read start counts for Ribo-Seq data.

## Getting started

Pipeline is designed to be run using GLBRC's HTCondor compute cluster.
After logging into scarcity-submit activate the riboSeq conda environment.

conda activate /home/glbrc.org/mplace/.conda/envs/riboSeq

Please run this pipeline in a dedicated directory and create a subdirectory
called fastq. Put all of your Read 1 (fwd) fastq files in this directory.
If fastq files are .gz, decompress them. Create a text file containing the
fastq names.  This file will be used to run the pipeline.

ls fastq > fastqFiles.txt

MyProject/ -> fastq/*R1.fastq

## Running the Pipeline

RiboSeqPipeline.py -f fastqFiles.txt

## Results
The following directories will be created:
 
 alignments  -- alignment to reference genome results (.sam files)

 alignNonCodingRNA -- filtering step, alignment of intial fastq to Non-Coding RNA.

 condor    -- condor files

 cutadapt  -- fastq files filtered with cutadapt and have additionally had the 1st base
              removed if quality score <= 10.

 log       -- log files for each step

 results   -- contain the count tables for each sample.

 ## Requirements

 See RiboSeq_condaEnv.yml

 Pipeline designed to run as HTCondor Dagman job.

 A bowtie2 built reference genome.
 
*******************************************************************************


## Preliminary pipeline for development and testing using Jupyter Notebook (do not use.)

Designed used with MicroSoft VS Code and to be used with on GLBRC scarcity with HTCondor.  

1) Download VS Code : https://code.visualstudio.com

2) Make sure you have an OpenSSH compatible SSH Client installed locally. On MacOS this should be already be installed.

3) In VS Code install the Remote-SSH extension. If you plan to work with other remote extensions in VS Code, you may choose to install the Remote Development extension pack.

4) In VS Code, select Remote-SSH: Connect to Host... from the Command Palette (F1, Ctrl+Shift+P) and use your user@hostname.

    If you have trouble with the above instructions take a look here: https://code.visualstudio.com/docs/remote/ssh
    for more details on how to get remote ssh development with VS Code started.

6) Open the folder you want to use.  I recommend using git to clone this repository. 
    git clone git@gitpub.wei.wisc.edu:mplace/riboSeqPipeline.git

7) Set up your input fastq files within a subdirectory called "data".

8) Try to use my riboSeq conda environment in VS Code.

 If this doesn't work used conda to create your own conda riboSeq environment.
 Use this file: RiboSeq_condaEnv.yml  located in the git repo.

 Edit the last line in the yml file to reflect your directory.
 prefix: /home/glbrc.org/USER/.conda/envs/riboSeq

 Then to create the environment:  conda create -f RiboSeq_condaEnv.yml
 This can then be used within your VS Code. 

 9) Open RiboSeq_pipeline.ipynb and walk through the pipeline.
 Note that most steps will create condor submit files along with a shell script.
 Just use scarcity-submit.glbrc.org and run condor_submit mySubmit.txt and wait for 
 the job to finish before running the next step.

 

