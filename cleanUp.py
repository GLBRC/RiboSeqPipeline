#!/usr/bin/env python
"""cleanUp.py

Clean up the directory after running the Ribo-Seq pipeline.

Parameters
----------
None

Example
-------
    usage:

        cleanUp.py
"""
import argparse 
import glob
import os
import re
import shutil
import sys


def main():
    cmdparser = argparse.ArgumentParser(description="Ribo-Seq pipeline clean up step.",
                                        usage='%(prog)s ' ,prog='cleanUp.py'  )
    
    cwd = os.getcwd() + '/'
    
    if not os.path.exists(cwd + 'condor'):
        os.mkdir('condor')
    csDir = cwd + 'condor'
    for sb in os.listdir():
        if sb.endswith('.jtf'):
            os.rename((cwd + '/' + sb), (csDir + '/' + sb))
    
    if not os.path.exists(cwd + 'results'):
        os.mkdir('results')        
    resDir = cwd + 'results' + '/'
    for res in os.listdir('alignments'):
        if res.endswith('.txt'):
            os.rename(cwd + 'alignments/' + res, resDir + res)
            
        
    
        

if __name__ == "__main__":
    main()
