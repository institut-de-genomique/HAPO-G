#!/usr/bin/env python3
from lib import mapping
from lib import pipeline

import argparse
import os
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="hapog", 
        description="\n\nHAPoG uses alignments produced by BWA (or any other aligner that produces SAM files) to polish the consensus of a genome assembly.",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True)

    mandatory_args = parser.add_argument_group("Mandatory arguments")
    mandatory_args.add_argument("--genome", "-g",  
        action="store", 
        dest="input_genome", 
        help="Input genome file to map reads to",
        default=None,
        required=True)
    mandatory_args.add_argument("--pe1",  
        action="append", 
        dest="pe1", 
        help="Fastq.gz paired-end file (pair 1, can be given multiple times)",
        default=None,
        required=True)
    mandatory_args.add_argument("--pe2", 
        action="append", 
        dest="pe2", 
        help="Fastq.gz paired-end file (pair 2, can be given multiple times)",
        default=None,
        required=True)

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("--output", "-o",  
        action="store", 
        dest="output_dir", 
        help="Output directory name",
        default="HAPoG_results",
        required=False)
    optional_args.add_argument("--threads", "-t",  
        action="store", 
        dest="threads", 
        help="Number of threads (used in BWA, Samtools and HAPoG)",
        default="8",
        required=False)

    args = parser.parse_args()
    pipeline.check_dependencies()

    args.input_genome = os.path.abspath(args.input_genome)
    args.output_dir = os.path.abspath(args.output_dir)

    pe1 = []
    for pe in args.pe1 :
        pe1.append(os.path.abspath(pe))
    pe2 = []
    for pe in args.pe2 :
        pe2.append(os.path.abspath(pe))  

    try:
        os.mkdir(args.output_dir)
    except:
        pass
    os.chdir(args.output_dir)

    try:
      os.mkdir("bam")
    except:
      pass

    try:
      os.mkdir("logs")
    except:
      pass

    mapping.launch_mapping(args.input_genome, pe1, pe2, args.threads)
    pipeline.create_chunks(args.input_genome, args.threads)
    pipeline.extract_bam(int(args.threads))
    pipeline.launch_hapog()
    pipeline.merge_results()
