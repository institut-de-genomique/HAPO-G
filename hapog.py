#!/usr/bin/env python3
from lib import mapping
from lib import pipeline

import argparse
import os
import sys
import time


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
        required=False)
    mandatory_args.add_argument("--pe2", 
        action="append", 
        dest="pe2", 
        help="Fastq.gz paired-end file (pair 2, can be given multiple times)",
        default=None,
        required=False)

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-b",  
        action="store", 
        dest="bam_file", 
        help="Skip mapping step and provide a sorted bam file",
        default="",
        required=False)
    optional_args.add_argument("-u",  
        action="store_true", 
        dest="include_unpolished", 
        help="Include unpolished sequences in final output",
        default=False,
        required=False)
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
    optional_args.add_argument("--bin",  
        action="store", 
        dest="hapog_bin", 
        help="Use a different HAPoG binary (for debug purposes)",
        default=None,
        required=False)

    args = parser.parse_args()
    pipeline.check_dependencies()

    args.input_genome = os.path.abspath(args.input_genome)
    args.output_dir = os.path.abspath(args.output_dir)

    pe1 = []
    pe2 = []
    if args.bam_file:
        args.bam_file = os.path.abspath(args.bam_file)
    else:
        if not args.pe1 or not args.pe2:
            print("You need to specify the paths to paired-end read files.")
            sys.exit(-1)

        for pe in args.pe1 :
            pe1.append(os.path.abspath(pe))
        for pe in args.pe2 :
            pe2.append(os.path.abspath(pe))  

    try:
        os.mkdir(args.output_dir)
    except:
        print(f"\nOutput directory {args.output_dir} can't be created, please erase it before launching HAPoG.\n")
        sys.exit(1)
    os.chdir(args.output_dir)

    os.mkdir("bam")
    os.mkdir("logs")
    os.mkdir("cmds")

    global_start = time.perf_counter()

    non_alphanumeric_chars = False
    if not args.bam_file:
        non_alphanumeric_chars = pipeline.check_fasta_headers(args.input_genome)
        if non_alphanumeric_chars:
            print("\nNon alphanumeric characters detected in fasta headers. Renaming sequences.", flush=True)
            pipeline.rename_assembly(args.input_genome)
        else:
            os.system(f"ln -s {args.input_genome} assembly.fasta")
        mapping.launch_mapping("assembly.fasta", pe1, pe2, args.threads)
    else:
        if pipeline.check_fasta_headers(args.input_genome):
            print("\nERROR: Non-alphanumeric characters detected in fasta headers will cause samtools view to crash.", flush=True)
            print("Please remove these characters before launching Hapo-G with -b or let Hapo-G do the mapping by itself.", flush=True)
            print("Authorized characters belong to this list: 'a-z', 'A-Z', '0-9', '_-'.", flush=True)
            sys.exit(-1)
        os.system(f"ln -s {args.input_genome} assembly.fasta")
        os.system(f"ln -s {args.bam_file} bam/aln.sorted.bam")
        mapping.index_bam()

    if int(args.threads) > 1:
        pipeline.create_chunks("assembly.fasta", args.threads)
        pipeline.extract_bam(int(args.threads))
    else:
        os.system(f"ln -s {args.input_genome} chunks/chunks_1.fasta")
        os.system(f"ln -s bam/aln.sorted.bam chunks_bam/chunks_1.bam")

    pipeline.launch_hapog(args.hapog_bin)
    pipeline.merge_results(int(args.threads))
        
    if non_alphanumeric_chars:
        pipeline.rename_results()
    else:
        os.system("mv HAPoG_results/hapog.changes.tmp HAPoG_results/hapog.changes")
        os.system("mv HAPoG_results/hapog.fasta.tmp HAPoG_results/hapog.fasta")

    if args.include_unpolished:
        pipeline.include_unpolished(args.input_genome)

    print("\nResults can be found in the HAPoG_results directory")
    print(f"Total running time: {int(time.perf_counter() - global_start)} seconds")
    print("\nThanks for using HAPoG, have a great day :-)\n")     
