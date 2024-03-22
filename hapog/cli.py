#!/usr/bin/env python3
from hapog import mapping
from hapog import pipeline

import argparse
import os
import sys
import time


def main():
    parser = argparse.ArgumentParser(
        prog="hapog",
        description="\n\nHapo-G uses alignments produced by BWA (or any other aligner that produces SAM files) to polish the consensus of a genome assembly.",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )

    mandatory_args = parser.add_argument_group("Mandatory arguments")
    mandatory_args.add_argument(
        "--genome",
        "-g",
        action="store",
        dest="input_genome",
        help="Input genome file to map reads to",
        default=None,
        required=True,
    )
    mandatory_args.add_argument(
        "--pe1",
        action="append",
        dest="pe1",
        help="Fastq.gz paired-end file (pair 1, can be given multiple times)",
        default=None,
        required=False,
    )
    mandatory_args.add_argument(
        "--pe2",
        action="append",
        dest="pe2",
        help="Fastq.gz paired-end file (pair 2, can be given multiple times)",
        default=None,
        required=False,
    )
    mandatory_args.add_argument(
        "--single",
        action="store",
        dest="long_reads",
        help="Use long reads instead of short reads (can only be given one time, please concatenate all read files into one)",
        default=None,
        required=False,
    )

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument(
        "-b",
        action="store",
        dest="bam_file",
        help="Skip mapping step and provide a sorted bam file. Important: the BAM file must not contain secondary alignments, please use the 'secondary=no' option in Minimap2.",
        default="",
        required=False,
    )
    optional_args.add_argument(
        "-u",
        action="store_true",
        dest="include_unpolished",
        help="Include unpolished sequences in final output",
        default=False,
        required=False,
    )
    optional_args.add_argument(
        "--output",
        "-o",
        action="store",
        dest="output_dir",
        help="Output directory name",
        default="hapog_results",
        required=False,
    )
    optional_args.add_argument(
        "--threads",
        "-t",
        action="store",
        dest="threads",
        help="Number of threads (used in BWA, Samtools and Hapo-G)",
        default="8",
        required=False,
    )
    optional_args.add_argument(
        "--hapog-threads",
        action="store",
        dest="hapog_threads",
        help="Maximum number of Hapo-G jobs to launch in parallel (Defaults to the same value as --threads)",
        default=0,
        type=int,
        required=False,
    )
    optional_args.add_argument(
        "--bin",
        action="store",
        dest="hapog_bin",
        help="Use a different Hapo-G binary (for debug purposes)",
        default=None,
        required=False,
    )
    optional_args.add_argument(
        "--samtools-mem",
        action="store",
        dest="samtools_mem",
        help="Amount of memory to use per samtools thread (Default: '5G')",
        default="5G",
        required=False,
    )

    args = parser.parse_args()
    pipeline.check_dependencies()

    args.input_genome = os.path.abspath(args.input_genome)
    args.output_dir = os.path.abspath(args.output_dir)
    if args.hapog_threads == 0:
        args.hapog_threads = args.threads

    pe1 = []
    pe2 = []
    single = []
    use_short_reads = False

    if args.bam_file:
        args.bam_file = os.path.abspath(args.bam_file)
    else:
        if not args.long_reads and (not args.pe1 or not args.pe2):
            print("You need to specify the paths to paired-end or long reads files.")
            sys.exit(-1)

        if not args.long_reads:
            for pe in args.pe1:
                pe1.append(os.path.abspath(pe))
            for pe in args.pe2:
                pe2.append(os.path.abspath(pe))
            use_short_reads = True
        else:
            args.long_reads = os.path.abspath(args.long_reads)
            if not os.path.exists(args.long_reads):
                print("Long reads not found: %s" % (args.long_reads))
                sys.exit(-1)

    try:
        os.mkdir(args.output_dir)
    except:
        print(
            f"\nOutput directory {args.output_dir} can't be created, please erase it before launching Hapo-G.\n"
        )
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
            print(
                "\nNon alphanumeric characters detected in fasta headers. Renaming sequences.",
                flush=True,
            )
            pipeline.rename_assembly(args.input_genome)
        else:
            os.system(f"ln -s {args.input_genome} assembly.fasta")

        if use_short_reads:
            mapping.launch_PE_mapping("assembly.fasta", pe1, pe2, args.threads, args.samtools_mem)
        else:
            mapping.launch_LR_mapping("assembly.fasta", args.long_reads, args.threads, args.samtools_mem)

    else:
        if pipeline.check_fasta_headers(args.input_genome):
            print(
                "\nERROR: Non-alphanumeric characters detected in fasta headers will cause samtools view to crash.",
                flush=True,
            )
            print(
                "Please remove these characters before launching Hapo-G with -b or let Hapo-G do the mapping by itself.",
                flush=True,
            )
            print(
                "Authorized characters belong to this list: 'a-z', 'A-Z', '0-9', '_-'.",
                flush=True,
            )
            sys.exit(-1)

        os.system(f"ln -s {args.input_genome} assembly.fasta")
        os.system(f"ln -s {args.bam_file} bam/aln.sorted.bam")
        mapping.index_bam()

    if int(args.hapog_threads) > 1:
        pipeline.create_chunks("assembly.fasta", args.threads)
        pipeline.extract_bam(int(args.threads))
    else:
        os.mkdir("chunks")
        os.mkdir("chunks_bam")
        if non_alphanumeric_chars:
            os.system("ln -s ../assembly.fasta chunks/chunks_1.fasta")
        else:
            os.system(f"ln -s {args.input_genome} chunks/chunks_1.fasta")
        os.system(f"ln -s ../bam/aln.sorted.bam chunks_bam/chunks_1.bam")

    pipeline.launch_hapog(args.hapog_bin, args.hapog_threads)
    pipeline.merge_results(int(args.threads))

    if non_alphanumeric_chars:
        pipeline.rename_results()
    else:
        os.system("mv hapog_results/hapog.changes.tmp hapog_results/hapog.changes")
        os.system("mv hapog_results/hapog.fasta.tmp hapog_results/hapog.fasta")

    if args.include_unpolished:
        pipeline.include_unpolished(args.input_genome)

    print("\nResults can be found in the hapog_results directory")
    print(f"Total running time: {int(time.perf_counter() - global_start)} seconds")
    print("\nThanks for using Hapo-G, have a great day :-)\n")
