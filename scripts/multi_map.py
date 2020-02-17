from lib import mapping
from lib import masking
from lib import misc

import argparse
import os
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="mask_and_map.py", 
        description="\n\nUses bedtools to extract < 5-cov regions, masks the other regions in the genome, maps reads on it and merges bams",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)

    mandatory_args = parser.add_argument_group("Mandatory arguments")
    mandatory_args.add_argument("--bam", "-b",  
        action="store", 
        dest="input_bam", 
        help="Input BAM file that contains mappings",
        default=None,
        required=True)
    mandatory_args.add_argument("--genome", "-g",  
        action="store", 
        dest="input_genome", 
        help="Input genome file to map reads to",
        default=None,
        required=True)
    mandatory_args.add_argument("--pe1",  
        action="append", 
        dest="pe1", 
        help="Fastq.gz paired-end file (pair 1)",
        default=None,
        required=True)
    mandatory_args.add_argument("--pe2", 
        action="append", 
        dest="pe2", 
        help="Fastq.gz paired-end file (pair 2)",
        default=None,
        required=True)
    mandatory_args.add_argument("--output", "-o",  
        action="store", 
        dest="output_dir", 
        help="Output directory name",
        default="multi_map",
        required=False)
    mandatory_args.add_argument("--threads", "-t",  
        action="store", 
        dest="threads", 
        help="Number of threads to run BWA and samtools on",
        default="8",
        required=False)

    args = parser.parse_args()
    misc.check_dependencies()

    args.input_genome = os.path.abspath(args.input_genome)
    args.input_bam = os.path.abspath(args.input_bam)

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

    misc.create_genome_file(args.input_genome)
    masking.mask_genome(args.input_bam, args.input_genome)	
    mapping.launch_mapping(args.input_bam, pe1, pe2, args.threads)