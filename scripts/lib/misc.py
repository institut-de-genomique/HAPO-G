from Bio import SeqIO

import os
import subprocess


def check_dependencies():
    missing_dependency = False
    print("Checking dependencies...", flush=True)
    FNULL = open(os.devnull, 'w')

    #Looking for BWA
    try:
        subprocess.call(["bwa", "mem"], stdout=FNULL, stderr=FNULL)
        print("\tFound BWA.", flush=True)
    except OSError:
        print("\tWARNING : BWA not found.", flush=True)
        missing_dependency = True

    #Looking for Samtools
    try:
        subprocess.call(["samtools"], stdout=FNULL, stderr=FNULL)
        print("\tFound Samtools.", flush=True)
    except OSError:
        print("\tWARNING : Samtools not found.", flush=True)
        missing_dependency = True

    #Looking for Bedtools
    try:
        subprocess.call(["bedtools", "-h"], stdout=FNULL, stderr=FNULL)
        print("\tFound Bedtools.", flush=True)
    except OSError:
        print("\tWARNING : Bedtools not found.", flush=True)
        missing_dependency = True

    if missing_dependency:
    	exit(-1)


def create_genome_file(input_genome):
    with open("genome.txt", "w") as out:
        for record in SeqIO.parse(open(input_genome), "fasta"):
            out.write("%s\t%s\n" % (record.id, len(record.seq)))