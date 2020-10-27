from Bio import SeqIO

import glob
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

    if missing_dependency:
    	exit(-1)


def get_genome_size(genome):
    cumul_size = 0
    for record in SeqIO.parse(open(genome), "fasta"):
        cumul_size += len(record.seq)
    return(cumul_size)


def create_chunks(genome, threads):
    cumul_size = get_genome_size(genome)

    chunk_size = cumul_size / int(threads)
    print(f"\nFragmenting the genome into {threads} chunks of {int(chunk_size):,} bases")
    try:
        os.mkdir("chunks")
    except:
        pass

    current_chunk = 1
    current_chunk_size = 0
    current_chunk_file = open("chunks/chunk_1.fasta", "w")
    current_bed_file = open("chunks/chunk_1.bed", "w")

    for record in SeqIO.parse(open(genome), "fasta"):
        if current_chunk_size >= chunk_size and current_chunk != threads:
            current_chunk_file.close()
            current_chunk_file = open(f"chunks/chunks_{current_chunk + 1}.fasta", "w")
            current_bed_file = open(f"chunks/chunks_{current_chunk + 1}.bed", "w")
            current_chunk += 1
            current_chunk_size = 0

        current_chunk_file.write(record.format("fasta"))
        current_bed_file.write(f"{record.id}\t0\t{len(record.seq)}\n")
        current_chunk_size += len(record.seq)
 
    current_chunk_file.close()
    current_bed_file.close()


def extract_bam():
    print("Extracting bam for each chunk")
    try:
        os.mkdir("chunks_bam")
    except:
        pass

    cmds = []
    nb_beds = 0
    for bed in glob.glob("chunks/*.bed"):
        bam = "chunks_bam/" + bed.split("/")[1].replace(".bed", ".bam")

        nb_beds += 1
        cmds.append(["samtools", "view", "-b", "bam/aln.sorted.bam", "-L", bed])

    for i in range(0, nb_beds, 5):
        procs = []
        
        for j in range(0, min(5, nb_beds-i-5)):
            cmd = cmds[i+j]
            bam = "chunks_bam/" + cmd[5].split("/")[1].replace(".bed", ".bam")
            print(" ".join(cmd), flush=True)
            procs.append(subprocess.Popen(cmd,
                stdout = open(bam, "w"),
                stderr = open("logs/samtools_sort.e", "a")))
        
        for p in procs:
            p.wait()
