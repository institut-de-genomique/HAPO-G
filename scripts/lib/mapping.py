import os
import subprocess

def launch_mapping(input_bam, pe1, pe2, threads):
    ########## BWA INDEX ##########
    print("\nGenerating bwa index...", flush=True)
    cmd = ["bwa", "index", "genome_masked.fasta"]
    print(" ".join(cmd), flush=True)
    _ = subprocess.run(cmd, 
        stdout = open("bwa_index.o", "w"),
        stderr = open("bwa_index.e", "w"))

    ########## BWA MEM ##########
    print("\nLaunching mapping on masked genome...", flush=True)
    cmd = "bash -c 'bwa mem -t %s %s " % (threads, "genome_masked.fasta")

    cmd += "<(zcat"
    for pe in pe1:
        cmd += " %s" % pe
    cmd += ") "
    cmd += "<(zcat"
    for pe in pe2:
        cmd += " %s" % pe
    cmd += ")"
    cmd += " 1> aln.sam 2> bwa_mem.e'"

    print(cmd, flush=True)
    os.system(cmd)

    ########## SAMTOOLS SORT ##########
    print("\nSorting and indexing sam file...", flush=True)
    cmd = ["samtools", "sort", "-m", "10G", "-@", threads, "-o", "aln.sorted.bam", "aln.sam"]
    print(" ".join(cmd), flush=True)
    _ = subprocess.run(cmd, 
	    stdout = open("samtools_sort.o", "w"), 
	    stderr = open("samtools_sort.e", "w"))

    ########## SAMTOOLS INDEX ##########
    cmd = ["samtools", "index", "aln.sorted.bam"]
    print(" ".join(cmd), flush=True)
    _ = subprocess.run(cmd, 
	    stdout = open("samtools_index.o", "w"), 
	    stderr = open("samtools_index.e", "w"))


    ########## SAMTOOLS MERGE ##########
    print("\nMerging bams...", flush=True)
    cmd = ["samtools", "merge", "-@", threads, "merged.bam", input_bam, "aln.sorted.bam"]
    print(" ".join(cmd), flush=True)
    _ = subprocess.run(cmd, 
	    stdout = open("samtools_index.o", "w"), 
	    stderr = open("samtools_index.e", "w"))