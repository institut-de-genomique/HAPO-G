import os
import subprocess

def launch_mapping(genome, pe1, pe2, threads):
    ########## BWA INDEX ##########
    print("\nGenerating bwa index...", flush=True)
    cmd = ["bwa", "index", genome]
    print(" ".join(cmd), flush=True)
    _ = subprocess.run(cmd, 
        stdout = open("logs/bwa_index.o", "w"),
        stderr = open("logs/bwa_index.e", "w"))

    ########## BWA MEM ##########
    print("\nLaunching mapping on genome...", flush=True)
    cmd = "bash -c 'bwa mem -t %s %s " % (threads, genome)

    cmd += "<(zcat"
    for pe in pe1:
        cmd += " %s" % pe
    cmd += ") "
    cmd += "<(zcat"
    for pe in pe2:
        cmd += " %s" % pe
    cmd += ")"
    cmd += " 1> bam/aln.sam 2> logs/bwa_mem.e'"

    print(cmd, flush=True)
    os.system(cmd)

    ########## SAMTOOLS SORT ##########
    print("\nSorting and indexing sam file...", flush=True)
    cmd = ["samtools", "sort", "-m", "10G", "-@", threads, "-o", "bam/aln.sorted.bam", "bam/aln.sam"]
    print(" ".join(cmd), flush=True)
    _ = subprocess.run(cmd, 
	    stdout = open("logs/samtools_sort.o", "w"), 
	    stderr = open("logs/samtools_sort.e", "w"))

    ########## SAMTOOLS INDEX ##########
    cmd = ["samtools", "index", "bam/aln.sorted.bam"]
    print(" ".join(cmd), flush=True)
    _ = subprocess.run(cmd, 
	    stdout = open("logs/samtools_index.o", "w"), 
	    stderr = open("logs/samtools_index.e", "w"))
