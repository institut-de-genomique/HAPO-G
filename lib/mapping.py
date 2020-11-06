import os
import subprocess
import time


def launch_mapping(genome, pe1, pe2, threads):
    ########## BWA INDEX ##########
    print("\nGenerating bwa index...", flush=True)
    cmd = ["bwa", "index", genome]
    
    start = time.perf_counter()
    print(" ".join(cmd), flush=True, file=open("cmds/bwa_index.cmds", "a"))
    _ = subprocess.run(cmd, 
        stdout = open("logs/bwa_index.o", "w"),
        stderr = open("logs/bwa_index.e", "w"))
    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)

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

    start = time.perf_counter()
    print(cmd, flush=True, file=open("cmds/bwa_mem.cmds", "a"))
    os.system(cmd)
    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)

    ########## SAMTOOLS SORT ##########
    print("\nSorting and indexing sam file...", flush=True)
    cmd = ["samtools", "sort", "-m", "10G", "-@", threads, "-o", "bam/aln.sorted.bam", "bam/aln.sam"]

    start = time.perf_counter()
    print(" ".join(cmd), flush=True, file=open("cmds/samtools_sort.cmds", "a"))
    _ = subprocess.run(cmd, 
        stdout = open("logs/samtools_sort.o", "w"), 
        stderr = open("logs/samtools_sort.e", "w"))
    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)

    ########## SAMTOOLS INDEX ##########
    cmd = ["samtools", "index", "bam/aln.sorted.bam"]

    start = time.perf_counter()
    print(" ".join(cmd), flush=True, file=open("cmds/samtools_index.cmds", "a"))
    _ = subprocess.run(cmd, 
        stdout = open("logs/samtools_index.o", "w"), 
        stderr = open("logs/samtools_index.e", "w"))
    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)
