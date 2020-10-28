from Bio import SeqIO

import glob
import os
import shutil
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
    print(f"\nFragmenting the genome into {threads} chunks of {int(chunk_size):,} bases (if scaffolds sizes permit it)")
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


def extract_bam(processes):
    print("\nExtracting bam for each chunk")
    try:
        os.mkdir("chunks_bam")
    except:
        pass

    cmds = []
    nb_beds = 0
    for bed in glob.glob("chunks/*.bed"):
        bam = "chunks_bam/" + bed.split("/")[1].replace(".bed", ".bam")

        nb_beds += 1
        cmds.append(["samtools", "view", "-ML", bed, "-b", "bam/aln.sorted.bam"])

    while len(cmds) > 0:
        procs = []
        
        for j in range(0, min(min(processes, 8), len(cmds))):
            cmd = cmds.pop(0)
            bam = "chunks_bam/" + cmd[3].split("/")[1].replace(".bed", ".bam")
            print(" ".join(cmd), flush=True)
            procs.append(subprocess.Popen(cmd,
                stdout = open(bam, "w"),
                stderr = open("logs/samtools_split.e", "a")))
        
        for p in procs:
            p.wait()


def launch_hapog():
    print(f"\nLaunching HAPoG on each chunk")
    try:
        os.mkdir("HAPoG_chunks")
    except:
        pass

    script_path = os.path.realpath(__file__).replace("/lib/misc.py", "")
    procs = []
    for chunk in glob.glob("chunks/*.fasta"):
        chunk_prefix = chunk.split("/")[-1].replace(".fasta", "")
        cmd = [
            f"{script_path}/build/hapog", 
            "-b", f"chunks_bam/{chunk_prefix}.bam", 
            "-f", chunk, 
            "-o", f"HAPoG_chunks/{chunk_prefix}.fasta",
            "-c", f"HAPoG_chunks/{chunk_prefix}.changes"
        ]
        print(" ".join(cmd), flush=True)
        procs.append(subprocess.Popen(cmd,
            stdout = open(f"logs/hapog_{chunk_prefix}.o", "w"),
            stderr = open(f"logs/hapog_{chunk_prefix}.e", "w")))
      
    for p in procs:
        p.wait()


def merge_results():
    print("\nMerging results")
    try:
        os.mkdir("HAPoG_results")
    except:
        pass

    with open("HAPoG_results/hapog.fasta", "wb") as out:
        for f in glob.glob("HAPoG_chunks/*.fasta"):
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, out)
                out.write(b"\n")

    with open("HAPoG_results/hapog.changes", "wb") as out:
        for f in glob.glob("HAPoG_chunks/*.changes"):
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, out)
                out.write(b"\n")

    print("Done.")
    print("Results can be found in the HAPoG_results directory\n")
    print("Thanks for using HAPoG, have a great day :-)")
