from Bio import SeqIO

import glob
import os
import shutil
import subprocess
import time


def check_dependencies():
    missing_dependency = False
    print("\nChecking dependencies...", flush=True)
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


def check_fasta_headers(genome):
    authorized_chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-_"
    for line in open(genome):
        if line.startswith(">"):
            header = line[1:].rstrip("\n")
            for char in header:
                if char not in authorized_chars:
                    return True
    return False


def get_genome_size(genome):
    cumul_size = 0
    for record in SeqIO.parse(open(genome), "fasta"):
        cumul_size += len(record.seq)
    return(cumul_size)


def rename_assembly(genome):
    correspondance_file = open("correspondance.txt", "w")
    with open("assembly.fasta", "w") as out:
        counter = 0
        for line in open(genome):
            if line.startswith(">"):
                out.write(f">Contig{counter}\n")
                correspondance_file.write(f"Contig{counter}\t{line[1:]}")
                counter += 1
            else:
                out.write(line)
    correspondance_file.close()


def create_chunks(genome, threads):
    cumul_size = get_genome_size(genome)

    chunk_size = cumul_size / int(threads)
    print(f"\nFragmenting the genome into {threads} chunks of {int(chunk_size):,} bases (if scaffolds sizes permit it)", flush=True)
    try:
        os.mkdir("chunks")
    except:
        pass

    current_chunk = 1
    current_chunk_size = 0
    current_chunk_file = open("chunks/chunks_1.fasta", "w")
    current_bed_file = open("chunks/chunks_1.bed", "w")

    start = time.perf_counter()
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
    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)

    current_chunk_file.close()
    current_bed_file.close()


def extract_bam(processes):
    print("\nExtracting bam for each chunk", flush=True)
    try:
        os.mkdir("chunks_bam")
    except:
        pass

    start = time.perf_counter()

    cmds = []
    nb_beds = 0
    for bed in glob.glob("chunks/*.bed"):
        nb_beds += 1
        cmds.append(["samtools", "view", "-ML", bed, "-b", "bam/aln.sorted.bam"])

    procs = []
    for j in range(0, len(cmds)):
        cmd = cmds.pop(0)
        bam = "chunks_bam/" + cmd[3].split("/")[1].replace(".bed", ".bam")
        print(" ".join(cmd), flush=True, file=open("cmds/extract_bam.cmds", "a"))
        procs.append(subprocess.Popen(cmd,
            stdout = open(bam, "w"),
            stderr = open("logs/samtools_split.e", "a")))
    
    for p in procs:
        p.wait()
        if p.returncode != 0:
            print(f"ERROR: Samtools view didn't finish correctly, return code: {p.returncode}")
            print("Faulty command: {p.args}")
            exit(1)

    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)


def launch_hapog(hapog_bin):
    print(f"\nLaunching HAPoG on each chunk", flush=True)
    try:
        os.mkdir("HAPoG_chunks")
    except:
        pass

    start = time.perf_counter()

    script_path = os.path.realpath(__file__).replace("/lib/pipeline.py", "")
    if not hapog_bin:
        hapog_bin = f"{script_path}/bin/hapog"
    else:
        print(f"Using this bin: {hapog_bin}" )

    procs = []
    for chunk in glob.glob("chunks/*.fasta"):
        chunk_prefix = chunk.split("/")[-1].replace(".fasta", "")
        cmd = [
            hapog_bin, 
            "-b", f"chunks_bam/{chunk_prefix}.bam", 
            "-f", chunk, 
            "-o", f"HAPoG_chunks/{chunk_prefix}.fasta",
            "-c", f"HAPoG_chunks/{chunk_prefix}.changes"
        ]
        print(" ".join(cmd), flush=True, file=open("cmds/hapog.cmds", "a"))
        procs.append(subprocess.Popen(cmd,
            stdout = open(f"logs/hapog_{chunk_prefix}.o", "w"),
            stderr = open(f"logs/hapog_{chunk_prefix}.e", "w")))
      
    for p in procs:
        p.wait()
        return_code = p.returncode
        if return_code != 0:
            print(f"ERROR: HAPoG didn't finish successfully, exit code: {return_code}")
            print("Faulty command: %s" % (" ".join(p.args)))
            exit(1)

    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)


def merge_results(threads):
    print("\nMerging results", flush=True)
    try:
        os.mkdir("HAPoG_results")
    except:
        pass

    start = time.perf_counter()

    with open("HAPoG_results/hapog.fasta.tmp", "wb") as out:
        for i in range(1, threads+1):
            f = f"HAPoG_chunks/chunks_{i}.fasta"
            if os.path.exists(f):
                with open(f,'rb') as fd:
                    shutil.copyfileobj(fd, out)
                    out.write(b"\n")

    with open("HAPoG_results/hapog.changes.tmp", "wb") as out:
        for i in range(1, threads+1):
            f = f"HAPoG_chunks/chunks_{i}.changes"
            if os.path.exists(f):
                with open(f,'rb') as fd:
                    shutil.copyfileobj(fd, out)
                    out.write(b"\n")

    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)


def rename_results():
    correspondance_file = open("correspondance.txt")
    dict_correspondance = {}
    for line in correspondance_file:
        new, original = line.strip("\n").split("\t")
        dict_correspondance[new] = original
    correspondance_file.close() 

    with open("HAPoG_results/hapog.fasta", "w") as out:
        for record in SeqIO.parse(open("HAPoG_results/hapog.fasta.tmp"), "fasta"):
            out.write(f">{dict_correspondance[str(record.id).replace('_polished', '')]}\n{record.seq}\n")

    with open("HAPoG_results/hapog.changes", "w") as out:
        for line in open("HAPoG_results/hapog.changes.tmp"):
            line = line.strip("\n").split("\t")
            try:
                line[0] = dict_correspondance[line[0]]
            except:
                continue
            line = "\t".join(line)
            out.write(f"{line}\n")
    
    for f in glob.glob("assembly.fasta*"):
        os.remove(f)
    os.remove("HAPoG_results/hapog.fasta.tmp")
    os.remove("HAPoG_results/hapog.changes.tmp")


def include_unpolished(genome):
    print("\nWriting unpolished contigs to final output...")
    start = time.perf_counter()

    initial_contig_names = set()
    if os.path.exists("correspondance.txt"):
        for line in open("correspondance.txt"):
            _, initial = line.strip("\n").split("\t")
            initial_contig_names.add(initial)
    else:
        for line in open(genome):
            if line.startswith(">"):
                initial_contig_names.add(line[1:].strip("\n"))

    polished_contig_names = set()
    for line in open("HAPoG_results/hapog.fasta"):
        if line.startswith(">"):
            contig_name = line[1:].strip("\n").replace("_polished", "")
            polished_contig_names.add(contig_name)

    with open("HAPoG_results/hapog.fasta", "a") as out:
        for record in SeqIO.parse(open(genome), "fasta"):
            if record.description.replace("_polished", "") not in polished_contig_names:
                out.write(record.format("fasta"))

    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)
