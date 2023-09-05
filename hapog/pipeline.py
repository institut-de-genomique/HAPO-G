from Bio import SeqIO

import glob
import os
import shutil
import subprocess
import time
import warnings


def is_in_path(tool):
    return shutil.which(tool) is not None


def check_dependencies():
    missing_dependency = False
    print("\nChecking dependencies...", flush=True)

    tools = ["bwa", "samtools"]
    for tool in tools:
        if not is_in_path(tool):
            print(f"\tWARNING: {tool} not found.", flush=True)
            missing_dependency = True
        else:
            print(f"\tFound {tool}", flush=True)

    if missing_dependency:
        exit(-1)


def check_fasta_headers(genome):
    authorized_chars = (
        "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-_"
    )
    with open(genome) as genome_file:
        for line in genome_file:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                for char in header:
                    if char not in authorized_chars:
                        return True
    return False


def get_genome_size(genome):
    cumul_size = 0
    with open(genome) as genome_file:
        for record in SeqIO.parse(genome_file, "fasta"):
            cumul_size += len(record.seq)
    return cumul_size


def rename_assembly(genome):
    correspondance_file = open("correspondance.txt", "w")
    with open("assembly.fasta", "w") as out:
        counter = 0
        with open(genome) as genome_file:
            for line in genome_file:
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
    print(
        f"\nFragmenting the genome into {threads} chunks of {int(chunk_size):,} bases (depending of scaffold sizes)",
        flush=True,
    )
    try:
        os.mkdir("chunks")
    except:
        pass

    current_chunk = 1
    current_chunk_size = 0
    current_chunk_file = open("chunks/chunks_1.fasta", "w")
    current_bed_file = open("chunks/chunks_1.bed", "w")

    start = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

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
        with open("cmds/extract_bam.cmds", "a") as cmd_file:
            print(" ".join(cmd), flush=True, file=cmd_file)

        # Ignore the ResourceWarning about unclosed files
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            procs.append(
                subprocess.Popen(
                    cmd, stdout=open(bam, "w"), stderr=open("logs/samtools_split.e", "a")
                )
            )

    has_failed = False
    for p in procs:
        p.wait()

        if p.returncode != 0:
            print(
                f"ERROR: Samtools view didn't finish correctly, return code: {p.returncode}"
            )
            print("Faulty command: {p.args}")
            has_failed = True
    
    if has_failed:
        exit(1)

    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)


def launch_hapog(hapog_bin, parallel_jobs):
    print(f"\nLaunching Hapo-G on each chunk", flush=True)
    try:
        os.mkdir("hapog_chunks")
    except:
        pass

    start = time.perf_counter()

    script_path = os.path.realpath(__file__).replace("/hapog/pipeline.py", "")
    if not hapog_bin:
        if not is_in_path("hapog_bin"):
            hapog_bin = f"{script_path}/hapog_build/hapog"
        else:
            hapog_bin = "hapog_bin"
    else:
        print(f"Using this bin: {hapog_bin}")

    procs = []
    for chunk in glob.glob("chunks/*.fasta"):
        chunk_prefix = chunk.split("/")[-1].replace(".fasta", "")
        cmd = [
            hapog_bin,
            "-b",
            f"chunks_bam/{chunk_prefix}.bam",
            "-f",
            chunk,
            "-o",
            f"hapog_chunks/{chunk_prefix}.fasta",
            "-c",
            f"hapog_chunks/{chunk_prefix}.changes",
        ]
        with open("cmds/hapog.cmds", "a") as cmd_file:
            print(" ".join(cmd), flush=True, file=cmd_file)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            procs.append(
                subprocess.Popen(
                    cmd,
                    stdout=open(f"logs/hapog_{chunk_prefix}.o", "w"),
                    stderr=open(f"logs/hapog_{chunk_prefix}.e", "w"),
                )
            )

        # Only launch a job if there is less than 'parallel_jobs' running
        # Otherwise, wait for any to finish before launching a new one
        while len([p for p in procs if p.poll() is None]) >= parallel_jobs:
            time.sleep(10)


    has_failed = False
    for p in procs:
        p.wait()

        return_code = p.returncode
        if return_code != 0:
            print(f"ERROR: Hapo-G didn't finish successfully, exit code: {return_code}")
            print("Faulty command: %s" % (" ".join(p.args)))
            has_failed = True

    if has_failed:
        exit(1)

    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)


def merge_results(threads):
    print("\nMerging results", flush=True)
    try:
        os.mkdir("hapog_results")
    except:
        pass

    start = time.perf_counter()

    with open("hapog_results/hapog.fasta.tmp", "w") as out:
        for i in range(1, threads + 1):
            f = f"hapog_chunks/chunks_{i}.fasta"
            if os.path.exists(f):
                with open(f, "r") as fd:
                    shutil.copyfileobj(fd, out)
                    out.write("\n")

    with open("hapog_results/hapog.changes.tmp", "w") as out:
        for i in range(1, threads + 1):
            f = f"hapog_chunks/chunks_{i}.changes"
            if os.path.exists(f):
                with open(f, "r") as fd:
                    shutil.copyfileobj(fd, out)
                    out.write("\n")

    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)


def rename_results():
    correspondance_file = open("correspondance.txt")
    dict_correspondance = {}
    for line in correspondance_file:
        new, original = line.strip("\n").split("\t")
        dict_correspondance[new] = original
    correspondance_file.close()

    with open("hapog_results/hapog.fasta", "w") as out:
        hapog_tmp = open("hapog_results/hapog.fasta.tmp")
        for record in SeqIO.parse(hapog_tmp, "fasta"):
            out.write(
                f">{dict_correspondance[str(record.id).replace('_polished', '')]}\n{record.seq}\n"
            )
        hapog_tmp.close()

    with open("hapog_results/hapog.changes", "w") as out:
        hapog_tmp = open("hapog_results/hapog.changes.tmp")
        for line in hapog_tmp:
            line = line.strip("\n").split("\t")
            try:
                line[0] = dict_correspondance[line[0]]
            except:
                continue
            line = "\t".join(line)
            out.write(f"{line}\n")
        hapog_tmp.close()

    for f in glob.glob("assembly.fasta*"):
        os.remove(f)
    os.remove("hapog_results/hapog.fasta.tmp")
    os.remove("hapog_results/hapog.changes.tmp")


def include_unpolished(genome):
    print("\nWriting unpolished contigs to final output...")
    start = time.perf_counter()

    initial_contig_names = set()
    if os.path.exists("correspondance.txt"):
        with open("correspondance.txt") as corr:
            for line in corr:
                _, initial = line.strip("\n").split("\t")
                initial_contig_names.add(initial)
    else:
        with open(genome) as genome_file:
            for line in genome_file:
                if line.startswith(">"):
                    initial_contig_names.add(line[1:].strip("\n"))

    polished_contig_names = set()
    with open("hapog_results/hapog.fasta") as fasta:
        for line in fasta:
            if line.startswith(">"):
                contig_name = line[1:].strip("\n").replace("_polished", "")
                polished_contig_names.add(contig_name)

    with open("hapog_results/hapog.fasta", "a") as out:
        genome_file = open(genome)
        for record in SeqIO.parse(genome_file, "fasta"):
            if record.description.replace("_polished", "") not in polished_contig_names:
                out.write(record.format("fasta"))
        genome_file.close()

    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)
