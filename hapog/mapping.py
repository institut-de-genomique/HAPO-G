import os
import subprocess
import time
import warnings


def launch_PE_mapping(genome, pe1, pe2, threads, samtools_memory):
    print("\nGenerating bwa index...", flush=True)
    cmd = ["bwa", "index", genome]

    start = time.perf_counter()
    with open("cmds/bwa_index.cmds", "w") as cmd_file:
        print(" ".join(cmd), flush=True, file=cmd_file)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _ = subprocess.run(
                cmd,
                stdout=open("logs/bwa_index.o", "w"),
                stderr=open("logs/bwa_index.e", "w"),
                check=True,
            )
    except Exception as e:
        print("\nERROR: Couldn't index genome", flush=True)
        print(e)
        exit(1)

    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)

    print("\nLaunching mapping on genome...", flush=True)
    cmd = "bash -c 'bwa mem -t %s %s " % (threads, genome)

    streamer = "cat"
    if pe1[0].endswith(".gz"):
        streamer = "zcat"

    if len(pe1) > 1:
        cmd += f"<({streamer}"
        for pe in pe1:
            cmd += " %s" % pe
        cmd += ") "
        cmd += f"<({streamer}"
        for pe in pe2:
            cmd += " %s" % pe
        cmd += ") 2> logs/bwa_mem.e"
    else:
        cmd += pe1[0] + " " + pe2[0] + " 2> logs/bwa_mem.e"

    cmd += f" | samtools sort -m {samtools_memory} -@ {threads} -o bam/aln.sorted.bam - 2> logs/samtools_sort.e'"

    start = time.perf_counter()

    with open("cmds/mapping.cmds", "w") as cmd_file:
        print(cmd, flush=True, file=cmd_file)

    return_code = os.system(cmd)
    if return_code != 0:
        print(f"Error in bwa mem and samtools sort, return code: {return_code}")
        print(f"Faulty command: {cmd}")
        exit(1)
    else:
        print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)

    index_bam()


def launch_LR_mapping(genome, long_reads, threads, samtools_memory):
    print("\nLaunching mapping on genome...", flush=True)
    cmd = f"bash -c 'minimap2 -t {threads} -a --secondary=no -x map-pb {genome} {long_reads} 2> logs/minimap2.e"
    cmd += f" | samtools sort -m {samtools_memory} -@ {threads} -o bam/aln.sorted.bam - 2> logs/samtools_sort.e'"

    start = time.perf_counter()

    with open("cmds/mapping.cmds", "w") as cmd_file:
        print(cmd, flush=True, file=cmd_file)

    return_code = os.system(cmd)
    if return_code != 0:
        print(f"Error in minimap2 and samtools sort, return code: {return_code}")
        print(f"Faulty command: {cmd}")
        exit(1)
    else:
        print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)

    index_bam()


def remove_secondary_alignments(bam, output_dir):
    print("\nRemoving secondary alignments from BAM file...", flush=True)
    cmd = ["samtools", "view", "-b", "-h", "-F", "0x900", bam]

    start = time.perf_counter()
    with open(f"{output_dir}/cmds/samtools_view.cmds", "w") as cmd_file:
        print(" ".join(cmd), flush=True, file=cmd_file)
        
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _ = subprocess.run(
                cmd,
                stdout=open(f"{output_dir}/bam/aln.sorted.bam", "w"),
                stderr=open(f"{output_dir}/logs/samtools_view.e", "w"),
                check=True,
            )
    except Exception as e:
        print("\nERROR: Couldn't remove secondary alignments")
        print(e)
        exit(1)
    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)


def index_bam():
    print("\nIndexing the BAM file...", flush=True)
    cmd = ["samtools", "index", "bam/aln.sorted.bam"]

    start = time.perf_counter()
    with open("cmds/samtools_index.cmds", "w") as cmd_file:
        print(" ".join(cmd), flush=True, file=cmd_file)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _ = subprocess.run(
                cmd,
                stdout=open("logs/samtools_index.o", "w"),
                stderr=open("logs/samtools_index.e", "w"),
                check=True,
            )
    except Exception as e:
        print("\nERROR: Couldn't index bam file")
        print(e)
        exit(1)
    print(f"Done in {int(time.perf_counter() - start)} seconds", flush=True)
