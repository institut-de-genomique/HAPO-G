from Bio import SeqIO
from collections import defaultdict

import subprocess


def extract_subx_coverage_regions(coverage_file, min_coverage):
	with open("sup_%s_coverage_regions.txt" % (min_coverage), "w") as out:
		for line in open(coverage_file):
			_, _, coverage = line.rstrip("\n").split("\t")
			coverage = int(coverage)
			if coverage > min_coverage:
				out.write(line)


def replace_high_cov_regions_by_N(input_genome, min_coverage):
	high_cov_dict = defaultdict(list)
	print("\tLoading sup-%s coverage regions..." % (min_coverage), flush=True)
	for line in open("sup_%s_coverage_regions.txt" % (min_coverage)):
		line = line.split("\t")
		high_cov_dict[line[0]].append(int(line[1]))

	genome_dict = {}
	print("\tLoading genome...", flush=True)
	for record in SeqIO.parse(open(input_genome), "fasta"):
		genome_dict[str(record.id)] = list(str(record.seq))

	print("\tMasking high-coverage regions...", flush=True)
	for contig in high_cov_dict:
		for pos in high_cov_dict[contig]:
			genome_dict[contig][pos-1] = "N"

	print("\tWriting masked genome...", flush=True)
	with open("genome_masked.fasta", "w") as out:
		for contig in genome_dict:
			out.write(">%s\n%s\n" % (contig, "".join(genome_dict[contig])))


def mask_genome(input_bam, input_genome):
	print("\nLaunching Bedtools...", flush=True)
	coverage_file = "coverage.txt"
	min_coverage = 5

	cmd = ["bedtools", "genomecov", "-ibam", input_bam, "-d", "-g", "genome.txt"]
	print(" ".join(cmd), flush=True)
	_ = subprocess.run(cmd, 
	    stdout = open(coverage_file, "w"), 
	    stderr = open("bedtools.e", "w"))

	print("Parsing coverage file...", flush=True)
	extract_subx_coverage_regions(coverage_file, min_coverage)

	print("\nMasking genome...")
	replace_high_cov_regions_by_N(input_genome, min_coverage)