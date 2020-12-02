# HAPO-G - Haplotype-Aware Polishing of Genomes

HAPO-G (pronounced like apogee) is a tool that aims to improve the quality of genome assemblies by polishing the consensus with accurate reads.

Biorxiv preprint : [link](https://www.biorxiv.org/ "HAPO-G Biorxiv preprint")

In case of troubles when using or installing the software, please open up an issue by clicking [here](https://github.com/institut-de-genomique/HAPoG/issues/new "Github issue page").


## Dependencies

HAPO-G depends on some software and libraries:
- GCC and G++ (HAPoG has been tested with GCC 4.9.2 and GCC 7.3.0)
- Python3 (minimum version 3.6)
- [HTSlib](https://github.com/samtools/htslib "HTSlib github") (Automatically downloaded and built with HAPoG)
- [BioPython](https://biopython.org/wiki/Download "BioPython")
- [BWA](https://github.com/lh3/bwa "BWA")
- [Samtools](https://github.com/samtools/samtools "Samtools")


## Installation
First, clone this repository:
```
git clone https://github.com/institut-de-genomique/HAPO-G hapog
```

Go into the created directory and run the build script:
```
cd hapog
bash build.sh
```

If everything went as expected, a binary of HAPoG was created in `build/` and a symlink was created in the `bin/` folder


## Using HAPO-G
Before running HAPO-G, you should make sure that BWA and Samtools are in your `$PATH`:
```
which bwa
which samtools
```

Then, you can launch HAPO-G by using the Python3 script in its root directory:
```
python3 HAPOG_ROOT/hapog.py \
  --genome assembly.fasta \   # Fasta file of the genome to polish
  --pe1 R1.fastq.gz \         # Illumina R1 reads in .fastq or .fastq.gz, can be given multiple times
  --pe2 R2.fastq.gz \         # Illumina R2 reads in .fastq or .fastq.gz, can be given multiple times
  -o polishing \              # Output directory
  -t 36                       # Number of threads to use
```


