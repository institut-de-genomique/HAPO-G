# Hapo-G - Haplotype-Aware Polishing of Genomes

Hapo-G (pronounced like apogee) is a tool that aims to improve the quality of genome assemblies by polishing the consensus with accurate reads.

Biorxiv preprint : [link](https://www.biorxiv.org/content/10.1101/2020.12.14.422624v1 "Hapo-G Biorxiv preprint")

In case of troubles when using or installing the software, please open up an issue by clicking [here](https://github.com/institut-de-genomique/Hapo-G/issues/new "Github issue page").


## Dependencies

Hapo-G depends on some software and libraries:
- GCC and G++ (Hapo-G has been tested with GCC 4.9.2 and GCC 7.3.0)
- Autoconf with minimum 2.69 version (to build HTSlib)
- Python3 (minimum version 3.6)
- [HTSlib](https://github.com/samtools/htslib "HTSlib github")
- [BioPython](https://biopython.org/wiki/Download "BioPython")
- [BWA](https://github.com/lh3/bwa "BWA")
- [Samtools](https://github.com/samtools/samtools "Samtools")


## Installation
First, clone this repository:
```
git clone https://github.com/institut-de-genomique/HAPO-G hapog
```

If htslib is already installed on your system, go to the next point `Build with existing htslib`. If you want Hapo-G to download and install htslib for you, go to the point `Build with a new htslib install`

### Build with existing htslib
Building with an existing htslib ensures that Hapo-G and Samtools are using the same version of the library and should reduce compatibility issues. To build with an existing htslib, do:
```
cd hapog
bash build.sh -l path_to_htslib
```
If samtools is already installed on your system at `/home/user/samtools`, htslib is probably installed at `/home/user/samtools/htslib`.

### Build with a new htslib
Hapo-G can download and compile htslib for you, to do so, please run:
```
cd hapog
bash build.sh
```

If everything went as expected, a binary of Hapo-G was created in `build/` and a symlink was created in the `bin/` folder


## Using Hapo-G
Before running Hapo-G, you should make sure that BWA and Samtools are in your `$PATH`:
```
which bwa
which samtools
```

### Standard pipeline
You can launch Hapo-G by using the Python3 script in its root directory:
```
python3 HAPOG_ROOT/hapog.py \
  --genome assembly.fasta \   # Fasta file of the genome to polish
  --pe1 R1.fastq.gz \         # Illumina R1 reads in .fastq or .fastq.gz, can be given multiple times
  --pe2 R2.fastq.gz \         # Illumina R2 reads in .fastq or .fastq.gz, can be given multiple times
  -o polishing \              # Output directory
  -t 36 \                     # Number of threads to use
  -u                          # Include unpolished sequences in the output
```

### Skipping the mapping step
The mapping step can be skipped if a sorted BAM file is provided via the `-b` switch. Please verify that your fasta headers don't contain any non-alphanumerical characters (`-`and `_`are accepted) before launching Hapo-G.
A typical command line with a bam file would look like this:
```
python3 HAPOG_ROOT/hapog.py \
  --genome assembly.fasta \   # Fasta file of the genome to polish
  -b mapping.sorted.bam       # Sorted BAM file
  -o polishing \              # Output directory
  -t 36 \                     # Number of threads to use
  -u                          # Include unpolished sequences in the output
```

## Acknowledgements
Some Cmake files have been taken and/or modified from several projects. We would like to thank:
- [panguangze](https://delta.cs.cityu.edu.hk/gzpan2) for his/her `FindHTSLIB.cmake` library
- [L. Kärkkäinen](https://github.com/Tronic) for his `LibFindMacros.cmake` library
