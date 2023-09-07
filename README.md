# Hapo-G - Haplotype-Aware Polishing of Genomes

Hapo-G (pronounced like apogee) is a tool that aims to improve the quality of genome assemblies by polishing the consensus with accurate reads.

Publication : [link](https://academic.oup.com/nargab/article/3/2/lqab034/6262629 "Hapo-G publication")

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
### Installation with conda
```
conda create -n hapog
conda activate hapog
conda install -c bioconda hapog
```

### Installing from Github
First, clone this repository:
```
git clone https://github.com/institut-de-genomique/HAPO-G hapog
```

If htslib is already installed on your system, go to the next point `Build with existing htslib`. If you want Hapo-G to download and install htslib for you, go to the point `Build with a new htslib install`

#### Build with existing htslib
Building with an existing htslib ensures that Hapo-G and Samtools are using the same version of the library and should reduce compatibility issues. To build with an existing htslib, do:
```
cd hapog
bash build.sh -l path_to_htslib
```
If samtools is already installed on your system at `/home/user/samtools`, htslib is probably installed at `/home/user/samtools/htslib`.

#### Build with a new htslib
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

**NOTE**: If you installed Hapo-G using conda, you can invoke it by directly running `hapog`.

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

## Output files
#### `hapog_results/hapog.fasta`
The corrected sequences. Hapo-G will parse the read alignments to the genome and focus on phasing errors (i.e the assembly switched from one haplotype to the other) and base errors (insertions, deletions, mismatches) that may be related or not to phasing errors. Remember to include the `-u` flag to tell Hapo-G to output sequences with no reads mapped and thus could not need to be changed.

Hapo-G will not add any new contigs or scaffolds to the assembly if, as an example, one of the haplotype is missing in the input assembly file. Instead, it will correct the haplotype that is present in the input file and output a corrected version of the sequence that is phased as best as we could with the data at hand. 

As an example, let’s consider the following heterozygous genome:
```text
maternal hap.: ACCGTTA
paternal hap.: ATCGTGA
```
If the assembler outputted a version with one phasing error (the C in 2nd position is associated with a G in 6th position) and one deletion in 4th position:
```text
assembly: ACC-TGA
```
Then, if Hapo-G was able to correct all the errors, the `hapog.fasta` file will contain:
```text
hapog.fasta: ACCGTTA
```

## `hapog_results/hapog.changes`
This file is a tabulated file that gives more information on what Hapo-G did to the input assembly. It has eight columns that show:
Name of the input sequence where the change was made
Position in the input sequence where the change was made
Nucleotide(s) at the position in column 2
Nucleotide(s) that will replace the nucleotide(s) shown in column 3
Name of the read that is used as the current template. In the Hapo-G algorithm, we try to follow a read for as long as possible to not switch haplotypes. If an error is found in the template read, we switch to a different read of the same haplotype
`hetero` if the change is only present in one of the two possible haplotypes (i.e a phasing error). `homo` if the change is present in both haplotypes
Ratio of reads from the same haplotypes as the template read that validate the changes
Ratio of reads from both haplotypes that validate the changes

Here is an examples:
```text
Contig_1	1000	ref=A	read=TA	readname=read_2	hetero	ratio1=0.7419	ratio2=0.4237
Contig_1  2000  ref=T read=G  readname=read_2 homo    ratio1=0.8142 ratio2=0.8323
```
We can see that on the contig `Contig_1`, Hapo-G found a phasing error (`hetero` on the first line) and replaced a A at position 1000 by TA. This change was validated by 74.19% of reads of the same haplotype as the template read (ratio1) and by 42.37% of reads if we do not discriminate on which haplotype they belong to. It also found a mismatch at position 2000 (`homo`) and replaced a T by a G. This change was validated by 83% of reads of no matter which haplotype (ratio2).


## Acknowledgements
Some Cmake files have been taken and/or modified from several projects. We would like to thank:
- [panguangze](https://delta.cs.cityu.edu.hk/gzpan2) for his/her `FindHTSLIB.cmake` library
- [L. Kärkkäinen](https://github.com/Tronic) for his `LibFindMacros.cmake` library
