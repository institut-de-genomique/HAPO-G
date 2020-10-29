# HAPoG - Haplotype-Aware Polishing of Genomes

HAPoG (pronounced like apogee) is a tool that aims to improve the quality of genome assemblies by polishing the consensus with accurate reads.

Biorxiv preprint : [link](https://www.biorxiv.org/ "BiSCoT Biorxiv preprint")

In case of troubles when using or installing the software, please open up an issue by clicking [here](https://github.com/institut-de-genomique/HAPoG/issues/new "Github issue page").


## Dependencies

HAPoG depends on some software and libraries:
- GCC and G++ (HAPoG has been tested with GCC 4.9.2 and GCC 7.3.0)
- Python3 (minimum version 3.6)
- [HTSlib](https://github.com/samtools/htslib "HTSlib github") (Automatically downloaded and built with HAPoG)
- [BioPython](https://biopython.org/wiki/Download "BioPython")
- [BWA](https://github.com/lh3/bwa "BWA")
- [Samtools](https://github.com/samtools/samtools "Samtools")


## Installation
First, clone this repository:
```
git clone https://github.com/institut-de-genomique/HAPoG
```

Go into the created directory and run the build script:
```
cd HAPoG
bash install.sh
```

If everything went as expected, a binary of HAPoG was created in `build/` and a symlink was created in the `bin/` folder
