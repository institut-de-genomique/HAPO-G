import distutils.command.build

from setuptools import setup
from os import path

# Override build command
class BuildCommand(distutils.command.build.build):
    def initialize_options(self):
        distutils.command.build.build.initialize_options(self)
        self.build_base = 'package_build'

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md")) as f:
    long_description = f.read()

setup(
    name="hapog",
    packages=["hapog"],
    version="1.3.5",
    license="CeCILL",
    description="Haplotype-Aware Polishing of Genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=["Jean-Marc Aury", "Benjamin Istace"],
    author_email=["jmaury@genoscope.cns.fr", "bistace@genoscope.cns.fr"],
    url="https://github.com/institut-de-genomique/HAPO-G",
    download_url="https://github.com/institut-de-genomique/HAPO-G",
    keywords=[
        "bioinformatics",
        "genomics",
        "genome",
        "polishing",
    ],
    install_requires=["biopython"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        "Programming Language :: Python :: 3.6",
    ],
    entry_points={"console_scripts" : ["hapog = hapog.cli:main"]},
    cmdclass={"build": BuildCommand},
)
