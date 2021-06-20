# 5R_SMURF2
SMURF2 reconstruction of sequencing results based on the 5R protocol

5R package – A software package accompanying a novel method of bacterial 16S rRNA profiling based on the 5R protocol. Five short regions along the 16S rRNA are amplified in multiplex and subsequently sequenced by an Illumina machine. Resulting reads from these regions are combined into a coherent solution using our Short MUltiple Regions Framework (SMURF2) algorithm enabling identificaiton of both known and undocumented bacteria.

## Installation and Requirments
Matlab
Python version = 2.7
python modules:
- sys, os, glob, shutil, time, logging
- optparse
- ast
- gmpy2
- itertools
- math
- numpy
- random
- difflib
- Bio 
- pandas

Download the 5R_SMURF2 package

## Usage
- change the directory to the code directory. 
- load matlab (not mandatory, for taxonomy)
- activate pyton environment
- python smurf2.py WORKDIR -1 FULL_PATH_FASTQ_FW -2 FULL_PATH_FASTQ_BW -f PATH_TO_5R_DATABASE -l kmer_len -r AMPLIFIED_REGIONS 

where the inputs are:

WORKDIR - the (full) path to SMURF2 working dirctory

FULL_PATH_FASTQ_FW - the (full) path to the fastq file of forward sequencing file of samples profiled by the 5R primers.

FULL_PATH_FASTQ_BW - the (full) path to the fastq file of backward sequencing file of samples profiled by the 5R primers.

kmer_len - the length of the k-mer to be applied for reconstruction (126 in the attached example).

AMPLIFIED_REGIONS - the amount of amplified regions (5)

## Output

- new_bacteria_*.fasta - Fasta files with the undocumented bacteria found by SMURF2 (file per amplified region).
- SMURF2_results.csv - A table with the abundance and the full sequences of the reconstructed bacteria.
- 
# optional:
- taxonomy_smurf2.txt - A table of species' relative abundances (rows), in case of undocumented bacteria, the output will be the closest one in the database. 
- SMURF2_results.mat - A matlab file contains the ids of the amplified bacteria. 

## Example
The directory "example_fastq" contains 5R sequencing files of in silico mixture. To perform SMURF2 reconstruction of these samples use the following command:

## Contact us
For questions please email:

Vered Bar-Moshe: veredhir@gmail.com

Noam Shental: shental@openu.ac.il


