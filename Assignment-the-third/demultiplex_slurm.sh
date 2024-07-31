#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=wesg@uoregon.edu
#SBATCH --mail-type=END

# this script should be ran from the main git repo "Demultiplex"

# for test run:
./Assignment-the-third/demultiplex.py --read1 TEST-input_FASTQ/testwkg_R1.fq.gz --read2 TEST-input_FASTQ/testwkg_R2.fq.gz --read3 TEST-input_FASTQ/testwkg_R3.fq.gz --read4 TEST-input_FASTQ/testwkg_R4.fq.gz -i TEST-input_FASTQ/testwkg_indexes.txt -c 30

# for final run: 

# ./Assignment-the-third/demultiplex.py --read1 TEST-input_FASTQ/testwkg_R1.fq.gz --read2 TEST-input_FASTQ/testwkg_R2.fq.gz --read3 TEST-input_FASTQ/testwkg_R3.fq.gz --read4 TEST-input_FASTQ/testwkg_R4.fq.gz -i TEST-input_FASTQ/testwkg_indexes.txt -c 30
