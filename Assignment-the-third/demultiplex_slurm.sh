#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=wesg@uoregon.edu
#SBATCH --mail-type=END

# this script should be ran from the main git repo "Demultiplex"

# ./Assignment-the-third/demultiplex.py \ # for test run:
	# --read1 TEST-input_FASTQ/testwkg_R1.fq.gz \
	# --read2 TEST-input_FASTQ/testwkg_R2.fq.gz \
	# --read3 TEST-input_FASTQ/testwkg_R3.fq.gz \
	# --read4 TEST-input_FASTQ/testwkg_R4.fq.gz \
	# -i TEST-input_FASTQ/testwkg_indexes.txt \
	# -c 30

./Assignment-the-third/demultiplex.py \ # for final run:
	--read1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
	--read2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
	--read3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
	--read4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
	-i /projects/bgmp/shared/2017_sequencing/indexes.txt \
	-c 35 \
