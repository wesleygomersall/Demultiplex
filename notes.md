# Demultiplex

## 2024-07-25

The files for this assignment are at the following path: 
    ```
    /projects/bgmp/shared/2017_sequencing/
    ```

Some initial data exploration on an interactive node: 

```
$ ls -lah /projects/bgmp/shared/2017_sequencing/
total 46G
drwxrwsr-x+  3 coonrod  is.racs.pirg.bgmp 8.0K Apr 23 13:48 .
drwxrwsr-x+ 42 sdwagner is.racs.pirg.bgmp 8.0K Jan 22  2024 ..
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  20G Jul 30  2018 1294_S1_L008_R1_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp 2.6G Jul 30  2018 1294_S1_L008_R2_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp 2.8G Jul 30  2018 1294_S1_L008_R3_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  21G Jul 30  2018 1294_S1_L008_R4_001.fastq.gz
drwxrws---+  2 coonrod  is.racs.pirg.bgmp 8.0K Jul  1  2022 demultiplexed
-rwxrwxr-x+  1 sdwagner is.racs.pirg.bgmp  631 Aug  9  2021 indexes.txt
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  327 Aug 16  2017 README.txt
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1
GNCTGGCATTCCCAGAGACATCAGTACCCAGTTGGTTCAGACAGTTCCTCTATTGGTTGACAAGGTCTTCATTTCTAGTGATATCAACACGGTGTCTACAA
+
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 2:N:0:1
NCTTCGAC
+
#AA<FJJJ
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 3:N:0:1
NTCGAAGA
+
#AAAAJJF
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1
NTTTTGATTTACCTTTCAGCCAATGAGAAGGCCGTTCATGCAGACTTTTTTAATGATTTTGAAGACCTTTTTGATGATGATGATGTCCAGTGAGGCCTCCC
+
#AAFAFJJ-----F---7-<FA-F<AFFA-JJJ77<FJFJFJJJJJJJJJJAFJFFAJJJJJJJJFJF7-AFFJJ7F7JFJJFJ7FFF--A<A7<-A-7--
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | wc -l
1452986940
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | wc -l
1452986940
$ cat /projects/bgmp/shared/2017_sequencing/indexes.txt 
sample  group   treatment       index   index sequence
1       2A      control B1      GTAGCGTA
2       2B      control A5      CGATCGAT
3       2B      control C1      GATCAAGG
4       2C      mbnl    B9      AACAGCGA
6       2D      mbnl    C9      TAGCCATG
7       2E      fox     C3      CGGTAATC
8       2F      fox     B3      CTCTGGAT
10      2G      both    C4      TACCGGAT
11      2H      both    A11     CTAGCTCA
14      3B      control C7      CACTTCAC
15      3C      mbnl    B2      GCTACTCT
16      3D      mbnl    A1      ACGATCAG
17      3E      fox     B7      TATGGCAC
19      3F      fox     A3      TGTTCCGT
21      3G      both    B4      GTCCTAAG
22      3H      both    A12     TCGACAAG
23      4A      control C10     TCTTCGAC
24      4A      control A2      ATCATGCG
27      4C      mbnl    C2      ATCGTGGT
28      4D      mbnl    A10     TCGAGAGT
29      4E      fox     B8      TCGGATTC
31      4F      fox     A7      GATCTTGC
32      4G      both    B10     AGAGTCCA
34      4H      both    A8      AGGATAGC
```

| File name                     | label         | Read length   | Phred encoding    |
| ---                           | ---           | ---           | ---               |
| 1294_S1_L008_R1_001.fastq.gz  | Forward read  | 101bp         | Phred+33          |
| 1294_S1_L008_R2_001.fastq.gz  | Forward index | 8bp           | Phred+33          |
| 1294_S1_L008_R3_001.fastq.gz  | Reverse index | 8bp           | Phred+33          |
| 1294_S1_L008_R4_001.fastq.gz  | Reverse read  | 101bp         | Phred+33          |

- Encoding is Phred+33, evidence is that there are `-` and `<` characters. 
- There are 1452986940 lines, there are 363246735 records in each file.
- Records in each file have nearly identical headers, the only difference is the read number corresponding to the file.


Made pseudo code [pseudo_wkg.txt](./Assignment-the-first/pseudo_wkg.txt) for demultiplexing algorithm which, given four input FASTQ files (2 with biological reads, 2 with index reads) and the 24 known indexes above, demultiplexes reads by index-pair, outputting:

- One R1 FASTQ file and one R2 FASTQ file per matching index-pair
- Two FASTQ files for non-matching index-pairs (index-hopping): R1 and R2
- Two additional FASTQ files when one or both index reads are unknown or low quality (do not match the 24 known indexes [this includes indexes with 'N's in them] or do not meet a quality score cutoff)
- Add the sequence of the index-pair to the header of BOTH reads in all of your FASTQ files for all categories (e.g. add “AAAAAAAA-CCCCCCCC” to the end of headers of every read pair that had an index1 of AAAAAAAA and an index2 of CCCCCCCC; this pair of reads would be in the unknown category as one or both of these indexes do not match the 24 known indexes).

Additionally, your algorithm should report:
- The number of read-pairs with properly matched indexes (per index-pair)
- The number of read pairs with index-hopping observed
- The number of read-pairs with unknown index(es)
- You should strive to report values for each possible pair of indexes (both swapped and dual matched). 

While developing pseudocode, determine high level functions: 
- Function headers (name and parameters)
- Description/doc string – What does this function do?
- Test examples for individual functions
- Return statement

## 2024-07-29

### (Assignment 1: Part 1) – Quality Score Distribution per-nucleotide
Create a python script [mean_qual_fq.py](./Assignment-the-first/mean_qual_fq.py) to generate per base mean of quality scores for each file: read1, read2, index1, and index2. This is the same thing I did in part 1 of PS4 in Bi621 (no variance necessary in this plot). I cannot use a 2-dimensional array in numpy to calculate this because there wont be enough memory!  
Link all 4 plots generated into markdown file `./Assignment-the-first/Answers.md`. 

I made a slurm script [demultiplexp1.sh](./Assignment-the-first/demultiplexp1.sh) to run this job. [Output](./Assignment-the-first/slurm-7776055.out)

Create conda environment with matplotlib: 
```
conda create --name bgmp_py.mplib --clone bgmp_py312
conda activate bgmp_py.mtplib
conda install matplotlib
```

Create unit tests: 4 fastq files input (at least 1 entry of each: dual matched, index-hopped, unknown index), at least 6 fastq files output (R1.fq and R2.fq for each type of entry). Differently-barcoded reads should end up in different fastq files. 

Created python functions `check_n` and `reverse_complement` to import into my script. I added it to [bioinfo.py](./Assignment-the-first/bioinfo.py).  

## 2024-07-30

Leslie says we should not be using these functions in bioinfo.py. Oops. I will just copy and paste them into my new script and comment them out of bioinfo.py for now.  

Begin scripting for Assignment-the-third. [demultiplex.py](./Assignment-the-third/demultiplex.py)

This is the command to test out the code. 

```
./Assignment-the-third/demultiplex.py --read1 TEST-input_FASTQ/testwkg_R1.fq --read2 TEST-input_FASTQ/testwkg_R2.fq --read3 TEST-input_FASTQ/testwkg_R3.fq --read4 TEST-input_FASTQ/testwkg_R4.fq -i TEST-input_FASTQ/testwkg_indexes.txt -c 30
```

This code demultiplexes the test files. It will not work on real data in its current state because it is not reading compressed files. 
The code to open these files is commented to swap between gzip.open and open

Comment out the `open` for `gzip.open` lines and try the series of commands:

	gzip ~/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/testwkg_R*.fq
	./Assignment-the-third/demultiplex.py --read1 TEST-input_FASTQ/testwkg_R1.fq.gz --read2 TEST-input_FASTQ/testwkg_R2.fq.gz --read3 TEST-input_FASTQ/testwkg_R3.fq.gz --read4 TEST-input_FASTQ/testwkg_R4.fq.gz -i TEST-input_FASTQ/testwkg_indexes.txt -c 30

And it works! I will add a print out of the count dictionary (the large one) and this script should be just about ready for the full data set. 

## 2024-07-31

Pylance is going crazy. I will just use `# type: ignore` to silence them.  

Change naming scheme: Use the barcode sequences themselves, not the shorter names.  

Also need to change to not writing to compressed files.  

Change the expected output file names to match the new file naming system. 

Create [slurm script](./Assignment-the-third/demultiplex_slurm.sh) for running demultiplex.py

First run of script, I use cutoff 35. I intend on reconsidering this value and re-running this later. See [slurm-7797951.out](./slurm-7797951.out). 

```
Hopped read count: 352378
Matched-index read count: 239750646
Unknown-indexed read count: 123143711
```

## 2024-08-06

Probility of an error in a given sequence is the sum of phred scores. So given average quality score, multiply by length to get an estimate of errors in the sequence.

Divide 0.05 by 8, and the closest quality score associated with .00625 is 22.

I am going to re-run the bash script to demultiplex with a more sensible cutoff of 22. See [slurm-7923319.out](./slurm-7923319.out).

```
Hopped read count: 662066
Matched-index read count: 329566234
Unknown-indexed read count: 33018435
```
