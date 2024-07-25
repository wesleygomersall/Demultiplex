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
```


`R2` file is FORWARD? 
`R3` file is REVERSE? 

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Forward read?  | 101bp | Phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | Forward index? | 8bp | Phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | Reverse index? | 8bp | Phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | Reverse read?  | 101bp | Phred+33 |

Encoding is Phred+33, evidence is that there are `-` and `<` characters. 
There are 1452986940 lines, there are 363246735 records in each file.
Records in each file have nearly identical headers, the only difference is the read number corresponding to the file.



## todo

Part 1 – Quality Score Distribution per-nucleotide

Generate a per base distribution of quality scores for read1, read2, index1, and index2.
Average the quality scores at each position for all reads and generate a per nucleotide mean distribution as you did in part 1 of PS4 in Bi621. 
(NOTE! Do NOT use the 2D array strategy from PS9 - you WILL run out of memory!)

In file `./Assignment-the-first/Answers.md`
- Turn in the 4 histograms.
- What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.
- How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command(s) you used. CHALLENGE: use a one-line command)