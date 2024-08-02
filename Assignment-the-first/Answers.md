# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:




| File name                     | label         | Read length   | Phred encoding    |
| ---                           | ---           | ---           | ---               |
| 1294_S1_L008_R1_001.fastq.gz  | Forward read  | 101bp         | Phred+33          |
| 1294_S1_L008_R2_001.fastq.gz  | Forward index | 8bp           | Phred+33          |
| 1294_S1_L008_R3_001.fastq.gz  | Reverse index | 8bp           | Phred+33          |
| 1294_S1_L008_R4_001.fastq.gz  | Reverse read  | 101bp         | Phred+33          |

Per-base NT distribution
    1. [meanqual_Read1.png](./meanqual_Read1.png)
    2. [meanqual_Index1.png](./meanqual_Index1.png)
    3. [meanqual_Index2.png](./meanqual_Index2.png)
    4. [meanqual_Read2.png](./meanqual_Read2.png)
    
My first run of the demultiplexing algorithm [../slurm-7797951.out](../slurm-779751.out) I used 35 as a cutoff. 
I know that is bad. I hadn't given it much thought. Now that I am considering it I am going to choose based on this criteria:
- I want to allow a maximum 5% chance that the barcode has an error.  

I know this: Probability of an error in a given sequence is the sum of phred scores. So given average quality score, multiply by length to get an estimate of errors in the sequence.  
Divide 0.05 by 8, and the closest quality score associated with .00625 is 22.  

I will use 22 on the next run of this program. 
If I wanted to allow 1% chance of error then I would use a cutoff of 29. 35 is overkill, that is demanding no more than 0.256% error rate in barcodes.

To determine how many indexes have 'N' base call I used the following commands.
```
[wesg@n0349 2017_sequencing]$ zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep 'N' | wc -l
[wesg@n0349 2017_sequencing]$ zcat 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep 'N' | wc -l
```
and got 3976613 and 3328051 respectively. This is not extremely helpful, as all I know from this is are upper (sum of those two) and lower bounds (larger of the two numbers) about how many read pairs have 'N' base calls. It would be very useful to know exactly how many entries I will categorize as unknown strictly because of 'N' base calls in the index sequences. 

I can do this with the following command: 

```
[wesg@n0349 2017_sequencing]$ zcat 1294_S1_L008_R2_001.fastq.gz 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep -n 'N' | cut --delim
iter=":" -f 1 | sort | uniq | wc -l
7304664
```

## Part 2
1. Define the problem  
Given four files: two reads and two indexes, separate the read files out into multiple output files based on index sequences.  
There are three categories of reads: Unknown index reads, hopped barcode reads, and (many) dual-indexed reads. 
2. Describe output  
The output for my program is going to be 4 + indexes * 2 separate files, each type of read will get a R1 and an R2 file containing the reads from the original R1 and R4 files respectively. 
3. Test files: 
- [../TEST-input_FASTQ/testwkg_R1.fq](../TEST-input_FASTQ/testwkg_R1.fq)
- [../TEST-input_FASTQ/testwkg_R2.fq](../TEST-input_FASTQ/testwkg_R2.fq)
- [../TEST-input_FASTQ/testwkg_R3.fq](../TEST-input_FASTQ/testwkg_R3.fq)
- [../TEST-input_FASTQ/testwkg_R4.fq](../TEST-input_FASTQ/testwkg_R4.fq)

- [../TEST-output_FASTQ/testwkg_CGATCGAT_R1.fq](../TEST-output_FASTQ/twkg_CGATCGAT_R1.fq)
- [../TEST-output_FASTQ/testwkg_CGATCGAT_R2.fq](../TEST-output_FASTQ/testwkg_CGATCGAT_R2.fq)
- [../TEST-output_FASTQ/testwkg_CTAGCTCA_R1.fq](../TEST-output_FASTQ/testwkg_CTAGCTCA_R1.fq)
- [../TEST-output_FASTQ/testwkg_CTAGCTCA_R2.fq](../TEST-output_FASTQ/testwkg_CTAGCTCA_R2.fq)
- [../TEST-output_FASTQ/testwkg_GATCAAGG_R1.fq](../TEST-output_FASTQ/testwkg_GATCAAGG_R1.fq)
- [../TEST-output_FASTQ/testwkg_GATCAAGG_R2.fq](../TEST-output_FASTQ/testwkg_GATCAAGG_R2.fq)
- [../TEST-output_FASTQ/testwkg_GTAGCGTA_R1.fq](../TEST-output_FASTQ/testwkg_GTAGCGTA_R1.fq)
- [../TEST-output_FASTQ/testwkg_GTAGCGTA_R2.fq](../TEST-output_FASTQ/testwkg_GTAGCGTA_R2.fq)
- [../TEST-output_FASTQ/testwkg_hopped_R1.fq](../TEST-output_FASTQ/testwkg_hopped_R1.fq)
- [../TEST-output_FASTQ/testwkg_hopped_R2.fq](../TEST-output_FASTQ/testwkg_hopped_R2.fq)
- [../TEST-output_FASTQ/testwkg_unk_R1.fq](../TEST-output_FASTQ/testwkg_unk_R1.fq)
- [../TEST-output_FASTQ/testwkg_unk_R2.fq](../TEST-output_FASTQ/testwkg_unk_R2.fq)

4. Pseudocode: [./pseudo_wkg.txt](./pseudo_wkg.txt)
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement  

```
def rev_complement(seq: str)-> revcomp: str:
	"""Returns the reverse complement of a DNA sequence"""
	input: "ACTGTGACA"
        expected output: "TGTCACAGT"
        return revcomp
def check_n_bases(seq: str) -> nbases: bool:
        """checks if there is an 'N' in a DNA sequence""""
        input: "ACTGTGACA"
        expected output: False
        input 2: "ACTGTCANA"
        expected output 2: True
        return nbases
def mean_qual(quals: str) -> meanq: float:
        """calculates mean quality of a quality score string"""
        input: "FFHH"
        expected output: 38.0
        return meanq
```
