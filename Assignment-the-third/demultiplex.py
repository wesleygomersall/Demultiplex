#!/usr/bin/env python

import argparse
import gzip
import bioinfo

def get_args():
    parser = argparse.ArgumentParser(description="Program to demultiplex illumina reads. Takes in all 4 illumina files R1-R4 and file containing indexes. Output fastqs separated by index (for dual-indexed samples) as well as files for unknown indexes and for hopped indexes.")
    parser.add_argument("--read1", help="Illumina R1 fastq file.", type=str, required=True) 
    parser.add_argument("--read2", help="Illumina R2 fastq file.", type=str, required=True) 
    parser.add_argument("--read3", help="Illumina R3 fastq file.", type=str, required=True) 
    parser.add_argument("--read4", help="Illumina R4 fastq file.", type=str, required=True) 
    parser.add_argument("-i", "--indexes", help="Indexes to use for checking dual-indexes", type=str, required=True) 
    parser.add_argument("-c", "--cutoff", help="Lowest allowed mean quality of indexes", type=int, required=True) 
    return parser.parse_args()

READ1 = get_args().read1
READ2 = get_args().read2
READ3 = get_args().read3
READ4 = get_args().read4
INDEXES = get_args().indexes 
CUTOFF = get_args().cutoff 

def reverse_complement(sequence: str, RNAflag: bool = False) -> str:
    '''Returns the reverse complement of a DNA or an RNA string (default DNA).
    Input and output are both 5' -> 3'.
    '''
    assert bioinfo.validate_base_seq(sequence, RNAflag)
    DNA_comp = {'A':'T', 'a':'t',
                'T':'A', 't':'a',
                'C':'G', 'c':'g',
                'G':'C', 'g':'c',
                'N':'N', 'n':'n'}
    RNA_comp = {'A':'U', 'a':'u',
                'U':'A', 'u':'a',
                'C':'G', 'c':'g',
                'G':'C', 'g':'c',
                'N':'N', 'n':'n'}
    reversecomp: str = ""
    if RNAflag == False:
        for index in range(len(sequence)):
            reversecomp = reversecomp + DNA_comp[sequence[-1*(index + 1)]]
    elif RNAflag:
        for index in range(len(sequence)):
            reversecomp = reversecomp + RNA_comp[sequence[-1*(index + 1)]]
    return reversecomp

assert reverse_complement("ATCG") == "CGAT"
assert reverse_complement("ATCGAAA") == "TTTCGAT"
assert reverse_complement("AUCGAAA", True) == "UUUCGAU" #RNA
assert reverse_complement("ATCGCGTGCTCTCTCTTCTCT") == "AGAGAAGAGAGAGCACGCGAT"

def check_n(seq: str) -> bool: 
    '''Bool check if a string contains 'n' or 'N'.
    '''
    bases = set(seq.upper())
    if {'N'}.issubset(bases): 
        return True
    else:
        return False

assert check_n("ACTCGCT") == False
assert check_n("ACTCGNT") == True
assert check_n("ACNCGCT") == True

filenamingdict = dict() 
indexsequencelist = []

with open(INDEXES, 'r') as fin:
    for i, line in enumerate(fin): 
        if i == 0 :
            continue
        index = line.strip('\n').split()[3] # the 4th column from file holds index names
        sequence = line.strip('\n').split()[4] # the final column from file holds index sequences
        indexsequencelist.append(sequence) # simple list of sequences to loop through for dictionaries
        filenamingdict[sequence] = index # add entry to dictionary: key is sequence, value is the name

allpairsdict = dict() 
matchedindexlist = []

for i in indexsequencelist:
    for j in indexsequencelist:
        pair = (i, j) #tuple of indexes as they appear in R2 and (reverse complemented in) R3
        allpairsdict.setdefault(pair, 0) # use this dictionary to increment found (hopped and matched) pairs. 
        if i == j: 
            matchedindexlist.append(pair) #list of matched indexes

matchedcount: int = 0
hoppedcount: int = 0
unkcount: int = 0 

filenames = [] # create a list of filenames
for i in indexsequencelist: #for index-matched files, create filenames
    filenames.append(f"./output/{filenamingdict[i]}_R1.fq.gz")
    filenames.append(f"./output/{filenamingdict[i]}_R2.fq.gz")
for i in ["unk", "hopped"]: #need to make files for unknown reads and for hopped reads
    filenames.append(f"./output/{i}_R1.fq.gz")
    filenames.append(f"./output/{i}_R2.fq.gz")

filedata = {filename: gzip.open(filename, 'wt') for filename in filenames} #open all the files to write to

with gzip.open(READ1, 'rt') as r1, gzip.open(READ2, 'rt') as r2, gzip.open(READ3, 'rt') as r3, gzip.open(READ4, 'rt') as r4: #open to read all 4 files (.gz)
#with open(READ1, 'r') as r1, open(READ2, 'r') as r2, open(READ3, 'r') as r3, open(READ4, 'r') as r4: #open to read all 4 files (unzipped)
    linenum: int = 0



    while True: 
        unknown: bool = False
        matched: bool = False
        hopped: bool = False

        r1line, r2line, r3line, r4line = r1.readline().strip(), r2.readline().strip(), r3.readline().strip(), r4.readline().strip()

        if r1line == "": # end of file
            break # stop loop
        
        if linenum % 4 == 0: #header line
            r1head, r2head, r3head, r4head = r1line, r2line, r3line, r4line
        elif linenum % 4 == 1: #seq line
            r1seq, r2seq, r3seq, r4seq = r1line, r2line, r3line, r4line
        elif linenum % 4 == 2: #'+' line
            r1plus, r2plus, r3plus, r4plus = r1line, r2line, r3line, r4line
        elif linenum % 4 == 3: #quality score line
            r1phred, r2phred, r3phred, r4phred = r1line, r2line, r3line, r4line
            
            indextuple = (r2seq, reverse_complement(r3seq)) # type: ignore # this tuple will be compared to keys of allpairsdict first
            seqpair = f"{indextuple[0]}_{indextuple[1]}"

            if bioinfo.qual_score(r2phred) < CUTOFF or bioinfo.qual_score(r3phred) < CUTOFF or check_n(r2seq) or check_n(r3seq):# type: ignore # first, check for unknown
                unknown = True # bad quality indexes, these are unknowns
            else:
                if indextuple in allpairsdict.keys():  #this is either hopped or matched
                    allpairsdict[indextuple] += 1 # freq +1 for this pair
                    if indextuple[0] == indextuple[1]:  
                        matched = True #matched if this tuple is in both the allpairsdict and the matchedindexdict
                    elif indextuple[0] != indextuple[1]:  
                        hopped = True # hopped if it is in allpairsdict but not in the matchedindexdict

                    else: #reads should not be ending up here
                        print(f"Error: this case is impossible. check lines: {linenum-3}-{linenum}")

                elif indextuple not in allpairsdict.keys():
                    unknown = True #this is not a hopped or a matched index (it is unknown, passes qual, and no 'N', but not an index)

                else: #reads should not be ending up here
                    print(f"this is an error, no read should be here. check lines: {linenum-3}-{linenum}")

            if unknown:
                unkcount += 1
                fname = "./output/unk" # write to files ./output/unk_R1.fq and ./output/unk_R2.fq
            elif hopped:
                hoppedcount += 1
                fname = "./output/hopped" # write to files ./output/hopped_R1.fq and ./output/hopped_R2.fq
            elif matched:
                matchedcount += 1
                fname = f"./output/{filenamingdict[r2seq]}" # write to files ./output/<index>_R1.fq and ./output/<index>_R2.fq
            else: print(f"Error: unclassified read. check lines {linenum-3}-{linenum}")

            # write to files ./output/<indexname>_R1.fq and ./output/<indexname>_R1.fq
            # here
            filedata[f"{fname}_R1.fq.gz"].write(f"{r1head} {seqpair}\n") # type: ignore
            filedata[f"{fname}_R1.fq.gz"].write(f"{r1seq}\n") # type: ignore
            filedata[f"{fname}_R1.fq.gz"].write(f"{r1plus}\n")# type: ignore
            filedata[f"{fname}_R1.fq.gz"].write(f"{r1phred}\n")# type: ignore

            filedata[f"{fname}_R2.fq.gz"].write(f"{r4head} {seqpair}\n")# type: ignore
            filedata[f"{fname}_R2.fq.gz"].write(f"{r4seq}\n")# type: ignore
            filedata[f"{fname}_R2.fq.gz"].write(f"{r4plus}\n")# type: ignore
            filedata[f"{fname}_R2.fq.gz"].write(f"{r4phred}\n")# type: ignore

        linenum += 1

# filedata["./output/A5_R1.fq"].write("Hello World") #will write "hello world" to file "./output/A5_R1.fq"

for file in filedata.values(): # close all the writing files
    file.close()


print(f"Hopped read count: {hoppedcount}")
print(f"Matched-index read count: {matchedcount}")
print(f"Unknown-indexed read count: {unkcount}")

for i in allpairsdict.keys():
    print(f"Index pair {i[0]}-{i[1]} count: {allpairsdict[i]}")