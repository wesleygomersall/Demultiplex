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
    # parser.add_argument("-c", "--cutoff", help="Lowest allowed mean quality of indexes", type=str, required=True) 
    return parser.parse_args()

READ1 = get_args().read1
READ2 = get_args().read2
READ3 = get_args().read3
READ4 = get_args().read4
INDEXES = get_args().indexes #"../indexes.txt"
CUTOFF = 30 # get_args().cutoff ###################################CHANGE THIS ONE LATER! 

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
matchedindexdict = dict()

for i in indexsequencelist:
    for j in indexsequencelist:
        pair = (i, reverse_complement(j)) #tuple of indexes as they appear in R2 and (reverse complemented in) R3
        allpairsdict.setdefault(pair, 0) # use this dictionary to increment found (hopped and matched) pairs. 
        if i == j: 
            matchedindexdict.setdefault(pair, f"{i}_{j}") #if R2 and R3 sequences match tuple key, append value to headers. 

filenames = [] # create a list of filenames
for i in indexsequencelist: #for index-matched files, create filenames
    filenames.append(f"./output/{filenamingdict[i]}_R1.fq")
    filenames.append(f"./output/{filenamingdict[i]}_R2.fq")
for i in ["unk", "hopped"]: #need to make files for unknown reads and for hopped reads
    filenames.append(f"./output/{i}_R1.fq")
    filenames.append(f"./output/{i}_R2.fq")

filedata = {filename: open(filename, 'w') for filename in filenames} #open all the files to write to

with open(READ1, 'r') as r1, open(READ2, 'r') as r2, open(READ3, 'r') as r3, open(READ4, 'r') as r4: #open to read all 4 files
    linenum = 0
    while True: 
        unknown: bool = False
        matched: bool = False
        hopped: bool = False

        r1line, r2line, r3line, r4line = r1.readline().strip('\n'), r2.readline().strip('\n'), r3.readline().strip('\n'), r4.readline().strip('\n')

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

            # whenever we reach the end of an entry: 
            if bioinfo.qual_score(r2phred) < CUTOFF or bioinfo.qual_score(r3phred) < CUTOFF: # first, check for quality of r2phred and r3phred
                unknown = True # bad quality indexes, these are unknowns
            else:
                indextuple = (r2seq, r3seq)

            # next check if the sequences r2seq and r3seq are in keys of allpairsdict
                # not inside -> unknown = True, 
                # is inside -> unknown = False (definitively), must be hopped or matched, incriment here
                    #if (is inside), check if it is in matchedindexdict 
                        # is in matchedindexdict
                            # matched index = True
                            # add value of matchedindexdict to headers: r1head and r4head
            
            if unknown:
                print("unk read")
                # need to append to headers (reverse comp R3 seq)
                # write to files ./output/unk_R1.fq and ./output/unk_R2.fq
            if hopped:
                print("hopped read")
                # need to append to headers (reverse comp R3 seq)
                # write to files ./output/hopped_R1.fq and ./output/hopped_R2.fq
            if matched:
                print("matched read")
                # need to append to headers (made dictionary for this one)
                # write to files ./output/<indexname>_R1.fq and ./output/<indexname>_R1.fq
        linenum += 1

# filedata["./output/A5_R1.fq"].write("Hello World") #will write "hello world" to file "./output/A5_R1.fq"

for file in filedata.values(): # close all the writing files
    file.close()



