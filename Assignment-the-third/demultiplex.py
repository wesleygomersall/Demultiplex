#!/usr/bin/env python

import argparse
import gzip
import bioinfo



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


# open all files to write to
    # open all files to read from 
        # get 4 lines from each file: rin1, rin2, rin3, rin4
        # 
