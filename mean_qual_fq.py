#!/usr/bin/env python

import bioinfo
import matplotlib.pyplot as plt

def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". List length will be 101. If no value is passed, initializes list
    with 101 values of 0.0.'''
    for i in range(101):
        lst.append(value)
    return lst

my_list: list = []
my_list = init_list(my_list)

def populate_list(file: str) -> tuple[list, int]:
    """Takes a fastq file path as string input. 
    Does some counting and outputs the sum of quality scores at each postition 
    in the sequence, for read lengths of exactly 101 bases.
    Also outputs the total line number. 
    """
    list = []
    list = init_list(list)
    with open(file) as fqall:
        for index, line in enumerate(fqall):
            if (index+1) % 4 == 0:
                # print(line.strip('\n'))
                # print(len(line.strip('\n')))
                for col, pchar in enumerate(line.strip('\n')):
                    list[col] += bioinfo.convert_phred(pchar)
    return list, (index+1)

file1 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
file2 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'

my_list, num_lines = populate_list(file1)
my_list = [n/(num_lines/4) for n in my_list] # turns sum to average

print("# Base Pair", "Mean Quality Score", sep="\t")
for i in range(len(my_list)):
    print(f"{i}\t{my_list[i]}") # make table of mean quality scores for each position.

position = range(len(my_list))

plt.figure(figsize=(10,6))
plt.bar(position, height = my_list) # no need to have std dev on plots 

plt.title('Average Phred Quality Scores for Read 1')
plt.xlabel('Position in sequence')
plt.ylabel('Converted quality score')