#!/usr/bin/env python

import bioinfo
import matplotlib.pyplot as plt
import gzip

def init_list(lst: list, n: int, value: float=0.0) -> list:
    '''This function takes an empty list and populates it with
    the value passed in "value" for a total length of n.
    '''
    for i in range(n):
        lst.append(value)
    return lst

def populate_list(file: str, n:int) -> tuple[list, int]:
    """Takes a fastq file path as string input. Also needs to know n, the length of untrimmed sequences.
    Does some counting and outputs the sum of quality scores at each postition 
    in the sequence, for read lengths of exactly 101 bases.
    Also outputs the total line number. 
    """
    list = []
    list = init_list(list, n)
    with gzip.open(file, 'rt') as fqall:
        for index, line in enumerate(fqall):
            if (index+1) % 4 == 0:
                # print(line.strip('\n'))
                # print(len(line.strip('\n')))
                for col, pchar in enumerate(line.strip('\n')):
                    list[col] += bioinfo.convert_phred(pchar)
    return list, (index+1)

file1 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
file2 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz'
file3 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz'
file4 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz'

my_list1, num_lines1 = populate_list(file1, 101)
my_list1 = [n/(num_lines1/4) for n in my_list1] # turns sum to average

my_list2, num_lines2 = populate_list(file2, 8)
my_list2 = [n/(num_lines2/4) for n in my_list2] # turns sum to average

my_list3, num_lines3 = populate_list(file3, 8)
my_list3 = [n/(num_lines3/4) for n in my_list3] # turns sum to average

my_list4, num_lines4 = populate_list(file4, 101)
my_list4 = [n/(num_lines4/4) for n in my_list4] # turns sum to average

plt.figure(figsize=(10,6))
plt.bar(range(len(my_list1)), height = my_list1) # no need to have std dev on plots 
plt.title('Average Phred Quality Scores for Read 1')
plt.xlabel('Position in sequence')
plt.ylabel('Converted quality score')
plt.savefig(f"./meanqual_Read1.png")

plt.figure(figsize=(10,6))
plt.bar(range(len(my_list2)), height = my_list2) # no need to have std dev on plots 
plt.title('Average Phred Quality Scores for Index 1')
plt.xlabel('Position in sequence')
plt.ylabel('Converted quality score')
plt.savefig(f"./meanqual_Index1.png")

plt.figure(figsize=(10,6))
plt.bar(range(len(my_list3)), height = my_list3) # no need to have std dev on plots 
plt.title('Average Phred Quality Scores for Index 2')
plt.xlabel('Position in sequence')
plt.ylabel('Converted quality score')
plt.savefig(f"./meanqual_Index2.png")

plt.figure(figsize=(10,6))
plt.bar(range(len(my_list4)), height = my_list4) # no need to have std dev on plots 
plt.title('Average Phred Quality Scores for Read 2')
plt.xlabel('Position in sequence')
plt.ylabel('Converted quality score')
plt.savefig(f"./meanqual_Read2.png")
