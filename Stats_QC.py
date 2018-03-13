__author__ = 'Eleonor Cayford'

import sys
import agate
import subprocess

# Defines a function used to get statistics on the fastq file passed into it
# Prints the number of reads, number of bases, average sequence length, and GC content
def get_stats(file_name):
    # Opens the file in a with block to ensure that it's correctly closed afterwards
    with open(file_name) as data:
        filename = file_name.split('.')[0]
        with open(filename+'_stats.csv', 'w') as stats_file:
            print 'Collecting stats on ' + file_name
            # Used to keep count of the number of bases in the complete file
            base_counter = 0
            # Used to keep count of the number of reads in the complete file
            reads_counter = 0
            # Used to keep count of the number of G and C bases in the complete file
            gc_counter = 0
            # For loop runs through the file 4 lines at a time, starting with the sequence identifier beginning with "@"
            # This allows for each sequence and its information to be handled together
            for seq_ID in data:
                # The next 3 lines of code collect the sequence, its sequence ID beginning with "+" , and its quality
                sequence = next(data).strip()
                seq_ID_plus = next(data).strip()
                quality = next(data).strip()
                # Adds the bases in the current sequence to the base counter
                base_counter = base_counter + len(sequence.strip())
                # Adds the count of Gs and Cs in the current sequence to the GC counter
                gc_counter = gc_counter + sequence.count('G') + sequence.count('C')
                # Increases the read count by one since each group of 4 lines contains one read
                reads_counter += 1
            # For naming and clarity, the number of reads is set equal to the total read counter
            num_reads = reads_counter
            # The numBases is the average number of bases per sequence in this file,
            # found with the total number of bases divided by the total number of reads
            num_bases = base_counter/num_reads
            # Finds the decimal GC content by dividing the count of GC bases by the total number of bases
            gc_content = float(gc_counter)/base_counter
            stats_file.write('Number_of_Reads,Number_of_Bases,Average_Read_Length,GC_Content,'+ '\n')
            stats_file.write(str(num_reads) + ',' + str(base_counter) + ',' + str(num_bases) + ',' + str(gc_content) + ',' + '\n')
            print 'Stat export complete on ' + file_name

# Creates a function to find and export the quality and base content information on passed in files
def export_quality_bases(file_name):
    # Creates an empty list for storing each possible base and their count
    base_counts = []
    # Makes base_counts a list of dictionaries with the value 0 assigned to each of the possible base bases
    for x in range(100):
        base_counts.append({base: 0 for base in 'ATCGN'})
    # Opens the processed data file that we want
    with open(file_name) as data:
        filename = file_name.split('.')[0]
        with open(filename+'_quality_scores.csv', 'w') as average_quality_score_file:
            print 'Collecting quality scores on ' + file_name
            #Creates a string of numbers from 1 to 100
            average_quality_score_file.write(','.join(map(str, range(1,101))) + '\n')
            for seq_ID in data:
                sequence = next(data).strip()
                seq_ID_plus  = next(data).strip()
                quality = next(data).strip()
                # Increments the count for the current base in the dictionary mapping to the current read location
                for read_location, base in enumerate(sequence):
                    base_counts[read_location][base] += 1
                # Translates the ASCII character squence into a series of phred quality scores
                quality_scores = [str(q-33) for q in map(ord, quality)]
                # Writes out the quality scores to a CSV file
                average_quality_score_file.write(','.join(quality_scores) + '\n')
            print 'Quality score export complete on ' + file_name
    # Opens a csv which has the per read location base contents written to it
    with open(filename+'_base_content.csv', 'w') as base_file:
        print 'Collecting base content information on ' + file_name
        base_file.write('A,T,C,G,N\n')
        for row in base_counts:
            base_file.write(','.join([str(row[base]/float(sum(row.values()))) for base in 'ATCGN'])+'\n')
        print 'Base content export complete on ' + file_name
"""---------------------------------------------------------------------------------------------------------------"""

subprocess.call(['java', '-jar', 'trimmomatic-0.35.jar', 'SE', '-phred33', 'raw_data.fastq', 'processed_data.fastq', 'ILLUMINACLIP:TruSeq3-SE:2:30:10', 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:80'])
get_stats('raw_data.fastq')
get_stats('processed_data.fastq')
export_quality_bases('raw_data.fastq')
export_quality_bases('processed_data.fastq')