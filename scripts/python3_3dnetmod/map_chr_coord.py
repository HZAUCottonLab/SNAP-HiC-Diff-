"""
This module creates a dictionary of 5C bin locations.
Input: .bed file of the format
#Chromosome     Start Site      End Site        Bin ID
Ouptut: nested dictionary of bin_locations:
bin_location[region][bin_index] = [bin_name, chr, start, end]
"""

def map_chr_coord(bed_filename):
        # Create a bin location dictionary
        input = open(bed_filename, 'r')
        bin_location = {} # nested dictionary with keys as regions and bin numbers within regions and values as a list of genomic coordinates.
        bin_names = []
        bin_number = 0
        previous_chr = 0
        for line in input:
                if line[0:1] == '#':
                        continue
                else:
                        chr = line.strip().split('\t')[0]
                        start = int(line.strip().split('\t')[1])
                        end = int(line.strip().split('\t')[2])
                        if chr not in bin_location:
                                bin_location[chr] = {}
                        bin_name = chr + '_' + str(start) + '_' + str(end)
                        if bin_name not in bin_names:
                                bin_names.append(bin_name)
                                if previous_chr == chr:
                                        bin_location[chr][bin_number] = [bin_name, chr, start, end]
                                        bin_number +=1
                                else:
                                        bin_number = 0
                                        bin_location[chr][bin_number] = [bin_name, chr, start, end]
                                        bin_number +=1
                        previous_chr = chr
        input.close()

        return bin_location
