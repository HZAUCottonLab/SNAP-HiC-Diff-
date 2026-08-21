"""
This module creates a dictionary of 5C bin locations.
Input: .bed file of the format
#Chromosome     Start Site      End Site        Bin ID
Ouptut: nested dictionary of bin_locations:
bin_location[region][bin_index] = [bin_name, chr, start, end]
"""

def create_bin_location_dict(bed_filename):
        def _region_name(bin_name):
                if '_BIN_' in bin_name:
                        return bin_name.rsplit('_BIN_', 1)[0]
                return bin_name.split('_')[0]

        def _bin_index(bin_name):
                parts = bin_name.split('_')
                if len(parts) > 3 and parts[2] == 'BIN':
                        return int(parts[3])
                return int(parts[2])

        # Create a bin location dictionary
        input = open(bed_filename, 'r')
        bin_location = {} # nested dictionary with keys as regions and bin numbers within regions and values as a list of genomic coordinates.
        for line in input:
                if line[0:1] == '#':
                        continue
                else:
                        chr = line.strip().split('\t')[0]
                        start = int(line.strip().split('\t')[1])
                        end = int(line.strip().split('\t')[2])
                        bin_name = line.strip().split('\t')[3]
                        region = _region_name(bin_name)
                        index_bin = _bin_index(bin_name)
                        if region not in bin_location:
                                bin_location[region] = {}
                        bin_location[region][index_bin] = [bin_name, chr, start, end]
        input.close()

        return bin_location
