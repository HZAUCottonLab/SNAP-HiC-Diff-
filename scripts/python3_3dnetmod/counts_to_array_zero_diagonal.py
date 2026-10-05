
"""
05/21/16
Heidi Norton and Jennifer Phillips-Cremins
This module loads binned 5C .counts data into a dictionary where key = region and value = square symmetric matrix of counts.
Input:
binned, normalized .counts file of 5C data
bed file containing genomic coordinates of 5C data
num_diag: number of diagonals that should be zeroed. If num_diag = 1, the single diagonal is zeroed out.
Output:
counts_as_arrays: a dictionary where key = 5C region name, value = square symmetric matrix of counts. Elements along the diagonal zeroed.
bin_location: a nested dictionary with keys as regions and bin numbers within regions. Values are a list of genomic coordinates for each bin number.

The main script loads 5C counts files and writes square symmetric .csv matrices of the counts for each region.

"""

import sys
import numpy as np
import os


def counts_to_array_zero_diagonal(counts_filename, bed_filename, num_diag):
        def _region_name(bin_name):
                if '_BIN_' in bin_name:
                        return bin_name.rsplit('_BIN_', 1)[0]
                return bin_name.split('_')[0]

        def _bin_index(bin_name):
                parts = bin_name.split('_')
                if len(parts) > 3 and parts[2] == 'BIN':
                        return int(parts[3])
                return int(parts[2])


        # Create a list of 5C regions by parsing input .counts file
        regions = []
        input = open(counts_filename,'r')
        for line in input:
                if line[0:1] == '#':
                        continue
                else:
                        region = _region_name(line.strip().split('\t')[0])
                        if region not in regions:
                                regions.append(region)

        input.close()


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

        # Create a dictionary of number of bins per region
        num_bins = {}
        for region in regions:
                num_bins[region] = len(bin_location[region])


        # Initialize data structure
        counts_as_arrays = {}
        for region in regions:
                counts_as_arrays[region] = np.zeros((num_bins[region], num_bins[region]), dtype=np.float64)

        input = open(counts_filename, 'r')
        for line in input:
                if line[0:1] == '#':
                        continue
                else:
                        bin_f = line.strip().split('\t')[0]
                        bin_r = line.strip().split('\t')[1]
                        norm_count = float(line.strip().split('\t')[2])
                        region = _region_name(bin_f)
                        index_bin_f = _bin_index(bin_f)
                        index_bin_r = _bin_index(bin_r)
                        if abs(index_bin_f - index_bin_r) < num_diag:
                                counts_as_arrays[region][index_bin_f,index_bin_r] = 0
                                counts_as_arrays[region][index_bin_r,index_bin_f] = 0
                        else:
                                counts_as_arrays[region][index_bin_f,index_bin_r] = norm_count
                                counts_as_arrays[region][index_bin_r,index_bin_f] = norm_count

        input.close()

        return counts_as_arrays, bin_location





def main():
        counts_filename, bed_filename = sys.argv[1:]
        counts_as_arrays, bin_location = counts_to_array('input/' + counts_filename, 'input/' + bed_filename)
        dir = 'input/counts_as_matrix/'
        if not os.path.isdir(dir): os.makedirs(dir)
        for region in counts_as_arrays:
                np.savetxt(dir + 'Matrix_' + region +'_' + counts_filename + '.csv', counts_as_arrays[region], delimiter=',', fmt='%10.5f')


if __name__ == '__main__':
        main()
