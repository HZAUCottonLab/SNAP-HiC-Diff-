import sys
#import numpy as np
from read_settings import read_settings
from create_tags import create_tags
import multiprocessing
from processor_limits import worker_count
from counts_by_chromosome_parallel_hybrids_to_5c import parallelchr
import os
import sys


'''
   Heidi Norton
   Daniel Emerson
   Harvey Huang

   Split input beds and counts in input directory by chromosome

   Uses multiprocessing which can be manually adjusted in line
                 pool = multiprocessing.Pool()
'''



def process(bed_name,counts_name,processors,tag):
        bed = open('input/' + bed_name,'r')
        output = open('input/firstfile.txt' , 'w')
        chromosomes = []
        ends = {}
        chromosome_bins = {}
        chromosome_bins_adj = {}
        adjustment = 1
        for line in bed:
                if line.startswith('#'):
                        print(line, file=output)
                else:
                        pieces = line.strip().split('\t')
                        #print 'pieces', pieces
                        chr = pieces[0]
                        start = pieces[1]
                        end = pieces[2]
                        bin  = pieces[3]
                        if chr not in chromosomes:
                                chromosome_bins[chr] = []
                                chromosome_bins_adj[chr] = []
                                output.close()
                                chromosomes.append(chr)
                                adjustment = int(bin)
                                output = open('input/' + chr + '_' + tag + '_' + bed_name, 'w')
                        mod_bin = "HiC" + str(chr).strip(" ") + "_BIN_" + str(int(bin) - adjustment)  #adjust the bin to be relative to 0
                        temp = "%s\t%s\t%s\t%s" % (chr, start, end, mod_bin)
                        print(temp, file=output)
                        chromosome_bins[chr].append(int(bin)) #unmodified list for splitting counts
                        chromosome_bins_adj[chr].append((int(bin)-adjustment) + 1) #starts at 1 instead of 0

        output.close()
        bed.close()

        region_start_stop = [] # list of tuples of chromosome start and stop bins

        for chr in chromosomes:
                ends[chr] = []
                region_start_stop.append((chr, min(chromosome_bins[chr]), max(chromosome_bins[chr]),counts_name,tag))
                ends[chr].append(max(chromosome_bins_adj[chr]))
        pool = multiprocessing.Pool(worker_count({"processors": processors}, len(region_start_stop))) #can specify number of processors if using cluster
        results = pool.map(parallelchr,region_start_stop)
        #pool.join()
        #pool.close()

        return chromosomes,ends

def main():
        settings = read_settings(sys.argv[1])
        counts_file = settings['counts_file_1']
        bed_file = settings['bed_file']
        tags = create_tags(sys.argv[1])
        tags_preprocess = tags[0]
        processors = settings["processors"]
        process(bed_file,counts_file,processors,tags_preprocess)








if __name__ == '__main__':
        main()
