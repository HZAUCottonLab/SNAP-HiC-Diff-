import sys
#import numpy as np
from read_settings import read_settings
import multiprocessing
from counts_by_chromosome_parallel import parallelchr
import os



'''
   Heidi Norton
   Daniel Emerson
   Harvey Huang

   Split input beds and counts in input directory by chromosome

   Uses multiprocessing which can be manually adjusted in line
                 pool = multiprocessing.Pool()
'''



def process_HiC(element):
        #settings = read_settings("settings_genomewide.txt")
        #bed = open('input/' + settings["bed_file"],'r')
        bed = open('input/' + str(element[0]),'r')
        counts = element[1]
        output = open('input/firstfile.txt' , 'w')
        chromosomes = []
        ends = []
        chromosome_bins = {}
        chromosome_bins_adj = {}
        adjustment = 1
        tag = element[3]
        processors = element[2]
        for line in bed:
                if line.startswith('#'):
                        print(line, file=output)
                else:
                        pieces = line.strip().split('\t')
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
                                #output = open('input/' + chr + '_' + str(settings["bed_file"]), 'w')
                                output = open('input/' + chr + '_' + tag + '_' + str(element[0]), 'w')
                        mod_bin = (int(bin) - adjustment) + 1  #adjust the bin to be relative to 1
                        temp = "%s\t%s\t%s\t%s" % (chr, start, end, mod_bin)
                        print(temp, file=output)
                        chromosome_bins[chr].append(int(bin))
                        chromosome_bins_adj[chr].append((int(bin)-adjustment) + 1)

        output.close()
        bed.close()

        region_start_stop = [] # list of tuples of chromosome start and stop bins

        for chr in chromosomes:
                region_start_stop.append((chr, min(chromosome_bins[chr]), max(chromosome_bins[chr]),counts,tag))
                ends.append(max(chromosome_bins_adj[chr]))
        #pool = multiprocessing.Pool(processors) #can specify number of processors if using cluster
        #results = pool.map(parallelchr,region_start_stop)
        parallelchr(region_start_stop)

        return chromosomes,ends







if __name__ == '__main__':
        process()
