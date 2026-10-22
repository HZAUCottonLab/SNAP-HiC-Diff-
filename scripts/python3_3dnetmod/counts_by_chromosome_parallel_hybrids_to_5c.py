import numpy as np
from read_settings import read_settings
import sys

def parallelchr(element):

     #settings = read_settings('settings_genomewide.txt')
     tag = element[4]
     counts = open('input/' + element[3], 'r')
     counts_out = open('input/' + element[0] + '_' + tag + '_' + element[3], 'w')
     for line in counts:
          pieces = line.strip().split('\t')
          bin_1 = pieces[0]
          bin_2 = pieces[1]
          value = pieces[2]
          if int(bin_1) >= element[1] and int(bin_1) <= element[2]:
               modbin_1 = "HiC"+ str(element[0]).strip(" ") + "_BIN_" + str(int(bin_1) - element[1]) #adjust bin 1 to start at 1 per region
               modbin_2 = "HiC" + str(element[0]).strip(" ") + "_BIN_" + str(int(bin_2) - element[1]) #adjust bin 2 relative to bin 1
               temp = "%s\t%s\t%f" % (modbin_1, modbin_2, float(value))
               print(temp, file=counts_out)
          if int(line.strip().split('\t')[0]) > element[2]:
               break
     return 0


if __name__ == "__main__":
     parallelchr()
