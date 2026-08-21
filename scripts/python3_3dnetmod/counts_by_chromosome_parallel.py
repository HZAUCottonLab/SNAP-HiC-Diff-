import numpy as np
from read_settings import read_settings

def parallelchr(element):

     #settings = read_settings('settings_genomewide.txt')
     #counts = open('input/' + settings["counts_file_1"], 'r')
     tag = element[0][4]
     counts = open('input/' + str(element[0][3]), 'r')
     #counts_out = open('input/' + element[0] + '_' + str(settings["counts_file_1"]), 'w')
     counts_out = open('input/' + element[0][0] + '_' + tag + '_' + str(element[0][3]), 'w')
     lookat = 0
     for line in counts:
          pieces = line.strip().split("\t")
          bin_1 = int(pieces[0])
          bin_2 = int(pieces[1])
          value = pieces[2]
          while lookat < len(element) and bin_1 > element[lookat][2]:
               counts_out.close()
               lookat = lookat + 1
               if lookat < len(element):
                    counts_out = open("input/" + element[lookat][0] + "_" + tag + "_" + str(element[lookat][3]), "w")
          if lookat >= len(element):
               break
          if bin_1 < element[lookat][1]:
               continue
          modbin_1 = (bin_1 - element[lookat][1]) + 1
          modbin_2 = (bin_2 - element[lookat][1]) + 1
          temp = "%s\t%s\t%f" % (modbin_1, modbin_2, float(value))
          print(temp, file=counts_out)
     counts_out.close()
     return 0


if __name__ == "__main__":
     parallelchr()
