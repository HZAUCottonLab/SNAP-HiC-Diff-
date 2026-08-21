import os
import sys
import re
import random
import copy
import csv
import numpy as np


def _bin_index(bin_name):
     parts = bin_name.split('_')
     if len(parts) > 3 and parts[2] == 'BIN':
          return int(parts[3])
     return int(parts[2])

'''
     split counts into uniform regions with overlap and no deadzone removal.  Assumes 5C format
     for input counts and bed (region_BIN_###)

'''


class MatrixMake:
      Maximum = 0
      Minimum = 100000
      logged = False

      def __init__(self, filename,filename2,start,stop,first,last,number,region,logged,stride,overlap,qnorm,find_char_counts,final_start):
           self.filename = filename  #input counts
           self.filename2 = filename2
           self.start = int(start)
           self.stop = int(stop)
           self.first = int(first)
           self.last = int(last)
           self.pcount = []
           self.list_of_bins = []
           self.list_of_regions = []
           self.number = number
           MatrixMake.logged = logged
           self.stop_processing = False
           self.stride = stride
           self.overlap = overlap
           self.end = 0
           self.qnorm = qnorm
           self.final_start = final_start
           self.find_char_counts = find_char_counts
           self.current_char_counts = find_char_counts
           self.increase_current_char_counts = True

           for i in range(int(self.start),int(self.stop)+1):
               self.list_of_bins.append(i)

           MatrixMake.Maximum = int(self.last)
           MatrixMake.Minimum = int(self.first)

           self.Heatmap = np.zeros(((self.stop-self.start)+1,(self.stop-self.start)+1))
           self.region = region
           self.readfile()
           self.fillemptydiags()
           self.pCountMake()


      def readfile(self):

           f = open(self.filename,"r+")
           f.seek(self.find_char_counts)
           counter = 0
           testing_start = True
           for lines in f.readlines():
                if testing_start:
                    print(lines)
                    testing_start = False
                offset = len(lines) #must account for \n for character offset
                counter = counter + 1
                lines = lines.strip('\n')
                values = lines.split('\t')
                values = [i.strip(' ') for i in values]
                self.test(values)
                if self.increase_current_char_counts:
                     self.current_char_counts = self.current_char_counts  +  offset
                if (lines == ""): #EOF
                     self.end = 1
                if self.stop_processing:  #keep from wasting too much time processing lines if exceeds the end
                     counter = -10
                     break
           f.close()

           if counter >= 0:  #end
                self.end = 1

           f = open(self.filename2,"r+")
           processing = 1
           for lines in f.readlines():
               lines = lines.strip('\n')
               values = lines.split('\t')
               values = [i.strip(' ') for i in values]
               self.test2(values,processing)
           f.close()

      def test(self,values):
           #test to see if input line from counts is suitable for storage

           if eval(self.qnorm):
                firstcol = _bin_index(values[0])
                secondcol = _bin_index(values[1])
           else:
                firstcol = values[0]
                secondcol = values[1]


           if ((int(firstcol) >= int(self.start)) and (int(firstcol) <= int(self.stop))):
                if ((int(secondcol) >= int(self.start)) and (int(secondcol) <= int(self.stop))):
                       bin1 = int(firstcol) - self.start
                       bin2 = int(secondcol) - self.start
                       if ((int(firstcol) >= (self.stop - self.overlap)) or (int(firstcol) >= self.final_start)):
                            self.increase_current_char_counts = False
                       if eval(MatrixMake.logged):
                            if float(values[2]) >= 1:
                                 self.Heatmap[bin1][bin2] = np.log2(float(values[2]))
                                 self.Heatmap[bin2][bin1] = np.log2(float(values[2]))
                            else:
                                 self.Heatmap[bin1][bin2] = 0
                                 self.Heatmap[bin2][bin1] = 0
                       else:
                            if float(values[2]) >= 0: #remove neg value from consideration
                                 self.Heatmap[bin1][bin2] = float(values[2])
                                 self.Heatmap[bin2][bin1] = float(values[2])
                            else:
                                 self.Heatmap[bin1][bin2] = 0
                                 self.Heatmap[bin2][bin1] = 0

           if int(firstcol) > self.stop:
                self.stop_processing = True

      def test2(self,values,processing):

           if eval(self.qnorm):
                forthcol = _bin_index(values[3])
           else:
                forthcol = values[3]
           binstart = int(forthcol) - self.start
           digits = len(str(MatrixMake.Maximum))

           if (int(forthcol) in self.list_of_bins):
               self.list_of_regions.append([values[0],values[1],values[2],str(self.region) + str(self.number).zfill(5) + "_BIN_" + str(binstart).zfill(digits)])  #attach all rows of .bed input to matching bins

      def fillemptydiags(self):
          #fill in empty diagonals
          diag_of_interest = []
          #find empty diagonals
          for i in range(0,len(self.Heatmap)):
              #diagonal is all zeros
              if float(len(self.Heatmap.diagonal(i)[self.Heatmap.diagonal(i)==0]))/float(len(self.Heatmap.diagonal(i))) == 1:
                  diag_of_interest.append(i)
              else:
                  break #immediately exit, avoid looking at off diagonals far from diagonal
          if len(diag_of_interest) > 0:
               np.fill_diagonal(self.Heatmap,0.0001) #if center is all zeros reset to 0.0001
               for col in range(len(diag_of_interest),1,-1):  #start at furthest from center diagonal
                   col_of_interest = diag_of_interest[col-1]
                   col_gap = diag_of_interest[col-1]
                   #only evaluate if there is a neighbor to right (diagonals not equal in size)
                   while col_of_interest < (len(self.Heatmap) - 1):
                       value = self.Heatmap[col_of_interest-(1*col_gap)][col_of_interest + 1]
                       self.Heatmap[col_of_interest-(1*col_gap)][col_of_interest] = value
                       self.Heatmap[col_of_interest][col_of_interest-(1*col_gap)] = value
                       col_of_interest = col_of_interest + 1
               #adjust center diagonal
               col_of_interest = 0
               while col_of_interest < (len(self.Heatmap) - 1):
                   value = self.Heatmap[col_of_interest][col_of_interest+1]
                   self.Heatmap[col_of_interest][col_of_interest] = value
                   col_of_interest = col_of_interest + 1
               for col in range(len(diag_of_interest)):
                   self.Heatmap[len(self.Heatmap) - 1 - col][len(self.Heatmap) - 1 - col] = 0.0001 #prevent left over zeros from influencing bad region removal



      def pCountMake(self):
          #store pcount data

          digits = len(str(MatrixMake.Maximum))
          temppcount =  []
          b = self.start + 1
          for i in range(0,len(self.Heatmap)):
               for a in range(self.start,b):
                   self.pcount.append([str(self.region) + str(self.number).zfill(5) + "_BIN_" + str(i).zfill(digits),str(self.region) + str(self.number).zfill(5) + "_BIN_" + str(a-self.start).zfill(digits),str(self.Heatmap[i][a-self.start]).strip(" ")])
               b = b + 1


      def save(self,fileName,fileName2,fileName3):

               np.savetxt("input/" + fileName, self.Heatmap, delimiter=',',fmt='%f')

               f = open("input/" + fileName2,"w") # pCount
               f.write('\n'.join('\t'.join('{:1}'.format(item) for item in row) for row in self.pcount))
               f.write('\n')
               f.close();

               f = open("input/" + fileName3,"w") # Bed
               f.write('\n'.join('\t'.join('{:1}'.format(item) for item in row) for row in self.list_of_regions))
               f.write('\n')
               f.close();



def process4(element):
     #looks for qnormed counts

     chr = element[0]
     bed = element[1]
     counts = element[2]
     sample = element[3]
     overlap = int(element[4])
     stride = int(element[5])
     logged = element[6]
     size = int(element[7])
     qnorm = element[8]
     tag = element[9]
     copies = int(np.floor((size-(overlap +1))/stride))
     extension = (size - (overlap+1)) % stride

     if chr == 'chrM':
         return (chr,1)

     name='HiC'+ str(chr) # region name shared across replicates, appended to bin inside output counts and bed file
     file=str(chr) + "_" + tag + "_" + str(sample)# Update 'NAME' to be incorporated into output file name
     if eval(qnorm):
         counts_file='input/' + str(sample) + str(chr) + '_' + tag + '_qnormed.counts' #Update 'COUNTS' to the file path of specific chromosome counts file
     else:
         counts_file='input/' + str(chr) + '_' + tag + "_" + str(counts)
     bed_file='input/' + str(chr) + "_" + tag + "_" + bed #Update 'BED' to the .bed file path of specific bed counts file

     if eval(qnorm):
          begin = 0  #5C format
          start = 0 #5C format
     else:
          begin = 1 #HiC format (starts at 1)
          start = 1  #HiC format

     last = size
     end = True
     stop = start + stride
     rep = 1

     if overlap > stride:
          if eval(qnorm):
               begin = overlap #5C format
               start = overlap #5C format
          else:
               begin = 1 + overlap #HiC format
               start = 1 + overlap #HiC format
          stop = start + stride
          copies = int(np.floor((size-(overlap + 1))/stride))

     #number_maj_regions = 30
     maj_region_size = 3
     if copies >= maj_region_size:
          remainder = copies % maj_region_size
          #maj_region_size = int((copies - remainder)/number_maj_regions)
          number_maj_regions = int((copies - remainder)/maj_region_size)
     else:
          maj_region_size = 1
          number_maj_regions = copies
          remainder = 0



     if eval(qnorm):
          #final_start = size-((stride+overlap)+1)
          #final_end = size-1 #final bin value for 5C (starts at 0)
          final_start = (stride*copies)
          final_end = size - 1
     else:
          #final_start = size - (stride+overlap)
          final_start = stride*copies + 1
          final_end = size #final bin value for HiC (starts at 1)

     majreg = 1

     seek = 0
     rep = 1
     for i in range(0,number_maj_regions):
          for j in range(0,maj_region_size):
               file=str(chr) + '.' + str(majreg) + "_" + tag + "_" + str(sample)
               data = MatrixMake(counts_file,bed_file,start,stop,begin,last,rep,name,logged,stride,overlap,qnorm,seek,final_start)
               data.save(file + str(rep).zfill(5) + "_heatmap.csv",file + str(rep).zfill(5) + "_pvalues.counts",file + str(rep).zfill(5) +".bed")
               start = stop - overlap
               stop = stop + stride
               rep = rep+1
               seek = data.current_char_counts
               print(seek)
          majreg = majreg + 1

     majreg = number_maj_regions + 1

     for j in range(0,remainder):  #last super region is remainder plus left over after copy #
          file=str(chr) + '.' + str(majreg) + "_" + tag + "_" + str(sample)
          data = MatrixMake(counts_file,bed_file,start,stop,begin,last,rep,name,logged,stride,overlap,qnorm,seek,final_start)
          data.save(file + str(rep).zfill(5) + "_heatmap.csv",file + str(rep).zfill(5) + "_pvalues.counts",file + str(rep).zfill(5) +".bed")
          start = stop - overlap
          stop = stop + stride
          rep = rep+1
          seek = data.current_char_counts
          print(seek)

     if extension > 0:  #remainder after dividing into copies
          file=str(chr) + '.' + str(majreg) + "_" + tag + "_" + str(sample)
          data = MatrixMake(counts_file,bed_file,final_start,final_end,begin,last,rep,name,logged,stride,overlap,qnorm,seek,final_start)
          data.save(file + str(rep).zfill(5) + "_heatmap.csv",file + str(rep).zfill(5) + "_pvalues.counts",file + str(rep).zfill(5) +".bed")

     #number of super regions equals copy number plus (remainder + leftover)
     return (chr,number_maj_regions+1)  #return number of super regions per chromosome

if __name__ == "__main__":
     process4()
