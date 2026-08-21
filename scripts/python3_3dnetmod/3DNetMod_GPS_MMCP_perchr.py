from find_bad_regions_mean_var_simplified_per_chrom import find_bad_regions
from load_bad_regions_simplified_per_chrom import load_bad_regions
#from determine_max_gamma_genomewide_simplified import determine_max_gamma
from find_max_gamma_simplified_per_chrom import find_max_gamma
from counts_to_array import counts_to_array
from ThreeDNetMod_GPS_MMCP_per_chr_simplified_per_chrom import collective
from create_tags import create_tags
import multiprocessing
from read_settings import read_settings
from outdir_helper import OUTDIR
import glob
import os
import sys

'''
   Heidi Norton
   Daniel Emerson
   Harvey Huang

   3DNetMod_GPS_MMCP.py determines modularity and outputs communities.  Uses output
   from prelim_genomewide.py.  Uses user input settings file to retrieve input settings.
   Final community calls post threshold are placed in directory:
        output/unique_communities_thresholded_post/results_files/

   If running locally keep line: pool = multiprocessing.Pool()
   This will automatically set to number of cores available for multiprocessing

   If on cluster, multiprocessing.Pool() can be adjusted to number of threads desired
        pool = multiprocessing.Pool(#_of_threads)

   3DNetMod_GPS_MMCP.py retrieves from user input settings file the desired chr or
   genomewide from key settings['scale'] and runs 3DNetMod_GPS_MMCP_per_chr.py per chr

'''



def main():

     settings = read_settings(sys.argv[1])
     chrom_o = sys.argv[2]
     tags = create_tags(sys.argv[1])
     tags_GPS = tags[1]
     tags_MMCP = tags[2]
     tags_preprocess = tags[0]

     chromosomes = []

     key_name = settings['sample_1']
     values = glob.glob('input/' + chrom_o + '.*' + str(key_name) + '*' + tags_preprocess + '*final.bed')



     settings['scale'] = settings['scale'].strip("'")
     settings['scale'] = settings['scale'].strip('"')

     samples= []

     for sample_inst in settings:
         if sample_inst[:7] == 'sample_':
             samples.append(settings[sample_inst])

     samples = sorted(samples)


     #test against list of acceptable values for 'scale'
     acceptable = ["genomewide","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",\
                    "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",\
                    "chr22","chr23","chr24","chr25","chr26","chr27","chr28","chr29","chr30","chr31",\
                    "chr32","chr33","chr34","chr35","chr36","chr37","chr38","chr39","chr40","chr41",\
                    "chr42","chr43","chr44","chr45","chr46","chr47","chr48","chr49","chr50","chrX","chrY"]

     if settings['scale'] not in acceptable:
          print("User did not provide acceptable format for chromosome in ",acceptable)
          print("Must adjust to correct format in list and/or correct input bed and counts that have format")
          return 0
     if settings['scale'] == 'genomewide':
          print('scale is genomewide')
          for lines in values:
               pre = lines.strip().split('/')
               ext = pre[1].split('_')
               chr = ext[0]  #chr in prefix separated by "_"
               if chr not in chromosomes:
                    if chr != 'allchr':
                        chromosomes.append(chr)
     else:
          #print "found ", settings['scale']
          #chromosomes.append(settings['scale'])
          findchr = glob.glob('input/' + settings['scale'] + '*' + str(key_name) + '_' + tags_preprocess + '*final.bed')
          for element in findchr:
               chrom = element.split('/')[-1].split('_')[0]
               chromosomes.append(chrom)
     print("chromosomes = ", chromosomes)

     all_sample_tag = ''
     for sample in samples:
         all_sample_tag = all_sample_tag + sample + '_'

     if not os.path.isfile(OUTDIR + '/GPS/bad_region_removal/bad_regions_' + all_sample_tag + tags_GPS + '_' + chrom_o + '.txt'):
         #bad_regions = find_bad_regions(settings, samples,tags_preprocess,tags_GPS)
         bad_regions = find_bad_regions(settings, samples, tags_preprocess, tags_GPS, chromosomes, chrom_o)
     else:
         print('found bad region file')
         bad_regions = load_bad_regions(settings,tags_GPS,chrom_o)

     print('chromosomes: ')
     print(chromosomes)

     max_gamma_return = find_max_gamma(settings, tags_GPS, tags_preprocess, samples, chrom_o)


     input = []
     for i in range(len(samples)): #iterate through all celltypes
         for j in range(len(chromosomes)): #iterate through all chromosomes
             sample = samples[i]
             print('sample: ', sample)
             max_gamma = max(max_gamma_return[sample])
             print('max_gamma: ', max_gamma)
             #collective((chromosomes[j], settings,samples[i],tags_GPS,tags_MMCP,tags_preprocess, max_gamma))
             input.append((chromosomes[j], settings,samples[i],tags_GPS,tags_MMCP,tags_preprocess, max_gamma,chrom_o))

     os.system("taskset -p 0xff %d" % os.getpid())
     pool = multiprocessing.Pool(int(settings['processors']))

     print('entering MMCP')


     results = pool.map(collective,input) #run per chromosome

def wrapper(args):
     return collective(*args)


if __name__ == "__main__":
     main()
