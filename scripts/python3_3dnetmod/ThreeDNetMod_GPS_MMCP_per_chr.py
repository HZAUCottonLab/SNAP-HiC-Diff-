import sys
from read_settings import read_settings
#from load_bad_regions_simplified_per_chrom import load_bad_regions
from load_bad_regions import load_bad_regions
from gamma_plateau_sweep import GPS
from modularity_maximization_consensus_partition import MMCP

'''
   Heidi Norton
   Daniel Emerson
   Harvey Huang

   Determines gamma range per chr region from settings file parameters using GPS (gamma plateau sweep).
   Identifies communities using modularity maximization and consensus partition (MMCP).
'''

def collective(element):
     chr = element[0]
     settings = element[1]
     sample = element[2]
     tags_GPS = element[3]
     tags_MMCP = element[4]
     tags_preprocess = element[5]
     max_gamma = element[6]
     chrom = element[7]
     badregions = load_bad_regions(settings,tags_GPS,chrom)
     test = GPS(chr,settings,sample,badregions,tags_GPS,tags_MMCP,tags_preprocess, max_gamma)
     if test == 0:
         MMCP(chr,settings,sample,badregions,tags_GPS,tags_MMCP,tags_preprocess)
     print('completed MMCP')

     return 0

def main():
        settings = read_settings(sys.argv[1])
        chr = 'chrX'
        sample = settings['sample_1']
        element = [chr, settings, sample]
        collective(element)

if __name__ == '__main__':

        main()
