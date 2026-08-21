from find_bad_regions import find_bad_regions
from load_bad_regions import load_bad_regions
from find_max_gamma import find_max_gamma
from counts_to_array import counts_to_array
from ThreeDNetMod_GPS_MMCP_per_chr import collective
from create_tags import create_tags
import multiprocessing
from read_settings import read_settings
from outdir_helper import OUTDIR
import glob
import os
import sys

# 运行目录（含 input/）；集群上 cwd 可能不对，用 PY3_RUNDIR 保证绝对路径
RUN_DIR = os.environ.get("PY3_RUNDIR", os.getcwd())

## script
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)

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
     tags = create_tags(sys.argv[1])
     tags_GPS = tags[1]
     tags_MMCP = tags[2]
     tags_preprocess = tags[0]


     #initialize all possible directories
     dir = OUTDIR + '/GPS/bad_region_removal'
     if not os.path.isdir(dir): os.makedirs(dir)

     dir = OUTDIR + '/GPS/chaos_regions'
     if not os.path.isdir(dir): os.makedirs(dir)

     dir = OUTDIR + '/GPS/gamma_dynamic_range'
     if not os.path.isdir(dir): os.makedirs(dir)

     dir = OUTDIR + '/MMCP/communities_across_numpart'
     if not os.path.isdir(dir): os.makedirs(dir)

     dir = OUTDIR + '/MMCP/variance_barplots'
     if not os.path.isdir(dir): os.makedirs(dir)

     dir = OUTDIR + '/MMCP/unique_communities/results_files'
     if not os.path.isdir(dir): os.makedirs(dir)

     dir = OUTDIR + '/GPS/random_networks'
     if not os.path.isdir(dir): os.makedirs(dir)

     dir = OUTDIR + '/GPS/scatterplots'
     if not os.path.isdir(dir): os.makedirs(dir)

     dir = OUTDIR + '/MMCP/consensus_filtered_files'
     if not os.path.isdir(dir): os.makedirs(dir)

     dir = OUTDIR + '/MMCP/communities_pre_filters'
     if not os.path.isdir(dir): os.makedirs(dir)

     pre_chromosomes = []
     chromosomes = {}

     key_name = settings['sample_1']
     # 支持任意染色体前缀；用 RUN_DIR 保证集群上也能找到 input/
     input_dir = os.path.join(RUN_DIR, 'input')
     pattern = os.path.join(input_dir, '*.*' + str(key_name) + '*' + tags_preprocess + '*final.bed')
     values = glob.glob(pattern)


     settings['scale'] = settings['scale'].strip("'")
     settings['scale'] = settings['scale'].strip('"')

     samples= []
     chrom_o = []

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
               pre = lines.strip().split(os.sep)
               fname = pre[-1]  # 支持绝对路径时取文件名，如 HC04_A01.100_H1398_150_100_Truefinal.bed
               # 染色体.区域 在 _sample_ 前：xxx.H1398_... -> chr=xxx
               parts = fname.split('_' + str(key_name) + '_')
               chr = parts[0] if len(parts) >= 2 else fname.split('_')[0]
               if '.' not in chr:  #not a superregion.  Skip to next element in list
                   continue
               chr_main = chr.split('.')[0]  #actual main chr
               if chr not in pre_chromosomes: #only store one occurence of super region name from all samples
                    if chr != 'allchr':
                        pre_chromosomes.append(chr)
               if chr_main not in chrom_o:
                    chrom_o.append(chr_main)
     else:
          chromosomes[settings['scale']] = []
          chrom_o.append(settings['scale'])
          findchr = glob.glob(os.path.join(input_dir, settings['scale'] + '.*' + str(key_name) + '_' + tags_preprocess + '*final.bed'))
          print('findchr: ', findchr)
          for element in findchr:
               chrom = element.split('/')[-1].split('_')[0]
               if '.' not in chrom:  #not a superregion.  Skip to next element in list
                    continue
               chromosomes[settings['scale']].append(chrom)

     print("chromosomes = ", chrom_o)
     print("pre_chromosomes = ",pre_chromosomes)

     all_sample_tag = ''
     for sample in samples:
         all_sample_tag = all_sample_tag + sample + '_'

     if settings['scale'] == 'genomewide':
         for i in range(0,len(pre_chromosomes)):
             #print pre_chromosomes[i]
             superregion = pre_chromosomes[i]
             found_chr_main = superregion.split('.')[0]
             if found_chr_main not in chromosomes:
                 chromosomes[found_chr_main] = []
             chromosomes[found_chr_main].append(superregion)

     total_bad_regions = []
     print('chromosomes: ', chromosomes)

     for i in range(0,len(chrom_o)):
         if not os.path.isfile(OUTDIR + '/GPS/bad_region_removal/bad_regions_' + all_sample_tag + tags_GPS +  '_' + str(chrom_o[i]) +  '.txt'):
             bad_regions = find_bad_regions(settings, samples, tags_preprocess, tags_GPS, chromosomes[chrom_o[i]],chrom_o[i])
             total_bad_regions.append(bad_regions)
         else:
             print('found bad region file')
             print('file: ', OUTDIR + '/GPS/bad_region_removal/bad_regions_' + all_sample_tag + tags_GPS +  '_' + str(chrom_o[i]) +  '.txt')
             bad_regions = load_bad_regions(settings,tags_GPS,chrom_o[i])
             total_bad_regions.append(bad_regions)


     total_max_gamma_return = {}
     skipped_chroms = set()
     for i in range(0,len(chrom_o)):
         max_gamma_return = find_max_gamma(settings, tags_GPS, tags_preprocess, samples,chrom_o[i])
         if max_gamma_return is None:
             print('WARNING: no good regions for chromosome %s, skipping' % chrom_o[i])
             skipped_chroms.add(chrom_o[i])
         total_max_gamma_return[chrom_o[i]] = max_gamma_return

     input = []
     for i in range(len(chrom_o)): #iterate through all chromosomes
         if chrom_o[i] in skipped_chroms:
             continue
         for j in range(len(samples)): #iterate through all celltypes
             for k in range(len(chromosomes[chrom_o[i]])): #iterate through all chromosome super regions
                 sample = samples[j]
                 print('sample: ', sample)
                 if 'max_gamma' in list(settings.keys()):
                      max_gamma = float(settings['max_gamma'])
                 elif 'single_gamma' in list(settings.keys()):
                      max_gamma = float(settings['single_gamma'])
                 else:
                      max_gamma = max(total_max_gamma_return[chrom_o[i]][sample]) #max gamma per cell type per chr (across all super regions)
                 print('max_gamma: ', max_gamma)
                 #collective((chromosomes[j], settings,samples[i],tags_GPS,tags_MMCP,tags_preprocess, max_gamma))
                 input.append((chromosomes[chrom_o[i]][k], settings,samples[j],tags_GPS,tags_MMCP,tags_preprocess, max_gamma,chrom_o[i]))


     os.system("taskset -p 0xff %d" % os.getpid())
     pool = multiprocessing.Pool(int(settings['processors']))

     print('entering MMCP')


     results = pool.map(collective,input) #run per chromosome


def wrapper(args):
     return collective(*args)


if __name__ == "__main__":
     main()
