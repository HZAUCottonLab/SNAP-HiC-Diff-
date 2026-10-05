from read_settings import read_settings
from outdir_helper import OUTDIR
from load_community_dict import load_community_dict
from load_region_start_end_coord import load_region_start_end_coord
from create_tags import create_tags
from hierarchical_spatial_variance_minimization import HSVM
import sys
import glob
import os
from remove_chaos_communities import remove_chaos_communities
from write_consistent_communities import write_consistent_communities
import numpy as np

## script
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)

def main():

     settings = read_settings(sys.argv[1])
     tags = create_tags(sys.argv[1])
     tags_HSVM = tags[3]
     tags_MMCP = tags[2]
     tags_GPS = tags[1]
     tags_preprocess = tags[0]


     dir = OUTDIR + '/FINAL_DOMAIN_CALLS/'
     if not os.path.isdir(dir): os.makedirs(dir)

     dir = OUTDIR + '/HSVM/size_hist/'
     if not os.path.isdir(dir): os.makedirs(dir)

     samples= []
     for sample_inst in settings:
         if sample_inst[:7] == 'sample_':
             samples.append(settings[sample_inst])
     samples = sorted(samples)
     for sample in samples:
         pre_threshold_file = 'Communities_*' + sample + '*_' + tags_MMCP + '.txt'  #test to make sure unique communities are present

         dir = OUTDIR + '/MMCP/unique_communities/results_files/'
         print('intermediate files: ', dir + pre_threshold_file)
         if  not os.path.isdir(dir):
             print("No communities directory present")
             return 0
         existence = glob.glob(dir + pre_threshold_file)

         if len(existence) == 0:
             pre_threshold_file = 'Communities_*' + sample + '*_' + tags_MMCP + '.txt'
             pre_threshold_file = pre_threshold_file.replace('_' + 'genomewide', '_' + 'chr*')
             existence = glob.glob(dir + pre_threshold_file)
             print(existence)
             if len(existence) == 0:
                 print("No intermediate file present. exiting")
                 return 0

         master_community_dict = {}
         existence.sort() #sort subchrom regions

         last_length = 0
         prevent_inf_loop = 0
         while last_length < 3:  #ensure last community list is not empty
             last_length = len(np.genfromtxt(existence[-1]))
             backend = existence.pop()
             existence.insert(0,backend)
             if prevent_inf_loop > len(existence):
                 break
             prevent_inf_loop = prevent_inf_loop + 1

         print("length of existence = ", len(existence))

         for existing_file in existence:
             print('existing_file: ', existing_file)
             _basename = existing_file.split('/')[-1]
             _rest = _basename.replace('Communities_', '', 1)
             cell_type = _rest.split(sample)[0] + sample
             print('cell_type: ', cell_type)
             chr = cell_type.split('.')[0]
             print('chr: ', chr)
             sub_chrom = cell_type.split(sample)[0]
             print('sub_chrom: ', sub_chrom)
             if chr not in master_community_dict:
                    master_community_dict[chr] = {}
             cell_type_input = cell_type.replace(sample,'') + "_" + sample
             #chr18.5_CKKO1Neuronalchr18
             print('cell_type_input: ', cell_type_input)
             community_dict = load_community_dict(existing_file)
             if len(community_dict) == 0:
                 print("file empty")
                 #del(master_community_dict[chr][sub_chrom])
             else:
                 master_community_dict[chr][sub_chrom] = community_dict
                 print('community_dict: ', community_dict)
         print("")
         print("their is community_dict:")
         print(master_community_dict[chr])
         print("")

         for chr in master_community_dict:
             dir2 = OUTDIR + '/HSVM/variance_thresholded_communities/results_files/merged/'
             final_output = 'Merged_' + 'Communities_' + chr + cell_type  + tags_HSVM + '.txt'
             existence = glob.glob(dir2 + final_output)
             #print "going to HSVM variance threshold"
             #print "chr: ", chr
             #for keys in master_community_dict[chr]:
             #     print keys
             if len(existence) > 0: #already found output based on merged
                 continue
             else:
                 #variance_threshold = HSVM(settings,master_community_dict[chr],cell_type_input,cell_type,tags_HSVM,tags_preprocess)
                 variance_threshold = HSVM(settings,master_community_dict[chr],sample, tags_HSVM,tags_preprocess, chr)

     if settings['chaosfilter'] == 'True':
         remove_chaos_communities(settings, tags_GPS, tags_HSVM)

     if 'final_consistent_domains' in list(settings.keys()):
         if settings['final_consistent_domains'] == 'True':
             write_consistent_communities(settings, tags_preprocess, tags_HSVM)


if __name__ == "__main__":
      main()
