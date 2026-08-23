"""
Heidi Norton

"""
from multiprocessing import Pool
from processor_limits import worker_count
import os
import numpy as np
import sys
import glob
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import scipy
from find_consistent_communities_par import find_master_consistent_communities
from map_chr_coord import map_chr_coord
from load_community_list_DBR import load_community_list
from write_final_consistent_communities import write_final_consistent_communities
from read_settings import read_settings
from outdir_helper import OUTDIR
from create_tags import create_tags
from convert_tag_genomewide_reverse import convert_tag_genomewide_reverse

def write_consistent_communities(settings, tags_preprocessing, tags_HSVM):
        # 支持任意染色体前缀（chr* 或 HC04_A01 等）
        bed_candidates = glob.glob('input/*.*' + '_' + settings['sample_1'] + '_' + tags_preprocessing + 'final.bed')
        if not bed_candidates:
            raise FileNotFoundError("No input bed found for sample_1 and tags_preprocessing")
        bed_file = bed_candidates[0]
        chr_location = map_chr_coord(bed_file)

        print('made chr_location: ')
        chr_holder = list(chr_location.keys())[0]
        bin_size = chr_location[chr_holder][0][3] - chr_location[chr_holder][0][2]
        print('bin_size: ', bin_size)

        acceptable = ["genomewide","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",\
                                        "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",\
                                        "chr22","chr23","chr24","chr25","chr26","chr27","chr28","chr29","chr30","chr31",\
                                        "chr32","chr33","chr34","chr35","chr36","chr37","chr38","chr39","chr40","chr41",\
                                        "chr42","chr43","chr44","chr45","chr46","chr47","chr48","chr49","chr50","chrX","chrY", "chrT"]
        if settings['scale'] not in acceptable:
                  print("User did not provide acceptable format for chromosome in ",acceptable)
                  print("Must adjust to correct format in list and/or correct input bed and counts that have format")
                  return 0
        sample1_comp = []
        sample2_comp = []
        sample3_comp = []
        sample4_comp = []
        consistency1_comp = []
        consistency2_comp = []
        chromosomes = []
        d = 2*int(settings['boundary_buffer'])

        sizes = []

        for size_inst in settings:
                if size_inst[:6]=='size_s':
                         sizes.append(int(settings[size_inst]))
        sizes = sorted(sizes)

        print('sizes: ', sizes)

        available_samples = sorted([k for k in settings if k.startswith('sample_')])
        if len(available_samples) >= 4 and 'sample_3' in settings and 'sample_4' in settings:
                samples = ['sample_1', 'sample_2', 'sample_3', 'sample_4']
        elif len(available_samples) >= 2 and 'sample_2' in settings:
                samples = ['sample_1', 'sample_2']
        elif len(available_samples) >= 1:
                samples = ['sample_1']
        else:
                print('No samples found in settings')
                return None

        community_superdict = {}
        for sample in samples:
                community_superdict[sample] = {}
                if settings['chaosfilter'] == 'True':
                        dir = OUTDIR + '/HSVM/chaos_filtered_communities/'
                        fname = dir + 'Merged_Communities_*' + settings[sample] + tags_HSVM + '.txt'
                else:
                        dir = OUTDIR + '/HSVM/variance_thresholded_communities/merged/'
                        fname = dir + 'Merged_Communities_*' + settings[sample]  + '_' + tags_HSVM + '.txt'
                        print('fname: ', fname)
                if len(glob.glob(fname)) ==0:
                        tags_HSVM = convert_tag_genomewide_reverse(tags_HSVM, settings)
                        fname = dir + 'Communities_*' + settings[sample]  + '_' + tags_HSVM + '.txt'
                if settings['scale'] == 'genomewide':
                        files = glob.glob(fname)
                        print('files: ', files)
                        chromosomes = []
                        chrlist = glob.glob(fname)
                        for file in chrlist:
                                  basename = file.split('/')[-1]
                                  rest = basename.replace('Merged_Communities_', '').replace('Communities_', '')
                                  chr = rest.split(settings[sample])[0].rstrip('_')
                                  if chr not in chromosomes:
                                      chromosomes.append(chr)
                                  print('chr: ', chr)
                                  communities = load_community_list(file)
                                  print('communities: ', communities)
                                  community_superdict[sample][chr] = communities

                else:
                        print(fname)
                        if os.path.isfile(fname):
                                 print(fname)
                                 communities = load_community_list(fname)
                        else:
                                 print('fname', fname)
                                 fname = fname.replace('*', settings['scale'][3:])
                                 print('fname after second replacement: ', fname)
                                 if os.path.isfile(fname):
                                         print(fname)
                                         communities = load_community_list(fname)
                                 else:
                                         print('scale does not match that from GPS_MMCP')
                                         return None


                        community_superdict[sample][settings['scale']] = communities
                        chromosomes = []
                        chromosomes.append(settings['scale'])




        print('chromosomes: ', chromosomes)

        print('loaded communities into lists')
        pool = Pool(processes=worker_count(settings, len(community_superdict["sample_1"])))
        consistency1_comp = []
        consistency2_comp = []

        if len(samples) == 1:
                # Single sample: no consistency check; write all communities directly
                all_communities = []
                for chr in community_superdict['sample_1']:
                        all_communities += community_superdict['sample_1'][chr]
                pool.terminate()
                dir = OUTDIR + '/FINAL_DOMAIN_CALLS/'
                if not os.path.isdir(dir): os.makedirs(dir)
                write_final_consistent_communities(all_communities, dir + settings['sample_1'] + '_' + settings['scale'] + '_' + tags_HSVM + '_master_consistent_domains.txt')
        elif ('sample_3' in settings) and ('sample_4' in settings): #4 samples
                for chr in community_superdict['sample_1']:
                        if (chr in list(community_superdict['sample_2'].keys())) and (chr in list(community_superdict['sample_3'].keys())) and (list(community_superdict['sample_4'].keys())):
                                sample1_comp.append((chr,community_superdict['sample_1'][chr]))
                                sample2_comp.append((chr,community_superdict['sample_2'][chr]))
                                sample3_comp.append((chr,community_superdict['sample_3'][chr]))
                                sample4_comp.append((chr,community_superdict['sample_4'][chr]))
                                consistency1_comp.append((chr, community_superdict['sample_1'][chr], community_superdict['sample_2'][chr], d))
                                consistency2_comp.append((chr, community_superdict['sample_3'][chr], community_superdict['sample_4'][chr], d))


                results = pool.map(find_master_consistent_communities,consistency1_comp)


                results2 = pool.map(find_master_consistent_communities,consistency2_comp)

                sample1_2_master_consistent = []
                sample3_4_master_consistent = []

                for element in results:
                        sample1_2_master_consistent +=element

                for element in results2:
                        sample3_4_master_consistent +=element

                pool.terminate()


                dir = OUTDIR + '/FINAL_DOMAIN_CALLS/'
                if not os.path.isdir(dir): os.makedirs(dir)
                write_final_consistent_communities(sample1_2_master_consistent, dir + settings['sample_1'] + '_' + settings['scale'] + '_' + tags_HSVM + '_master_consistent_domains.txt')
                write_final_consistent_communities(sample3_4_master_consistent, dir + settings['sample_3'] + '_' + settings['scale'] + '_' + tags_HSVM + '_master_consistent_domains.txt')
        else:
                for chr in community_superdict['sample_1']: #only 2 samples
                        if (chr in list(community_superdict['sample_2'].keys())):
                                sample1_comp.append((chr,community_superdict['sample_1'][chr]))
                                sample2_comp.append((chr,community_superdict['sample_2'][chr]))
                                consistency1_comp.append((chr, community_superdict['sample_1'][chr], community_superdict['sample_2'][chr], d))


                results = pool.map(find_master_consistent_communities,consistency1_comp)

                sample1_2_master_consistent = []

                for element in results:
                        sample1_2_master_consistent +=element

                pool.terminate()


                dir = OUTDIR + '/FINAL_DOMAIN_CALLS/'
                if not os.path.isdir(dir): os.makedirs(dir)
                write_final_consistent_communities(sample1_2_master_consistent, dir + settings['sample_1'] + '_' + settings['scale'] + '_' + tags_HSVM + '_master_consistent_domains.txt')

def main():
        settings = read_settings(sys.argv[1])
        tags = create_tags(sys.argv[1])
        tags_HSVM = tags[3]
        tags_MMCP = tags[2]
        tags_GPS = tags[1]
        tags_preprocess = tags[0]

        write_consistent_communities(settings, tags_preprocess, tags_HSVM)











if __name__ == '__main__':
        main()
