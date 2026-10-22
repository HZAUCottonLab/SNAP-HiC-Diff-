from read_settings import read_settings
from outdir_helper import OUTDIR
from create_tags import create_tags
import sys
import os
import glob
from load_intervals import load_intervals
from load_intervals import load_stripped_intervals
from check_intersect import check_intersect
from get_overlap import get_overlap

def remove_chaos_communities(settings, tag_GPS, tag_HSVM):
        print('in module')
        final_communities = []
        buff = 4 # number of bins buffer for removing communities due to chaos


        samples = []
        for sample_inst in settings:
                if sample_inst[:7] == 'sample_':
                        samples.append(settings[sample_inst])
        samples = sorted(samples)
        for sample in samples:
                files = glob.glob(OUTDIR + '/HSVM/variance_thresholded_communities/merged/Merged_' + 'Communities_*' + sample  + '_' + tag_HSVM + '.txt')
                print('HSVM files:')
                print(files)
                print('looking for: ', OUTDIR + '/HSVM/variance_thresholded_communities/merged/Merged_' + 'Communities_*' + sample  + '_' + tag_HSVM + '.txt')
                # Merged_Communities_chr18WT2Neuronalchr18_100_50_True_False_mm10_cent_tel_genedesert_False_0.01_chr18_1_0.0048_None_3_1_20_0_0_4_s1_400kb_s2_800kb_s3_1600kb_s4_3000kb_s5_12000kb_v1_90_v2_100_v3_100_v4_100_v5_100_40000.txt


                # chaos regions
                #output/GPS/chaos_regions/chr18.23WT2Neuronalchr18100_50_True_False_mm10_cent_tel_genedesert_False_True_0.01_chr18_1_0.0048_None_3_1.txt
                for file in files:
                        final_communities = []
                        print('file: ')
                        print(file)
                        basename = file.split('/')[-1]
                        rest = basename.replace('Merged_Communities_', '').replace('Communities_', '')
                        chr = rest.split(sample)[0].rstrip('_')
                        community_list = load_intervals(file)
                        print('chr: ', chr)
                        chaos_files = glob.glob(OUTDIR + '/GPS/chaos_regions/' + chr + '*' + sample +  tag_GPS + '.txt')
                        if len(chaos_files) < 1:
                                if settings['scale'] == 'genomewide':
                                        pre_threshold_file = OUTDIR + '/GPS/chaos_regions/' + chr + '*' + sample +  tag_GPS + '.txt'
                                        pre_threshold_file = pre_threshold_file.replace('_' + 'genomewide', '_' + str(chr))
                                        chaos_files = glob.glob(pre_threshold_file)
                        #chaos_files = glob.glob(OUTDIR + '/GPS/chaos_regions/' + chr + '*' + sample +  tag_HSVM + '.txt')
                        print('looking for chaos file: ')
                        print(OUTDIR + '/GPS/chaos_regions/' + chr + '*' + sample +  tag_GPS + '.txt')
                        print('chaos_files: ', chaos_files)
                        chaos_regions = []
                        for chaos_file in chaos_files:
                                print('chaos file in for loop')
                                chrom = chaos_file.split('/')[-1].split(sample)[0]
                                if chrom == chr:
                                        chaos = load_stripped_intervals(chaos_file)
                                        chaos_regions += chaos
                        print('chaos_regions: ', chaos_regions)
                        for community in community_list:
                                counter = 0
                                for chaos_region in chaos_regions:
                                        overlap = get_overlap(community,chaos_region)
                                        #print 'buff: ', buff
                                        #print "int(settings['boundary_buffer'])", int(settings['boundary_buffer'])
                                        #print "buff*(int(settings['boundary_buffer']))", buff*(int(settings['boundary_buffer']))
                                        #print 'community: ', community
                                        #print 'chaos_region: ', chaos_region
                                        #print 'overlap: ', overlap

                                        if (abs(chaos_region['start'] - community['start']) <= buff*(int(settings['boundary_buffer']))) and (abs(chaos_region['end'] - community['end']) <= buff*(int(settings['boundary_buffer']))): # keep community if it has similar coordinates as the chaos region (it's likely an outer boundary to the chaos region)
                                                continue
                                        elif (community['start'] < chaos_region['start']) and (community['end'] > chaos_region['end']): # if chaos region is inside the community
                                                continue
                                        elif overlap >= buff*(int(settings['boundary_buffer'])): # if community overlaps with chaos region significantly, we want to remove it
                                                counter +=1
                                                break
                                        elif (community['end'] - community['start'] < buff*(int(settings['boundary_buffer']))) and check_intersect(community,chaos_region):
                                                print('lost community to final condition')
                                                print("(community['end'] - community['start'] < buff*(int(settings['boundary_buffer']))) and check_intersect(community,chaos_region)")
                                                print(community)
                                                print(chaos_region)
                                                counter +=1
                                                break
                                        else:
                                                continue
                                if counter == 0:
                                        final_communities.append(community)
                                        print('appended final community')
                                else:
                                        print('lost community to chaos')


                        dir = OUTDIR + '/HSVM/chaos_filtered_communities/'
                        if not os.path.isdir(dir): os.makedirs(dir)
                        out_file = open(dir + 'Merged_Communities_' + chr + sample + tag_HSVM + '.txt', 'w')
                        for item in final_communities:
                                output = item['chrom'] + '\t' + str(item['start']) + '\t' + str(item['end'])
                                print(output, file=out_file)
                        out_file.close()


def main():

        settings = read_settings(sys.argv[1])
        tags = create_tags(sys.argv[1])
        tags_HSVM = tags[3]
        tag_GPS = tags[1]
        remove_chaos_communities(settings, tag_GPS, tags_HSVM)


if __name__ == '__main__':
        main()
