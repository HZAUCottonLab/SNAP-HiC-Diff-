import os
from outdir_helper import OUTDIR
import glob
import numpy as np
import multiprocessing
from load_good_regions_per_chrom import load_good_regions
from determine_max_gamma_per_region import determine_max_gamma_per_region
from convert_tag_genomewide import convert_tag_genomewide


def find_max_gamma(settings, tags_GPS, tags_preprocessing, samples,chrom):
        dir = OUTDIR + '/GPS/gamma_dynamic_range/'
        if not os.path.isdir(dir): os.makedirs(dir)

        max_gamma_dict_rearranged = {}
        counter = 0
        tags_GPS_test = convert_tag_genomewide(tags_GPS, settings)
        for sample in samples:

                file = dir + 'max_gamma_genomewide_' + sample + '_'+ tags_GPS + '_' + chrom + '.txt'
                file2 = dir + 'max_gamma_genomewide_' + sample + '_' + tags_GPS_test + '_' + chrom + '.txt'
                print(file)
                if os.path.isfile(file):
                        print('already found max gamma')
                        max_gamma_dict_rearranged[sample] = np.loadtxt(file)
                        counter +=1
                elif os.path.isfile(file2):
                        max_gamma_dict_rearranged[sample] = np.loadtxt(file2)
                        counter +=1
                else:
                        print('did not find max gamma file')


        print('counter', counter)
        if counter == len(samples):
                return max_gamma_dict_rearranged

        #consideration_size = 30
        consideration_size = 5
        good_regions = load_good_regions(settings, tags_GPS,chrom)
        dist_regions_of_interest = []
        if len(good_regions) > consideration_size:
                for k in range(0,int(np.floor(consideration_size))):
                        index_of_interest = int(np.floor(k*(len(good_regions)/consideration_size)))
                        dist_regions_of_interest.append(good_regions[index_of_interest])
        elif len(good_regions) == 0:
                print('no regions to consider')
                return None
        else:
                for k in range(0,len(good_regions)):
                        dist_regions_of_interest.append(good_regions[k])

        print('dist_regions_of_interest: ', dist_regions_of_interest)
        file2 = dir + 'regions_for_determining_max_gamma_'  + tags_GPS
        output = open(file2, 'w')
        for region in dist_regions_of_interest:
                out_info = region[0] + '\t' + region[1]
                print(out_info, file=output)
        output.close()



        settings_superlist = []
        chr_superlist = []
        region_superlist = []
        tags_preprocess_superlist = []
        tags_GPS_superlist = []
        cell_type_superlist = []
        for sample in samples:
                for region in dist_regions_of_interest:
                        #pieces = region.split('0')
                        #chra = pieces[0]
                        #chr = chra[3:]
                        #chr = good_regions_dict[region]
                        chr = region[1]
                        print('chr: ', chr)
                        region_superlist.append(region[0])
                        settings_superlist.append(settings)
                        chr_superlist.append(chr)
                        tags_preprocess_superlist.append(tags_preprocessing)
                        tags_GPS_superlist.append(tags_GPS)
                        cell_type_superlist.append(sample)
        max_gamma_input = list(zip(settings_superlist, chr_superlist, region_superlist, tags_preprocess_superlist, tags_GPS_superlist, cell_type_superlist))
        print('max_gamma_input: ')
        print(max_gamma_input)

        os.system("taskset -p 0xff %d" % os.getpid())
        pool = multiprocessing.Pool(int(settings['processors']))
        max_gamma_dict_list  = pool.map(determine_max_gamma_per_region, max_gamma_input)
        #for i in range(len(cell_type_superlist)):
                #max_gamma_dict = determine_max_gamma_per_region((settings_superlist[i], chr_superlist[i], region_superlist[i], tags_preprocess_superlist[i], tags_GPS_superlist[i], cell_type_superlist[i]))
        pool.close()

        print('max_gamma_dict_list: ', max_gamma_dict_list)

        max_gamma_dict_rearranged = {}
        for sample in samples:
                max_gamma_dict_rearranged[sample] = []

        for element in max_gamma_dict_list:
                sample = list(element.keys())[0]
                max_gamma_dict_rearranged[sample].append(element[list(element.keys())[0]])

        print('max_gamma_dict_rearranged: ', max_gamma_dict_rearranged)


        for sample in max_gamma_dict_rearranged:
                file = dir + 'max_gamma_genomewide_' + sample + '_'+ tags_GPS + '_' + chrom + '.txt'
                np.savetxt(file, max_gamma_dict_rearranged[sample])

        return max_gamma_dict_rearranged
