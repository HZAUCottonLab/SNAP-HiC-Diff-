from create_bin_location_dict import create_bin_location_dict
from operator import itemgetter
from itertools import groupby
import os
from outdir_helper import OUTDIR

#def write_bad_subregions(bad_subregion_dict, bed_file, tags_GPS, sample):
def write_bad_subregions(bad_subregion_dict, chromosomes, tag_preprocess, tags_GPS, sample):
        for chr in chromosomes:
                counts_file = 'input/' + chr + '_' + sample  + '_' + tag_preprocess + 'finalpvalues.counts'
                bed_file = 'input/' + counts_file.split('/')[-1][:-14] + '.bed'

                bin_location = create_bin_location_dict(bed_file)

                dir = OUTDIR + '/GPS/chaos_regions/'
                if  not os.path.isdir(dir): os.makedirs(dir)
                fname = dir +  chr + sample + tags_GPS + '.txt'
                output = open(fname, 'w')
                header = '#chr\tstart\tstop'
                print(header, file=output)

                for region in bad_subregion_dict[chr]:
                        print('regions in bin_location_dict: ')
                        print(list(bin_location.keys()))
                        print('chr: ', chr)
                        print('region: ', region)
                        print('bins in bin_location: ')
                        print(list(bin_location[region].keys()))
                        a = bad_subregion_dict[chr][region]
                        for k, g in groupby(enumerate(a), lambda i_x:i_x[0]-i_x[1]):
                                group = list(map(itemgetter(1), g))
                                print('group: ', group)
                                start_bin = min(group)
                                end_bin = max(group)
                                #chr = bin_location[region][start_bin][1]
                                start_coord = bin_location[region][start_bin][2]
                                end_coord = bin_location[region][end_bin][3]

                                temp = chr + '\t' + str(start_coord) + '\t' + str(end_coord)
                                print(temp, file=output)


                output.close()
