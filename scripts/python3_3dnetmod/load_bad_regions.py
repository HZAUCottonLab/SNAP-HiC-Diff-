from convert_tag_genomewide import convert_tag_genomewide
from convert_tag_genomewide_reverse import convert_tag_genomewide_reverse
from outdir_helper import OUTDIR
import os
import glob


def load_bad_regions(settings,tag,chrom):

    tag_test = convert_tag_genomewide(tag, settings)
    tag_test2 = convert_tag_genomewide_reverse(tag, settings)
    print('tag_test: ', tag_test)
    samples= []
    for sample_inst in settings:
         if sample_inst[:7] == 'sample_':
             samples.append(settings[sample_inst])
    samples = sorted(samples)
    all_sample_tag = ''
    for sample in samples:
         all_sample_tag = all_sample_tag + sample + '_'

    bad_region_file = OUTDIR + '/GPS/bad_region_removal/bad_regions_' + all_sample_tag + tag + '_'+ chrom + '.txt'
    bad_region_file2 = OUTDIR + '/GPS/bad_region_removal/bad_regions_' + all_sample_tag + tag_test + '_' + chrom +  '.txt'
    bad_region_file3 = OUTDIR + '/GPS/bad_region_removal/bad_regions_' + all_sample_tag + tag_test2 + '_' + chrom +  '.txt'
    print('bad_region_file: ', bad_region_file)
    print('bad_region_file2: ', bad_region_file2)
    if os.path.isfile(bad_region_file):
        print('bad region file: ', bad_region_file)
        badregions = []
        input = open(bad_region_file, 'r')
        for line in input:
            if line.startswith('#'):
                continue
            else:
                reg = line.strip().split('\t')[0]
                chrom = line.strip().split('\t')[1]
                #badregions.append(line.strip().split('\t')[0])
                badregions.append((reg, chrom))
        print('loaded bad regions: ', badregions)
    elif os.path.isfile(bad_region_file2):
        print('bad region file: ', bad_region_file2)
        badregions = []
        input = open(bad_region_file2, 'r')
        for line in input:
            if line.startswith('#'):
                continue
            else:
                reg = line.strip().split('\t')[0]
                chrom = line.strip().split('\t')[1]
                #badregions.append(line.strip().split('\t')[0])
                badregions.append((reg, chrom))
    elif len(glob.glob(bad_region_file3)) > 0:
        bad_region_file = glob.glob(bad_region_file3)[0]
        print('bad_region_file', bad_region_file)
        badregions = []
        input = open(bad_region_file, 'r')
        for line in input:
            if line.startswith('#'):
                continue
            else:
                reg = line.strip().split('\t')[0]
                chrom = line.strip().split('\t')[1]
                #badregions.append(line.strip().split('\t')[0])
                badregions.append((reg, chrom))



        print('loaded bad regions: ', badregions)
    else:
        print('did not find bad region file: ')
        return None
    return badregions
