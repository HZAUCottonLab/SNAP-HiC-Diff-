from outdir_helper import OUTDIR


def load_good_regions(settings,tag,chrom):
    samples= []
    good_regions_dict = {}
    for sample_inst in settings:
         if sample_inst[:7] == 'sample_':
             samples.append(settings[sample_inst])
    samples = sorted(samples)
    all_sample_tag = ''
    for sample in samples:
         all_sample_tag = all_sample_tag + sample + '_'

    good_region_file = OUTDIR + '/GPS/bad_region_removal/good_regions_' + all_sample_tag + tag + '_' + chrom +'.txt'
    goodregions = []
    input = open(good_region_file, 'r')
    for line in input:
        if line.startswith('#'):
             continue
        else:
             region = line.strip().split('\t')[0]
             chrom = line.strip().split('\t')[1]
             if (region[-3:] !='001'):
                  goodregions.append((region, chrom))

    print('loaded good regions except for the first region of each chromosome ', goodregions)
    return goodregions
