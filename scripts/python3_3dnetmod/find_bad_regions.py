from heatmap_removal import heatmap_removal
from outdir_helper import OUTDIR
import numpy as np
import sys
import glob
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def find_bad_regions(settings, samples,tag_preprocess,tag_GPS, chromosomes,chrom):
     print('chromosomes: ', chromosomes)
     total_badregions = []
     total_badregions_coords = []
     total_regions = []
     for sample in samples:
          total_finalvar = []
          regions, badregions, badregion_coords, finalvar = heatmap_removal(sample, chromosomes, tag_preprocess, tag_GPS, settings)
          print('regions: ')
          print(regions)

          total_badregions = total_badregions + badregions
          total_regions = total_regions + regions
          total_badregions_coords = total_badregions_coords + badregion_coords
     masterlist_badregions = total_badregions
     # Python 3: dicts are not comparable; sort by (region_name, chrom)
     def _badregion_key(d):
         k, v = list(d.keys())[0], list(d.values())[0]
         return (k, v)
     masterlist_badregions = sorted(masterlist_badregions, key=_badregion_key)
     masterlist_regions = total_regions
     # Python 3: elements may be (dict, coords); sort by comparable key
     def _region_key(r):
         if isinstance(r, tuple):
             d = r[0]
         else:
             d = r
         return (list(d.keys())[0], list(d.values())[0])
     masterlist_regions = sorted(masterlist_regions, key=_region_key)
     total_badregions_coords = sorted(
         total_badregions_coords,
         key=lambda c: (c.get('chrom', c.get('chr', '')), c.get('start', 0), c.get('end', 0))
     )

     #print "masterlist_badregions =", masterlist_badregions

     unique_bad_regions = []
     for region in masterlist_badregions:
          if region not in unique_bad_regions:
              unique_bad_regions.append(region)

     unique_regions = []
     unique_region_coords = []
     for region in masterlist_regions:
          if region[0] not in unique_regions:
              unique_regions.append(region[0])
              unique_region_coords.append(region[1])

     print('number of unique regions: ', len(unique_regions))
     print('number of unique_region-coords: ', len(unique_region_coords))

     unique_good_regions = []
     unique_good_region_coords = []
     for i in range(len(unique_regions)):
         if unique_regions[i] not in unique_bad_regions:
              unique_good_regions.append(unique_regions[i])
              unique_good_region_coords.append(unique_region_coords[i])

     all_sample_tag = ''
     for sample in samples:
         all_sample_tag = all_sample_tag + sample + '_'

     dir = OUTDIR + '/GPS/bad_region_removal/'
     if not os.path.isdir(dir): os.makedirs(dir)
     fname = 'good_regions_' + all_sample_tag + tag_GPS + '_' + chrom + '.txt'
     outfile = open(dir + fname, 'w')
     for region in unique_good_regions:
         output = list(region.keys())[0] + '\t' + region[list(region.keys())[0]]
         print(output, file=outfile)
     outfile.close()

     fname = 'bad_regions_' + all_sample_tag + tag_GPS + '_' + chrom + '.txt'
     outfile = open(dir + fname, 'w')
     for region in unique_bad_regions:
          output = list(region.keys())[0] + '\t' + region[list(region.keys())[0]]
          print(output, file=outfile)
     outfile.close()

     fname = 'bad_region_genomic_coords_' + all_sample_tag + tag_GPS + '_' + chrom + '.txt'
     outfile = open(dir + fname, 'w')
     for coords in total_badregions_coords:
          output = coords['chrom'] + '\t' + str(coords['start']) + '\t' + str(coords['end'])
          print(output, file=outfile)
     outfile.close()

     fname = 'good_region_genomic_coords_' + all_sample_tag + tag_GPS + '_' + chrom + '.txt'
     outfile = open(dir + fname, 'w')
     for coords in unique_good_region_coords:
          output = coords['chrom'] + '\t' + str(coords['start']) + '\t' + str(coords['end'])
          print(output, file=outfile)
     outfile.close()



     return unique_bad_regions
