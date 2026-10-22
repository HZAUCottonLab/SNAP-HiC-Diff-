import numpy as np
import os
from compute_regional_spatial_variance import compute_regional_spatial_variance
from load_region_start_end_coord import load_region_start_end_coord
from counts_to_array import counts_to_array
from load_unqueriable import load_unqueriable
from check_intersect import check_intersect
from itertools import groupby
from write_bad_subregions_chaos import write_bad_subregions
from load_region_start_end_coord import load_region_start_end_coord


def get_consecutive_zeros(a):
        b = list(range(len(a)))
        zero_intervals = []
        for group in groupby(iter(b), lambda x: a[x]):
                lis = list(group[1])
                if float(a[max(lis)]) == 0:  #only consecutive that happens to be zero
                        zero_intervals.append(max(lis) - min(lis) + 1)
                else:
                        zero_intervals.append(0) #treat as insignificant otherwise

        return max(zero_intervals)


#def heatmap_removal(counts_file,settings):
def heatmap_removal(sample, chromosomes, tag_preprocess, tags_GPS, settings):

        badregion = []
        badregion_coords = []
        region_coords = []

        mean_variances = []   ### NOT sure if this is in the right place or not -- just putting here to get it running
        super_region_idx = {}

        if 'mean_var_threshold' not in settings:
                settings['mean_var_threshold'] = 'None' #default rate
        if 'diagonal_density' not in settings:
                settings['diagonal_density'] = 0.95
        if 'consecutive_diagonal_zero' not in settings:
                settings['consecutive_diagonal_zero'] = 3
        if 'resolution' not in settings: #default resolution set to 40000
                settings['resolution'] = 40000

        if settings['badregionfilter'] == 'False':
                settings['diagonal_density'] = 0
                settings['consecutive_diagonal_zero'] = 10000000000
                settings['badregionfilter'] = 'True'

        if settings['chaosfilter'] == 'False':
                settings['chaos_pct'] = 0

        for chr in chromosomes:

                counts_file = 'input/' + chr + '_' + sample  + '_' + tag_preprocess + 'finalpvalues.counts'
                bed_file = 'input/' + counts_file.split('/')[-1][:-14] + '.bed'
                region_idx, bin_size = load_region_start_end_coord(bed_file)
                super_region_idx[chr] = region_idx

                if 'region_size' in settings:
                        desired_region_size = (int(settings['region_size']) + 1)*bin_size
                else:
                        desired_region_size = (int(settings['stride']) + int(settings['overlap']) + 1)*bin_size
                for region in region_idx:
                        region_size = region_idx[region]['end'] - region_idx[region]['start']
                        if region_size != desired_region_size:
                                regiondict = {}
                                regiondict[region] = chr
                                badregion.append(regiondict)
                                coords = super_region_idx[chr][region]
                                badregion_coords.append(coords)
                                print('lost region for being incorrect size!')
                        else:
                                print('region is correct size')



                # REMOVE SMALL REGIONS

                if 'badregionfile' in list(settings.keys()):
                        if settings['badregionfile'] != 'None':
                                bf = settings['badregionfile']
                                if not os.path.isabs(bf) and not os.path.isfile(bf):
                                        code_dir = os.environ.get('PY3_CODEDIR', os.path.dirname(os.path.abspath(__file__)))
                                        bf_candidate = os.path.join(code_dir, bf)
                                        if os.path.isfile(bf_candidate):
                                                bf = bf_candidate
                                unqueriable = load_unqueriable(bf)
                                for region in region_idx:
                                        counter = 0
                                        for element in unqueriable:
                                                if check_intersect(element, region_idx[region]):
                                                        counter +=1

                                        if counter > 0:
                                                regiondict = {}
                                                regiondict[region] = chr
                                                badregion.append(regiondict)
                                                coords = super_region_idx[chr][region]
                                                badregion_coords.append(coords)


        if settings['badregionfilter'] == 'True' or settings['chaosfilter'] == 'True':
                window = 20
                slice_start = 2
                slice_stop = slice_start + 8

                if 'chaos_pct' in list(settings.keys()):
                        pct = float(settings['chaos_pct'])
                else:
                        pct = 0.90

                metrics = {}
                metrics['maxs']={}
                metrics['variances'] = {}
                final_metrics = {}
                final_metrics['maxs']= {}
                final_metrics['variances'] = {}

                for i in range(slice_start,slice_stop+1):
                        metrics['maxs'][i] = []
                        metrics['variances'][i] = []
                metrics['maxs']['slice'] = []
                metrics['variances']['slice'] = []


                regions = []
                for chr in chromosomes:



                        counts_file = 'input/' + chr + '_' + sample  + '_' + tag_preprocess + 'finalpvalues.counts'
                        bed_file = 'input/' + counts_file.split('/')[-1][:-14] + '.bed'

                        if not os.path.isfile(counts_file) or not os.path.isfile(bed_file):
                                print("skipping region without counts or bed file:", counts_file)
                                continue
                        counts_as_arrays,bin_location = counts_to_array(counts_file, bed_file)
                        for region in counts_as_arrays:
                                region_dict = {}
                                region_dict[region] = chr
                                A = counts_as_arrays[region]
                                if settings['mean_var_threshold'] != 'None':
                                        mean_var = compute_regional_spatial_variance(counts_as_arrays[region], 3)
                                else:
                                        mean_var = 'None'
                                mean_variances.append(mean_var)
                                diag = A.diagonal()
                                iter = 0
                                slice1 = diag[:-1]
                                slice2 = diag[:-1]
                                diag = slice2
                                greater_0 = np.count_nonzero(diag)
                                pct_greater_0 = float(greater_0)/float(len(diag))
                                consecutive_zeros = get_consecutive_zeros(diag)
                                if settings['mean_var_threshold'] != 'None':
                                        if mean_var > float(settings['mean_var_threshold']):
                                                coords = super_region_idx[chr][region]
                                                badregion_coords.append(coords)
                                                badregion.append(region_dict)
                                if pct_greater_0 < float(settings['diagonal_density']):
                                #if pct_greater_0 < 0.95:
                                        coords = super_region_idx[chr][region]
                                        badregion_coords.append(coords)
                                        badregion.append(region_dict)
                                if consecutive_zeros > int(settings['consecutive_diagonal_zero']):
                                #if consecutive_zeros > 3:
                                        coords = super_region_idx[chr][region]
                                        badregion_coords.append(coords)
                                        badregion.append(region_dict)
                                coords = super_region_idx[chr][region]
                                regions.append((region_dict, coords))

                                if region not in badregion:
                                        for i in range(slice_start,slice_stop+1):
                                                slice = A.diagonal(offset=i)
                                                slice = slice[slice>0] # don't want to consider zeros when computing mean or variance
                                                if len(slice) > 0:
                                                        metrics['maxs'][i].append(max(slice))
                                                        metrics['variances'][i].append(np.var(slice, ddof=1))


                for i in range(slice_start,slice_stop+1):
                        final_metrics['maxs'][i] = np.mean(metrics['maxs'][i])
                        final_metrics['variances'][i] = np.mean(metrics['variances'][i])

                bad_idx_dict = {}

                for chr in chromosomes:
                        #print 'chr: ', chr
                        bad_idx_dict[chr] = {}
                        counts_file = 'input/' + chr + '_' + sample  + '_' + tag_preprocess + 'finalpvalues.counts'
                        bed_file = 'input/' + counts_file.split('/')[-1][:-14] + '.bed'
                        if not os.path.isfile(counts_file) or not os.path.isfile(bed_file):
                                print("skipping region without counts or bed file:", counts_file)
                                continue
                        counts_as_arrays,bin_location = counts_to_array(counts_file, bed_file)
                        region_dict = {}
                        for region in counts_as_arrays:
                                A = counts_as_arrays[region]
                                #print 'region: ', region
                                region_dict[region] = chr
                                temp_metrics = {}
                                temp_metrics['maxs'] = {}
                                temp_metrics['variances'] = {}
                                start = 0
                                stop = window
                                bad_idxs = []
                                slices = {}
                                subslices = {}
                                for i in range(slice_start,slice_stop+1):
                                        slices[i] = A.diagonal(offset = i)
                                while stop<len(slices[slice_stop]): # furthest out slice is the shortest
                                        pass_count = 0
                                        for i in range(slice_start,slice_stop+1):
                                                subslices[i] = slices[i][start:stop]
                                                subslices[i] = subslices[i][subslices[i]>0]
                                                if len(subslices[i]) > 0:
                                                        temp_metrics['maxs'][i] = max(subslices[i])
                                                        temp_metrics['variances'][i] = np.var(subslices[i], ddof=1)
                                                        if temp_metrics['maxs'][i] > final_metrics['maxs'][i]*pct:
                                                                pass_count +=1
                                                        if temp_metrics['variances'][i] > final_metrics['variances'][i]*pct:
                                                                pass_count +=1

                                        if pass_count < 1: #1*(40000/int(settings['resolution'])):  #dynamic to bin size, old_output
                                                for j in range(start,stop):
                                                        if j not in bad_idxs:
                                                                bad_idxs.append(j)
                                        start = start + 1
                                        stop = stop + 1
                                #print stop
                                bad_idxs.sort()
                                bad_idx_dict[chr][region]= bad_idxs #### NEED TO UPDATE THIS!! REGION IS NOT UNIQUE



                write_bad_subregions(bad_idx_dict, chromosomes, tag_preprocess, tags_GPS, sample)


        else:
                regions = []
                mean_variances = []
                for chr in chromosomes:
                        counts_file = 'input/' + chr + '_' + sample  + '_' + tag_preprocess + 'finalpvalues.counts'
                        bed_file = 'input/' + counts_file.split('/')[-1][:-14] + '.bed'
                        if not os.path.isfile(counts_file) or not os.path.isfile(bed_file):
                                print("skipping region without counts or bed file:", counts_file)
                                continue
                        counts_as_arrays,bin_location = counts_to_array(counts_file, bed_file)
                        for region in counts_as_arrays:
                                region_dict = {}
                                region_dict[region] = chr
                                regions.append(region_dict)
                                mean_variances.append('None')



        return regions, badregion, badregion_coords, mean_variances
