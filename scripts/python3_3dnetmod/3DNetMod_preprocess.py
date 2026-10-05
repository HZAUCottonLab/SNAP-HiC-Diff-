from split_counts_by_chromosome_hybrids_to_5c import process
from split_counts_by_chromosome import process_HiC
#from demo_remove_primer_primer_pairs_genomewide import process2
#from demo_quantile_normalization_genomewide import process3
from slice_genomic_coordinates_genomewide_majregv import process4
from create_tags import create_tags
from processor_limits import worker_count
from read_settings import read_settings
import glob
import sys
import multiprocessing
import subprocess
import os
import numpy as np

## script
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)

def configured_excluded_chromosomes(settings):
    raw = settings.get("excluded_chromosomes", "")
    return {value.strip() for value in raw.split(",") if value.strip()}

def prelim_genomewide(settings,tag,scale):

    bed = settings['bed_file']
    overlap = settings['overlap']

    if 'region_size' in settings:
        stride = int(settings['region_size']) - int(overlap)
    else:
        stride = int(settings['stride'])
    logged = settings['logged']
    excluded_chromosomes = configured_excluded_chromosomes(settings)

    counts = []
    counts_names = []
    for counts_inst in settings:
        if counts_inst[:12] == 'counts_file_':
            counts_names.append(counts_inst)

    counts_names.sort()  #make sure properly sorted before assigning to counts list
    for counts_name_inst in counts_names:
            counts.append(settings[counts_name_inst])  #properly sorted so counts_file_1,counts_file_2,..

    samples= []
    samples_names = []
    increment = 0
    for sample_inst in settings:
        if scale == 'genomewide':
            if sample_inst[:7] == 'sample_':
                samples_names.append(sample_inst)
                test_files = glob.glob('input/*' + settings[sample_inst] + '_' + tag + 'finalpvalues.counts')
                if len(test_files) > 0:  #found input do not calculate
                    increment = increment + 1
        else:
            if sample_inst[:7] == 'sample_':
                samples_names.append(sample_inst)
                test_files = glob.glob('input/' + scale + '.*' + settings[sample_inst] + '_' + tag + 'finalpvalues.counts')
                if len(test_files) > 0:  #found input do not calculate
                    increment = increment + 1

    if increment == len(samples_names):
        print("already found stitched counts")
        return 0

    samples_names.sort()  #make sure properly sorted before assigning to sample names list
    for sample_name_inst in samples_names:
        samples.append(settings[sample_name_inst]) #properly sorted so sample_1,sample_2,sample_3,..

    #find all chromosomes and end bin values based on one of celltype counts and bed
    collection_prelim = []
    for i in range(len(counts)):
         collection_prelim.append([bed,counts[i],int(settings['processors']),tag])
    pool = multiprocessing.Pool(worker_count(settings, len(collection_prelim)))
    results_1 = pool.map(process_HiC, collection_prelim)
    #for i in range(len(counts)): #iterate through all celltypes with corresponding counts
    #     results_1 = process_HiC([bed,counts[i],int(settings['processors']),tag])  #phase 1

    print("results_1 = ",results_1)
    chr = results_1[0][0]
    ends = results_1[0][1]

    for mito in excluded_chromosomes:
        if mito in chr:
            index_remove = chr.index(mito)
            del chr[index_remove]
            del ends[index_remove]
            command_removal = 'rm input/' + mito + '_' + tag + '_' + settings["bed_file"]
            os.system(command_removal)
            for counts_inst in counts:
                test_files = glob.glob('input/' + mito + '_' + tag + '_' + counts_inst)
                if len(test_files) > 0:
                    command_removal = 'rm input/' + mito + '_' + tag + '_' + counts_inst
                    os.system(command_removal)

    settings_repeat = []
    for i in range(0,len(chr)):
        settings_repeat.append(settings)

    print("chr = ",chr)
    print("end = ",ends)

    complete = []
    for i in range(len(samples)): #iterate through all celltypes with corresponding counts
        if settings['scale'] == 'genomewide':
            for j in range(0,len(chr)):
                 complete.append((chr[j],bed,counts[i],samples[i],overlap,stride,logged,ends[j],'False',tag))
        else:
            for j in range(0,len(chr)):
                 if chr[j] == settings['scale']:
                     complete.append((chr[j],bed,counts[i],samples[i],overlap,stride,logged,ends[j],'False',tag))
    pool = multiprocessing.Pool(worker_count(settings, len(complete)))
    results = pool.map(process4, complete)
    #for i in range(0,len(complete)):
    #    process4(complete[i])
    pool.close()
    pool.join()

    #phase 3


    #chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']

    unique_results = set(results)


    output_beds = {}
    for sample_inst in samples:
        for chromosome in unique_results:
            for i in range(1,chromosome[1] + 1):
                print("removing chromosome meta files")
                if chromosome[0] not in excluded_chromosomes:
                    if str(chromosome[0]) not in list(output_beds.keys()):
                        output_beds[str(chromosome[0])] = []

                input_list = 'input/' + str(chromosome[0]) + '.' + str(i) + '_' + tag + '_' + sample_inst + '0*.bed' #assume < 100 regions per chromosome
                output_name = 'input/' + str(chromosome[0]) + '.' + str(i) + '_' + sample_inst  + '_' + tag + 'final.bed'

                command = 'cat '+ input_list + ' > ' + output_name
                if len(glob.glob('input/' + str(chromosome[0]) + '.' + str(i) + '_' + tag + '_' + sample_inst + '0*.bed')) > 0:  #ensure files exist before cat
                    if chromosome[0] not in excluded_chromosomes:
                        output_beds[str(chromosome[0])].append(output_name)
                        subprocess.call(command,shell=True)

                input_list = 'input/' + str(chromosome[0]) + '.' + str(i) + '_' + tag + '_' + sample_inst + '0*pvalues.counts'
                output_name = 'input/' + str(chromosome[0]) + '.' + str(i) + '_' + sample_inst + '_' + tag + 'finalpvalues.counts'
                command = 'cat '+ input_list + ' > ' + output_name
                if len(glob.glob('input/' + str(chromosome[0]) + '.' + str(i) + '_' + tag + '_' + sample_inst + '0*pvalues.counts')) > 0: #ensure files exist
                    if chromosome[0] not in excluded_chromosomes:
                        subprocess.call(command,shell=True)


    for sample_inst in samples:
        for chromosome in unique_results:
                print("concatenating final bed files across super-regions")
                concatenate_bed = 'input/' + str(chromosome[0]) + '.' + '*_' + sample_inst  + '_' + tag + 'final.bed'
                output_name = 'input/' + str(chromosome[0]) + '_' + sample_inst  + '_' + tag + 'final.bed'
                command = 'cat ' + concatenate_bed + ' > ' + output_name
                if len(glob.glob('input/' + str(chromosome[0]) + '.' + '*_' + sample_inst  + '_' + tag + 'final.bed')) > 0:
                     subprocess.call(command,shell=True)



    for sample_inst in samples:
        for chromosome in unique_results:
            for i in range(1,chromosome[1] + 1):
                print("removing chromosome meta files")

                if len(glob.glob('input/' + str(chromosome[0]) + '.' + str(i) + '_' + tag + '_' + sample_inst + '0*.bed')) > 0:  #ensure files exist before cat

                    remove_bed = 'input/' + str(chromosome[0]) + '.' + str(i) + '_' + tag + '_' + sample_inst + '0*.bed'
                    command_removal = 'rm ' +  remove_bed
                    subprocess.call(command_removal,shell=True)

                if len(glob.glob('input/' + str(chromosome[0]) + '.' + str(i) + '_' + tag + '_' + sample_inst + '0*pvalues.counts')) > 0: #ensure files exist

                    remove_counts = 'input/' + str(chromosome[0]) + '.' + str(i) + '_' + tag + '_' + sample_inst + '0*pvalues.counts'
                    command_removal = 'rm ' + remove_counts
                    subprocess.call(command_removal,shell=True)

                    remove_heatmap = 'input/' + str(chromosome[0]) + '.' + str(i) +  '_' + tag + '_' + sample_inst +  '*_heatmap.csv'
                    command_removal = 'rm ' + remove_heatmap
                    subprocess.call(command_removal,shell=True)

    for chromosome in list(output_beds.keys()):
        if len(output_beds[chromosome]) > 1:
            diff1 = np.genfromtxt(output_beds[chromosome][0]) #get first file
            diff2 = np.genfromtxt(output_beds[chromosome][1]) #get second file
            if len(diff1) != len(diff2):
                print("WARNING: bed sizes do not match")
    for chromosome in chr:
        remove_splice_bed = 'input/' + str(chromosome) + '_' + tag + '_' + settings["bed_file"]
        command_removal = 'rm ' + remove_splice_bed
        subprocess.call(command_removal,shell=True)

    for counts_inst in counts:
        for chromosome in chr:
            remove_splice_counts = 'input/' + str(chromosome) + '_' + tag + '_' + counts_inst
            command_removal = 'rm ' + remove_splice_counts
            subprocess.call(command_removal,shell=True)




def main():

    settings = read_settings(sys.argv[1]) # read in user specified settings file
    tags = create_tags(sys.argv[1])
    tag_preprocess = tags[0]
    scale = settings['scale']

    prelim_genomewide(settings,tag_preprocess,scale)

if __name__=="__main__":
    main()
