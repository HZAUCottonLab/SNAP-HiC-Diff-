

def load_region_start_end_coord(bed_path):
        def _region_name(bin_name):
                if '_BIN_' in bin_name:
                        return bin_name.rsplit('_BIN_', 1)[0]
                return bin_name.split('_')[0]

        """
        For a given .bed file of genomic regions of interest parsed into bins, this module returns a dictionary where the key = region, and the value is another dictionary with keys 'start'
        and 'end' which are the integer start and stop coordinates of the region of interest. It is assumed that the chromosome is known for a region of interest.
        """
        print(bed_path)
        input = open(bed_path, 'r')
        region_coords = {}
        region_idx = {}
        chr_of_region = {}
        bin_sizes = []
        for line in input:
                if line.startswith('#'):
                        continue
                else:
                        chr = line.strip().split('\t')[0]
                        start = int(line.strip().split('\t')[1])
                        end = int(line.strip().split('\t')[2])
                        size = end - start
                        bin_sizes.append(size)
                        bin_name = line.strip().split('\t')[3]
                        region = _region_name(bin_name)
                        if region not in region_coords:
                                region_coords[region] = []
                                region_idx[region] = {}
                                chr_of_region[region] = chr
                        region_coords[region].append(start)
                        region_coords[region].append(end)
        input.close()


        for region in region_coords:
                region_idx[region]['start'] = min(region_coords[region])
                region_idx[region]['end'] = max(region_coords[region])
                region_idx[region]['chrom'] = chr_of_region[region]

        return region_idx, bin_sizes[0]


def main():
        load_region_start_end_coord('input/RenPooledCortexChr12.001.bed')

if __name__=='__main__':
        main()
