
def load_unqueriable(filename):
        regions = []
        input = open(filename, 'r')
        for line in input:
                if line.startswith('#'):
                        continue
                else:
                        pieces = line.strip().split('\t')
                        region = {}
                        region['chrom'] =  pieces[0]
                        region['start'] = int(pieces[1])
                        region['end'] = int(pieces[2])
                        regions.append(region)

        input.close()
        return regions
