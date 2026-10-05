
# File should be in format
# chr \t start \t stop
# will be returned as list of dictionaries with keys 'chrom', 'start' 'end'
# if in format chr \t start1 \t start2 \t end1 \t end2
# interval will be start1 to end2

def load_intervals(results_path):
        input = open(results_path, 'r')
        communities = [] # this will be a list of dicts

        for line in input:
                if line.startswith('#'):
                        elements = line.strip().split('\t')
                        test = len(elements)
                else:
                        if test == 5:
                                community = {}
                                community['chrom'] = line.strip().split('\t')[0]
                                community['start'] = int(line.strip().split('\t')[1])
                                community['end'] = int(line.strip().split('\t')[4])
                                communities.append(community)
                        elif test == 3:
                                community = {}
                                community['chrom'] = line.strip().split('\t')[0]
                                community['start'] = int(line.strip().split('\t')[1])
                                community['end'] = int(line.strip().split('\t')[2])
                                communities.append(community)
                        else:

                                print('Community file is not in expected format!')
                                return None
        input.close()

        return communities


def load_stripped_intervals(results_path):
        ### Assumes chromosome is in format 'chr18.1'
        input = open(results_path, 'r')
        communities = [] # this will be a list of dicts

        for line in input:
                if line.startswith('#'):
                        elements = line.strip().split('\t')
                        test = len(elements)
                else:
                        if test == 5:
                                community = {}
                                community['chrom'] = line.strip().split('\t')[0].split('.')[0]
                                community['start'] = int(line.strip().split('\t')[1])
                                community['end'] = int(line.strip().split('\t')[4])
                                communities.append(community)
                        elif test == 3:
                                community = {}
                                community['chrom'] = line.strip().split('\t')[0].split('.')[0]
                                community['start'] = int(line.strip().split('\t')[1])
                                community['end'] = int(line.strip().split('\t')[2])
                                communities.append(community)
                        else:

                                print('Community file is not in expected format!')
                                return None
        input.close()

        return communities
