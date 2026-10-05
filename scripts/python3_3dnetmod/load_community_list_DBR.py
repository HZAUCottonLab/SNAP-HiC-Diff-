
# Community results must be in one of the following two formats:
#Region Gamma   Community_Assignment    Chr     Start_Pos       Stop_Pos        Left_Variance   Right_Variance
#Chr    Upstream_boundary_start Upstream_boundary_stop  Downstream_boundary_start       Downstream_boundary_stop


def load_community_list(results_path):
        input = open(results_path, 'r')
        communities = [] # this will be a list of dicts

        for line in input:
                elements = line.strip().split('\t')
                test = len(elements)
                if line.startswith('#'):
                        elements = line.strip().split('\t')
                        test = len(elements)
                        print('test: ', test)
                else:
                        if test == 8:
                                community = {}
                                comm_name = line.strip().split('\t')[2]
                                community['chr'] = line.strip().split('\t')[3]
                                community['start'] = line.strip().split('\t')[4]
                                community['end'] = line.strip().split('\t')[5]
                                communities.append(community)
                        elif test == 5:
                                community = {}
                                community['chr'] = line.strip().split('\t')[0]
                                community['start1'] = line.strip().split('\t')[1]
                                community['start2'] = line.strip().split('\t')[2]
                                community['start'] = line.strip().split('\t')[1]
                                community['end1'] = line.strip().split('\t')[3]
                                community['end2'] = line.strip().split('\t')[4]
                                community['end'] = line.strip().split('\t')[4]
                                communities.append(community)
                        elif test == 3:
                                community = {}
                                community['chr'] = line.strip().split('\t')[0]
                                community['start'] = line.strip().split('\t')[1]
                                community['end'] = line.strip().split('\t')[2]
                                communities.append(community)
                        else:

                                print('Community file is not in expected format!')
                                return None
        input.close()

        return communities
