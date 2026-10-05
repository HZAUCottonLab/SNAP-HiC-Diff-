

def write_final_consistent_communities(community_list, path):

        out_file = open(path, 'w')
        for element in community_list:
                output = element['chr'] + '\t' + element['start'] + '\t' + element['end']
                print(output, file=out_file)

        out_file.close()
