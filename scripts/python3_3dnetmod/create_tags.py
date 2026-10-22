def create_tags(filename):
        """
        General settings file reader used to provide other scripts with user-edited inputs from a text file

        INPUTS
        ------
        filename, str
                the filename and location of the settings file. The file should have the format of a category tab separated from string value

        OUTPUT
        ------
        var_dict, dict
                a dictionary constructed from the settings file where key = category and value = the stored values
        """

        preprocess_tag = ''
        gps_tag = ''

        #preprocess tag creation
        with open(filename, 'r') as input:
                for line in input:
                        if line.startswith('#GPS#'):
                                break # ignore empty lines and comments
                        if line.startswith('bed_file'):
                                continue
                        if line.startswith('sample_'):
                                continue
                        if line.startswith('counts_'):
                                continue
                        if line.startswith('processors'):
                                continue
                        pair = line.strip().split('\t')
                        if len(pair) == 1: #test to see if user put space separatation
                                a = line.strip().split(' ')
                                pair = list([x for x in a if x!= ''])
                        #print "pair = ", pair
                        if len(pair) > 1:
                                preprocess_tag  = preprocess_tag + pair[1] +'_'
        preprocess_tag = preprocess_tag[:-1]
        nextup = True
        gps_tag = preprocess_tag

        #gps tag creation
        with open(filename, 'r') as input:
                for line in input:
                        if line.startswith('#GPS#'):
                                nextup = False
                                continue
                        if line.startswith('switch'):
                                continue
                        if line.startswith('seed_file'):
                                continue
                        if line.startswith('hierarchies'):
                                continue
                        if line.startswith('plots'):
                                continue
                        if nextup:
                                continue
                        if line.startswith('#MMCP#'):
                                break
                        pair = line.strip().split('\t')
                        if len(pair) == 1: #test to see if user put space separatation
                                a = line.strip().split(' ')
                                pair = list([x for x in a if x!= ''])
                        if len(pair) > 1:
                                gps_tag  = gps_tag + '_' +  pair[1]

        nextup = True
        mmcp_tag = gps_tag

        #mmcp tag creation
        with open(filename, 'r') as input:
                for line in input:
                        if line.startswith('#MMCP#'):
                                nextup = False
                                continue
                        if line.startswith('switch'):
                                continue
                        if line.startswith('seed_file'):
                                continue
                        if line.startswith('hierarchies'):
                                continue
                        if line.startswith('plots'):
                                continue
                        if nextup:
                                continue
                        if line.startswith('#HSVM#'):
                                break
                        pair = line.strip().split('\t')
                        if len(pair) == 1: #test to see if user put space separatation
                                a = line.strip().split(' ')
                                pair = list([x for x in a if x!= ''])
                        if len(pair) > 1:
                                mmcp_tag  = mmcp_tag + '_' +  pair[1]


        nextup = True
        hsvm_tag = mmcp_tag
        size_list = []
        var_list = []

        #hsvm tag creation
        with open(filename, 'r') as input:
                for line in input:
                        if line.startswith('#HSVM#'):
                                nextup = False
                                continue
                        if nextup:
                                continue
                        if line.startswith('#DBR#'):
                                break
                        if line.startswith('size_s'):
                                line = line[5:]
                                pair = line.strip().split('\t')
                                if len(pair) == 1: #test to see if user put space separatation
                                        a = line.strip().split(' ')
                                        pair = list([x for x in a if x!= ''])
                                if len(pair) > 1:
                                        hsvm_tag = hsvm_tag + '_' + pair[0] + '_' + str(int(pair[1])/1000) + 'kb'
                        elif line.startswith('var_t'):
                                line = line[4:]
                                pair = line.strip().split('\t')
                                if len(pair) == 1: #test to see if user put space separatation
                                        a = line.strip().split(' ')
                                        pair = list([x for x in a if x!= ''])
                                if len(pair) > 1:
                                        hsvm_tag = hsvm_tag + '_v' + pair[0][-1] + '_' + pair[1]
                        else:
                                pair = line.strip().split('\t')
                                if len(pair) == 1: #test to see if user put space separatation
                                        a = line.strip().split(' ')
                                        pair = list([x for x in a if x!= ''])
                                if len(pair) > 1:
                                        hsvm_tag  = hsvm_tag + '_' +  pair[1]

        nextup = True
        dbr_tag = hsvm_tag

        #dbr tag creation
        with open(filename, 'r') as input:
                for line in input:
                        if line.startswith('#DBR#'):
                                nextup = False
                                continue
                        if nextup:
                                continue
                        pair = line.strip().split('\t')
                        if len(pair) == 1: #test to see if user put space separatation
                                a = line.strip().split(' ')
                                pair = list([x for x in a if x!= ''])
                        if len(pair) > 1:
                                dbr_tag  = dbr_tag + '_' +  pair[1]
        return preprocess_tag,gps_tag,mmcp_tag,hsvm_tag,dbr_tag
