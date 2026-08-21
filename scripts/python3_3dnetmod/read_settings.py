def read_settings(filename):
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

        var_dict = {}

        with open(filename, 'r') as input:
                for line in input:
                        if line.startswith('#') or not line.strip():
                                continue # ignore empty lines and comments
                        pair = line.strip().split('\t')
                        if len(pair) == 1: #test to see if user put space separatation
                                a = line.strip().split(' ')
                                pair = list([x for x in a if x!= ''])
                        var_dict[pair[0].strip()] = pair[1].strip() # key = category, value = stored string

        return var_dict
