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
                for line_number, line in enumerate(input, start=1):
                        record = line.rstrip('\r\n')
                        if not record.strip() or record.lstrip().startswith('#'):
                                continue # ignore empty lines and comments
                        if '\t' in record:
                                pair = record.split('\t')
                        else:
                                pair = record.split()
                        if len(pair) != 2 or not pair[0].strip():
                                raise ValueError(
                                        f"malformed setting on line {line_number} in {filename}: "
                                        "expected a key and a value"
                                )
                        var_dict[pair[0].strip()] = pair[1].strip() # key = category, value = stored string

        return var_dict
