def get_gammas(gamma_file):
        """
        Reads from a list of gamma and returns a list of float gammas
        INPUTs
        -----
        gamma_file, str
        name of file containing a separate gamma on each line
        OUTPUTs
        gamma, list of floats
        all the gammas from gamma_file
        """
        gammas = []
        input = open(gamma_file, 'r')
        for line in input:
                if line.startswith('#'):
                        continue
                else:
                        gamma = float(line.strip().split(' ')[0])
                        gammas.append(gamma)
        input.close()
        #print 'gammas: '
        #print gammas

        return gammas
