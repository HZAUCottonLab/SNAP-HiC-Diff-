from gamma_plateau_sweep import GPS
from modularity_maximization_consensus_partition import MMCP

'''
   Heidi Norton
   Daniel Emerson
   Harvey Huang

   Determines gamma range per chr region from settings file parameters using GPS (gamma plateau sweep).
   Identifies communities using modularity maximization and consensus partition (MMCP).
'''


def collective(element):
     chr = element[0]
     settings = element[1]
     sample = element[2]

     GPS(chr,settings,sample)
     MMCP(chr,settings,sample)

     return 0

def main():
     return 0

if __name__ == '__main__':
     main()
