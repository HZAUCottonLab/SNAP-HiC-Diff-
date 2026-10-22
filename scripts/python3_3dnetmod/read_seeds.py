"""
This module reads in a .txt file of random seeds and returns a list of integer seeds
"""
import os

def read_seeds(seed_file, num_part):
        if not os.path.isabs(seed_file) and not os.path.isfile(seed_file):
                code_dir = os.environ.get('PY3_CODEDIR', os.path.dirname(os.path.abspath(__file__)))
                candidate = os.path.join(code_dir, seed_file)
                if os.path.isfile(candidate):
                        seed_file = candidate
        seeds = []
        input = open(seed_file, 'r')
        num_part = int(num_part)
        for line in input:
                seed = int(line.strip())
                if seed not in seeds:
                        seeds.append(seed)

        if num_part <=100:
                seeds = seeds[0:num_part]
        input.close()
        print("number of seeds",  len(seeds))
        return seeds


def main():
        seed_file = 'gen_louvain_seed.txt'
        seeds = read_seeds(seed_file)
        filename = 'gen_louvain_100_seeds.txt'
        seeds_to_write = open(filename, 'w')
        for i in range(100):
                print(seeds[i], file=seeds_to_write)
        seeds_to_write.close()


if __name__ == '__main__':
        main()
