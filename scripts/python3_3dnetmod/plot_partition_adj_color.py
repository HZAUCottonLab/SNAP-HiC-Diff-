import os
from outdir_helper import OUTDIR
import matplotlib.pyplot as plt
from get_partition_colormap import get_colormap
import numpy as np
import itertools

def rev_map(m):
        r = {}
        for k, v in list(m.items()):
                if v not in r:
                        r[v] = []
                r[v].append(k)

        return r

def plot_partition(Ci, consensus, cmap, filename, vmin=0, vmax=None, subdir='', vlines=None):
        """ Plot a partition block Ci with a given cmap, and optional vmin, vmax kwargs for cmap

        Inputs
                Ci: ndarray, n_part by n_nodes
                cmap: matplotlib colormap
                filename: str, filename to save partition block as
                vmin: int/float, min normalization value for cmap
                vmax: int/float, max normalization value for cmap
                subdir: str, subfolder after output/communities_across_numpart/ to save file in
                vlines: iterable ints, x-coordinates of consensus boundary vertical lines

        Outputs
                none
        """


        plt.clf()
        #consensus = consensus[0]
        print('consensus for evaluting:')
        print(consensus)
        consensus = np.asarray(consensus)
        max_value = 2*np.amax(consensus)

        new_partition = np.zeros(Ci.shape)

        #for line in Ci:
        print('Ci: ')
        print(Ci)

        for j in range(len(Ci)):
                line = Ci[j]
                #print line
                communities = [list(g) for k, g in itertools.groupby(line)]
                matches = {} # dictionary where key is community and value is the best matched community number from consensus
                for community in communities:
                        idx_in_community = np.where(line==community[0])
                        idx_in_community = idx_in_community[0].tolist()
                        overlap_dict = {} # keeps track of the overlap each community within the partition has to the others. key is consensus community number
                        consensus_communities = [list(g) for k, g in itertools.groupby(consensus)]
                        #print 'consensus_communities: '
                        #print consensus_communities
                        for consensus_community in consensus_communities:
                                #print 'consensus:'
                                #print consensus
                                #print 'consensus_community: '
                                #print consensus_community
                                idx_in_consensus = np.where(consensus==consensus_community[0])
                                #print 'idx_in_consensus: ', idx_in_consensus
                                idx_in_consensus = idx_in_consensus[0].tolist()
                                overlap = len([w for w in idx_in_consensus if w in idx_in_community])
                                overlap_dict[consensus_community[0]] = overlap
                                #print 'idx_in_consensus: ', idx_in_consensus
                                #print 'overlap: ', overlap
                        best_match = max(overlap_dict, key=overlap_dict.get)
                        matches[community[0]] = best_match

                #print 'consensus: '
                #print consensus
                #print 'line: '
                #print line
                #print 'matches'
                #print matches
                rev_matches = rev_map(matches)
                for community in communities:
                        if community[0] not in rev_matches:
                                rev_matches[community[0]] = []
                #print 'rev_matches'
                #print rev_matches
                max_assignment = max(line)



                for k, v in list(rev_matches.items()):
                        while len(v) > 1:
                                #print 'v: ', v
                                sizes = []
                                assignments = []
                                for i in range(len(v)):
                                        #print 'i', i
                                        sizes.append(len(np.where(line == v[i])[0]))
                                        #print 'size: ', len(np.where(line == v[i])[0])
                                        assignments.append(v[i])
                                        #print 'assignment: ', v[i]
                                idx = np.argmin(sizes)
                                #print 'idx', idx
                                smallest = assignments[idx]
                                #print 'smallest: ', smallest
                                matches[smallest] = max_assignment + 1
                                max_assignment += 1
                                del v[idx]
                                #print 'v', v

                #print 'matches'
                #print matches
                #print 'rev_matches: '
                #print rev_matches

                #print 'max of consensus:'
                #print max(consensus[0])

                #print np.amax(consensus)
                #print '2*np.amax(consensus)'
                #print 2*np.amax(consensus)
                #print 'matches.items()'
                #print matches.items()
                line2 = np.copy(line)
                for k, v in list(matches.items()):
                        #print 'k: ', k
                        #print 'v: ', v
                        #print 'max_value: ', max_value
                        np.place(line2,line2 ==k , v + max_value)
                        #print 'replacing: ', k
                        #print 'with: ', v + max_value
                        #print 'line'
                        #print line
                        #print 'line2: '
                        #print line2
                        #print 'consensus'
                        #print consensus

                #print 'line2 before subtracting: '
                #print line2
                line2 = np.asarray(line2) - max_value
                #print 'line after subtracting '
                #print line2

                new_partition[j:] = line2
        Ci = new_partition

        #print Ci

        # Optionally specify the upper range of colormap
        if vmax is None:
                plt.imshow(Ci, interpolation = 'none', cmap=cmap, vmin=vmin)
                plt.axis('off')
        else:
                plt.imshow(Ci, interpolation = 'none', cmap=cmap, vmin=vmin, vmax=vmax)
                plt.axis('off')

        # plot consensus lines if they are inputted
        if not vlines is None:
                for xcoord in vlines:
                        plt.axvline(x=xcoord, linewidth='2', color='k')

        # Create directory
        dir = OUTDIR + '/MMCP/communities_across_numpart/' + subdir
        if not os.path.isdir(dir): os.makedirs(dir)


        # Save figure
        plt.savefig(os.path.join(dir,filename), dpi = 900, bbox_inches = 'tight')
        plt.clf()
        plt.close()

def main():
        """
        Ci = np.array([[1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 5], /
                       [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3], /
                       [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2], /
                       [1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 3, 4], /
                       [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3]])
        """

        Ci = np.array([[1, 2, 3, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7], [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3], [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2], [1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 3, 4], [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3]])
        #consensus = np.array([[1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3]])
        consensus = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3]
        plot_partition(Ci, consensus, get_colormap(), 'Test.png', vmin=0, vmax=None, subdir='test', vlines=None)



if __name__ == '__main__':
        main()
