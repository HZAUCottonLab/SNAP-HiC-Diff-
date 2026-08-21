import json
import random
import numpy as np
import argparse
#import null_model_und_signc as nmus #c version not in use
from randmio_und_signed import randmio_und_signed

'''
        NULL_MODEL_UND_SIGN     Random graphs with preserved weight, degree and
                        strength distributions

        W0 = null_model_und_sign(W);
        W0 = null_model_und_sign(W,bin_swaps);
        W0 = null_model_und_sign(W,bin_swaps,wei_freq);
        [W0 R] = null_model_und_sign(W,bin_swaps,wei_freq);

        This function randomizes an undirected network with positive and
        negative weights, while preserving the degree and strength
        distributions. This function calls randmio_und_signed.m

        Inputs: W,          Undirected weighted connection matrix
                bin_swaps,  Average number of swaps of each edge in binary randomization.
                           bin_swap=5 is the default (each edge rewired 5 times)
                           bin_swap=0 implies no binary randomization
                wei_freq,   Frequency of weight sorting in weighted randomization
                           wei_freq must be in the range of: 0 < wei_freq <= 1
                           wei_freq=1 implies that weights are resorted at each step
                                   (default in older [<2011] versions of MATLAB)
                           wei_freq=0.1 implies that weights are sorted at each 10th step
                                   (faster, default in newer versions of Matlab)

        Output:     W0,     Randomized weighted connection matrix
                    R,      Correlation coefficient between strength sequences
                           of input and output connection matrices

        Notes:
                The value of bin_swaps is ignored when binary topology is fully
        connected (e.g. when the network has no negative weights).
                Randomization may be better (and execution time will be slower) for
        higher values of bin_swaps and wei_freq. Higher values of bin_swaps may
        enable a more random binary organization, and higher values of wei_freq
        may enable a more accurate conservation of strength sequences.
                R are the correlation coefficients between positive and negative
        strength sequences of input and output connection matrices and are
        used to evaluate the accuracy with which strengths were preserved. Note
        that correlation coefficients may be a rough measure of
        strength-sequence accuracy and one could implement more formal tests
        (such as the Kolmogorov-Smirnov test) if desired.

        Example usage:

        %Create random weights matrix

        W=tril(randn(100),-1); W=W+W.';

        %Compute one instance of null model (slow execution time):
        %bin_swaps=5,   rewire each binary edge 5 times on average
        %wei_freq=1,    sort all edges at every step

        tic; [W0_slow R_slow]=null_model_und_sign(W,5,1); R_slow, toc

        R_slow =
                0.9720    0.9746
        Elapsed time is 2.112715 seconds.

        %Compute another instance of of null model (fast execution time):
        %bin_swaps=5,   rewire each binary edge 5 times on average
        %wei_freq=0.1,  sort all edges at every 10th step (10=1/0.1)

        tic; [W0_fast R_fast]=null_model_und_sign(W,5,0.1); R_fast, toc

        R_fast =
                0.9623    0.9789
        Elapsed time is 0.584797 seconds.


        Reference: Rubinov and Sporns (2011) Neuroimage 56:2068-79


        2011-2015, Mika Rubinov, U Cambridge

        Modification History
        Mar 2011: Original
        Sep 2012: Edge-sorting acceleration
        Dec 2015: Enforce preservation of negative degrees in sparse
             networks with negative weights (thanks to Andrew Zalesky).

        #ok<*ASGLU>
'''

def null_model_und_sign(HM,swaps,bin_swaps,wei_freq=1):

     #always use seed 31415 - same random rewiring for every attempt (ensure same max gamma is reproduced)
     random.seed(31415)

     n = HM.shape[0]  # number of nodes
     np.fill_diagonal(HM, 0) #clear diagonal
     Ap = HM > 0
     Ap = Ap*np.ones(HM.shape)              #positive adjacency matrix
     #print "Ap = ",Ap
     #print "len(Ap) = ",len(Ap)

     An = HM < 0
     An = An*np.ones(HM.shape)              #negative adjacency matrix

     #print "swaps = ", swaps
     #print "eval(swaps) = ", eval(swaps)
     if ((len(HM[HM> 0]) < n*(n-1)) and eval(swaps)):         #if Ap is not full (not implementing sparse matrix)
           #W_r  = nmus.randmio_und_signed(HM,bin_swaps); #c code and wrapper version of swaps
           W_r = randmio_und_signed(HM,bin_swaps) #python version of swaps
           #print "W_r[W_r > 0] = ",W_r[W_r>0]
           print("swapping turned on")
           Ap_r = W_r>0
           Ap_r = Ap_r*np.ones(W_r.shape)
           Ap_r = np.squeeze(np.array(Ap_r))

           An_r = W_r<0
           An_r = An_r*np.ones(W_r.shape)
           An_r=np.squeeze(np.array(An_r))
     elif ((len(HM[HM> 0]) >= n*(n-1)) and eval(swaps)):
           print('matrix too dense for swapping')
           W_r = HM
           Ap_r = Ap
           An_r = An
     else:
          print('swapping turned off')
          W_r = HM
          Ap_r = Ap
          An_r = An

     HM0 = np.zeros(HM.shape)                    #null model network
     #HM0 = HM0.flatten(order='F')

     for s in range(1,4,2):
           I = []
           J = []
           if s==1:
              #  S = np.sum(HM*Ap,axis=1)                    #positive strength
                S = np.sum(W_r*Ap_r,axis=1)
              #  pre = np.triu(Ap)*HM
                pre = np.triu(Ap_r)*W_r
                #print "pre  =",pre
                pre2 = pre[pre>0]
                wv = np.sort(pre2,kind='mergesort')         #sorted weights vector
                pre3 = np.triu(Ap_r*np.ones(HM.shape))
                pre3 = pre3.flatten(order='F')
                lij = np.where(pre3>0)                     #linear weights indices
                for i in range(0,len(lij)):
                    indexi = lij[i] % len(HM)
                    # Python 3: use integer division for indices
                    indexj = lij[i] // len(HM)
                    I.append(indexi)
                    J.append(indexj)
           elif s==3:
             #   S = np.sum(-1*HM*An,axis=1)                    #negative strength
                S = np.sum(-1*W_r*An_r,axis=1)                  #negative strength
             #   pre = np.triu(An)*HM
                pre = np.triu(An_r)*W_r
                pre2 = pre[pre<0]
                wv = np.sort(pre2,kind='mergesort')        #sorted weights vector
                pre3 = np.triu(An_r*np.ones(HM.shape))
                pre3 = pre3.flatten(order='F')
                lij = np.where(pre3>0)
                for i in range(0,len(lij)):
                    indexi = lij[i] % len(HM)
                    # Python 3: ensure integer index for column
                    indexj = lij[i] // len(HM)
                    I.append(indexi)
                    J.append(indexj)
           #print "over"

           I = np.array(I, dtype=int)
           J = np.array(J, dtype=int)
           p = np.outer(S,S.T)                   #expected weights matrix
           pflat = p.flatten(order='F')

           #print "over2"
           I = I[0]
           J = J[0]

           wei_period = np.round(1/wei_freq)      #convert frequency to period
           for m in range(len(wv)-1,-1,-1*int(wei_period)):       #iteratively explore at the given period
                #print "wei_freq != 1,m = ",m
                R = []
                dum = np.sort(pflat[lij])#,kind='mergesort')       #get indices of Lij that sort P
                Oind = np.argsort(pflat[lij])#,kind='mergesort')
                chosen = random.randint(0,m)
                R.append(chosen)
                #print "m = ",m
                #print "len(Oind) = ",len(Oind)
                while (len(R) < wei_period):  #inf loop if m==0
                     chosen = random.randint(0,m)
                     if ((chosen not in R) and (chosen < len(Oind))):
                          R.append(chosen)
                     if (m < wei_period) and len(R) >= m:
                          break

                R = np.squeeze(np.array(R))
                O = Oind[R]

                HM0[I[O],J[O]] = s*wv[R]   #assign corresponding sorted weight at this index

                if m <= int(wei_period):
                     break;

                WA = np.zeros(n)
                totalwv = np.concatenate([wv[R],wv[R]])
                subscripts = np.concatenate([I[O],J[O]])
                for u in range(0,n):
                    if u in subscripts:
                        lookup = np.where(subscripts==u)[0]
                        for v in range(0,len(lookup)):
                              WA[u] = WA[u] + totalwv[lookup[v]]  #%cumulative weight
                IJu = WA != 0
                IJu = IJu*np.ones(len(IJu))
                F = 1-WA[np.nonzero(WA)]/S[np.nonzero(WA)]
                adjustments = np.where(IJu==1)[0]
                for k in range(0,len(adjustments)):
                     p[adjustments[k],:] = p[adjustments[k],:]*F[k]
                     p[:,adjustments[k]] = p[:,adjustments[k]]*F[k]
                     S[adjustments[k]] = S[adjustments[k]]-WA[adjustments[k]]    #re-adjust strengths of nodes I(o) and J(o)

                O = Oind[R]
                lij = np.delete(lij,O)        #remove current index from further consideration
                I = np.delete(I,O)
                J = np.delete(J,O)
                wv = np.delete(wv,R)          #remove current weight from further consideration

     HM0 = HM0 + HM0.T

     v1 = np.sum(HM*(HM>0),axis=0)
     v2 = np.sum(HM0*(HM0>0),axis=0)
     rpos = np.corrcoef(v1,v2)

     v1 = np.sum(-1*HM*(HM<0),axis=0)
     v2 = np.sum(-1*HM0*(HM0<0),axis=0)
     rneg = np.corrcoef(v1,v2)

     return HM0,rpos,rneg


def main():

     parser = argparse.ArgumentParser()
     parser.add_argument("file", type=str, help="heatmap")
     parser.add_argument("no", type=str, help="copies")
     parser.add_argument("flips", type=str, help="turn on/off swaps between 0 and nonzero")
     args = parser.parse_args()

     if args.flips == "False":
          swap = False
     else:
          swap = True

     print("swap = ", swap)

     HM = np.genfromtxt(args.file,delimiter=',')

     for i in range(0,int(args.no)):
          #HM = np.genfromtxt("../grand_sim5.csv",delimiter=',')
          HM0 = null_model_und_sign(HM,swap,wei_freq=0.1)
          #print "HM0 = ", HM0
          #np.savetxt("input/random_networks/matrix_sim5no1_pNPC1_S4_nan_smoothed_obs_trimmed.counts.csv_random_" + str(i) + str(args.flips) + "py.csv",HM0,delimiter=',')
          #np.savetxt("input/random_networks/matrix_nestin2_pNPC1_S4_nan_smoothed_obs_trimmed.counts.csv_random_" + str(i) + str(args.flips) + "py.csv",HM0,delimiter=',')
          np.savetxt("input/random_networks/Dixon_cortex_Olig1-Olig2_40kb_5C_num_nodes_pvalues_log_2.csv_random_" + str(i) + str(args.flips) + "py.csv",HM0,delimiter=',')
if __name__ == "__main__":
     main()
