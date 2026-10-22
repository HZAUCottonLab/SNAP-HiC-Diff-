import numpy as np
import random

def randmio_und_signed(HM,iter):

      #always use seed 31415
      random.seed(31415)

      # Python 3: ensure ITER is an integer so it can be used in range()
      ITER = int(iter*len(HM)*(len(HM)-1)/2)
      maxAttempts = len(HM)/2;
      #output = np.zeros((len(HM),len(HM)))
      output = HM.copy()

      eff = 0
      a = 0
      b = 0
      c = 0
      d = 0

      for i in range(0, ITER):
           att=0
           while (att<=maxAttempts):

                a = random.randint(0,len(HM)-1)
                r = random.randint(0,len(HM)-1)
                if (r != a):
                     b = r

                r = random.randint(0,len(HM)-1)
                if ((r != a) and (r != b)):
                     c = r

                r = random.randint(0,len(HM)-1)
                if ((r != a) and (r != b) and (r != c)):
                     d = r


                r0_ab = output[a][b]
                r0_cd = output[c][d]
                r0_ad = output[a][d]
                r0_cb = output[c][b]

                if ((np.sign(r0_ab)==np.sign(r0_cd)) and (np.sign(r0_ad)==np.sign(r0_cb)) and (np.sign(r0_ab) != np.sign(r0_ad))):
                     output[a][d] = r0_ab;
                     output[d][a] = r0_ab;
                     output[c][b] = r0_cd;
                     output[b][c] = r0_cd;
                     output[a][b] = r0_ad;
                     output[b][a] = r0_ad;
                     output[c][d] = r0_cb;
                     output[d][c] = r0_cb;

                     eff = eff + 1
                     break

                att = att + 1

      return output
