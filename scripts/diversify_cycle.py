#!/usr/bin/env python3

# Copyright (C) 2019 Andreas Füglistaler <andreas.fueglistaler@epfl.ch>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import argparse
import sys
import os
import numpy as np

# usage: $ diversify pep_stat.sc > tags

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Diversify structures.')
    parser.add_argument("file", type=argparse.FileType('r'), nargs='?', default=sys.stdin,
                        help="Score file")

    args  = parser.parse_args()
    f    = np.loadtxt(args.file, usecols=tuple(range(19))) 
    print("# Total number of Structures", len(f[..., 0]))

    # Get rid of NaN (not a number), pep stats undetermined
    f  = f[~np.any(np.isnan(f), axis=1)]                          
    print("# Total number of non-nan Structures", len(f[..., 0])) 

    # Initial Bin sizes
    Nbin = 10
    
    s  = f[..., 2]      # I_sc + pep_sc, 3rd col of pep_stat
    c  = f[..., 3:6]    # Center of Mass com1, com2, com3
    w  = f[..., 6:9]    # eigenvalues w1, w2, w3

    # Filter factor (times standard deviation)
    ff = 3

    top = int(0.2*len(s))              # top X % of scores
    ms = np.median(s)                  # median score
    ss = np.std(s)                     # score s.d.
    s_sort = s[np.argsort(s)]          # sort by scores to filter top X %
    #print("Score-sorted array:", s)
    ii = np.where(s <= s_sort[top])[0] # filter score <= 'top' score
    print("# Filtering score      < %6.3f -> %6i"%(s_sort[top], len(ii))) # median score and number of decoys below median score

    vs = []
    aa = []
    for i in range(3):
        istart = 9 + i*3     # i=0 => v11,v12,v13
        istop  = 12 + i*3    # i=1 => v21,v22,v23
        ei     = np.zeros(3) # i=2 => v31,v32,v33
        ei[i]  = 1.          # reference vector (1,0,0) (0,1,0) (0,0,1)

        #v.append(f[..., istart:istop])
        #v[i] /= np.linalg.norm(v[i], axis=1)
        vs.append(np.r_[[vj/np.linalg.norm(vj) for vj in f[..., istart:istop]]]) # normalize all angles with principle axis (v1x)
        aa.append(np.arccos(np.clip(np.abs(np.dot(vs[i], ei)), 0, 1)))           # dot-product of above with reference vector ei
        ma = np.mean(aa[i][ii], axis=0)                                          # mean
        aa[i] -= ma                                                              # deviation from mean

    na = np.sqrt(aa[0]**2 + aa[1]**2 + aa[2]**2)                    # sd across 3 principle axes
    sa = np.std(na)                                                 # sd across all angles
    ii = ii[np.where(na[ii] < ff*sa)[0]]                            # subset of angles less than 3*sd filter factor
    print("# Filtering angle      < %6.3f -> %6i"%(ff*sa, len(ii))) # sd of angles * filter factor, number of decoys within filter of angle and score

    #COM
    mc = np.mean(c[ii], axis=0)                                     # mean COM of angle-filtered decoys
    c -= mc                                                         # deviation from mean COM
    nc = np.linalg.norm(c, axis=1)                                  # normalize
    sc = np.std(nc)                                                 # sd of COM
    ii = ii[np.where(nc[ii] < ff*sc)[0]]                            # subset of COM less than 3*sd filter factor
    print("# Filtering COM        < %6.3f -> %6i"%(ff*sc, len(ii))) # sd of COM * filter factor, number of decoys within filter of COM, and angle/score

    #Eigenvalues
    mw = np.mean(w[ii], axis=0)
    w -= mw
    nw = np.linalg.norm(w, axis=1)
    sw = np.std(nw)
    ii = ii[np.where(nw[ii] < ff*sw)[0]]
    print("# Filtering Eigenvalue < %6.3f -> %6i"%(ff*sw, len(ii))) # sd of eigenvalues * filter factor, number of decoys within filter of eigenvalues, COM/angle/score

    # Statistics
    sa  = np.std(na[ii]) # sd of angles among remaining decoys
    sc  = np.std(nc[ii]) # sd of COM among remaining decoys
    sw  = np.std(nw[ii]) # sd of eigenvalues among remaining decoys

    print("# Filtered standard deviations:")
    print("# Angles    : %6.3f"%sa)
    print("# COM       : %6.3f"%sc)
    print("# Eigenvalue: %6.3f"%sw)

##### Cycle to find binning to produce ~200 diverse structures

    Ng = 200
    Rmin = 0
    Rmax = 20
    for cycle in range(30):
        Npos = Nbin
        Nang = Nbin
        Neig = Nbin
        #Sort by score, prep CoM radius
        ii  = ii[np.argsort(s[ii])]
        #print("score-ordered decoys: ", ii)
        R   = ff*sc
        Rmi = 0.
        Rma = -1.

##### 1st binning, CoM
        for n in range(10):
            gi = np.array([ii[0]], dtype=int)                                 # add top-scoring decoy to diversified subset
            for i in ii[1:]:                                                  # for every filtered decoy
                nc = np.linalg.norm(c[gi] - c[i], axis=1)                     # sd angle / sd COM
                t  = np.sqrt(nc*nc)                                           # "radius" of difference between two decoys for all features
                if not np.any(t < R):
                    gi = np.append(gi, i)                                     # t>R, then decoy added to diversified subset
            N = len(gi)
            print("# CoM Iteration %2i, R = %6.3f, N = %6i"%(n+1, R, N))
            if abs(N/Npos - 1) < 0.02:
                break
            if N < Npos:
                Rma = R
            else:
                Rmi = R
            if Rma > 0:
                R = (Rmi+Rma)/2
            else:
                R *= 1.5
        
        Npos = len(gi)
        CoM_bins = [[] for x in range(Npos)]
        
        for i in ii:
            bin_dist = []
            for j in gi:
                nc = c[j] - c[i]
                dist = np.sqrt(nc[0]*nc[0] + nc[1]*nc[1] + nc[2]*nc[2])
                bin_dist.append(dist)
            min_dist = min(bin_dist)
            closest_bin = bin_dist.index(min_dist)
            CoM_bins[closest_bin].append(i)
        
        # Sort by score, prep angle radius
        #ii  = ii[np.argsort(s[ii])] # already sorted
        R   = ff*sw
        Rmi = 0.
        Rma = -1.

##### 2nd binning, shape
        for n in range(10):
            gi = np.array([ii[0]], dtype=int)                                 # add top-scoring decoy to diversified subset
            for i in ii[1:]:                                                  # for every filtered decoy
                nw = np.linalg.norm(w[gi] - w[i], axis=1)                     # eigenvalue comparison
                t  = np.sqrt(nw*nw)                                           # "radius" of difference between two decoys for all features
                if not np.any(t < R):
                    gi = np.append(gi, i)                                     # t>R, then decoy added to diversified subset
            N = len(gi)
            print("# Iteration %2i, R = %6.3f, N = %6i"%(n+1, R, N))
            if abs(N/Neig - 1) < 0.02:
                break
            if N < Neig:
                Rma = R
            else:
                Rmi = R
            if Rma > 0:
                R = (Rmi+Rma)/2
            else:
                R *= 1.5
        
        Neig = len(gi)
        CoM_Eigen_bins = []
        
        for bin in CoM_bins:
            Eigen_bins = [[] for x in range(Nang)]
            for i in bin:
                bin_eig = []
                for j in gi:
                    nw = w[j] - w[i]
                    diff = np.sqrt(nw[0]*nw[0] + nw[1]*nw[1] + nw[2]*nw[2])
                    bin_eig.append(diff)
                min_diff = min(bin_eig)
                closest_bin = bin_eig.index(min_diff)
                Eigen_bins[closest_bin].append(i)
            CoM_Eigen_bins = CoM_Eigen_bins + Eigen_bins

        # Sort by score, prep eigenvalue radius
        #ii  = ii[np.argsort(s[ii])] # already sorted
        R   = ff*sa
        Rmi = 0.
        Rma = -1.
        
##### 3rd binning, angle
        for n in range(10):
            gi = np.array([ii[0]], dtype=int)                                      # add top-scoring decoy to diversified subset
            for i in ii[1:]:                                                       # for every filtered decoy
                a0 = np.arccos(np.clip(np.abs(np.dot(vs[0][gi], vs[0][i])), 0, 1))
                a1 = np.arccos(np.clip(np.abs(np.dot(vs[1][gi], vs[1][i])), 0, 1))
                a2 = np.arccos(np.clip(np.abs(np.dot(vs[2][gi], vs[2][i])), 0, 1))
                t  = np.sqrt(a0*a0 + a1*a1 + a2*a2)                                # "radius" of difference between two decoys for all features
                if not np.any(t < R):
                    gi = np.append(gi, i)                                          # t>R, then decoy added to diversified subset
            N = len(gi)
            print("# Angle Iteration %2i, R = %6.3f, N = %6i"%(n+1, R, N))
            if abs(N/Nang - 1) < 0.02:
                break
            if N < Nang:
                Rma = R
            else:
                Rmi = R
            if Rma > 0:
                R = (Rmi+Rma)/2
            else:
                R *= 1.5
        
        Nang = len(gi)
        CoM_Eigen_Angle_bins = []
        
        for bin in CoM_Eigen_bins:
            Angle_bins = [[] for x in range(Nang)]
            for i in bin:
                bin_ang = []
                for j in gi:
                    a0 = np.arccos(np.clip(np.abs(np.dot(vs[0][j], vs[0][i])), 0, 1))
                    a1 = np.arccos(np.clip(np.abs(np.dot(vs[1][j], vs[1][i])), 0, 1))
                    a2 = np.arccos(np.clip(np.abs(np.dot(vs[2][j], vs[2][i])), 0, 1))
                    theta = np.sqrt(a0*a0 + a1*a1 + a2*a2)
                    bin_ang.append(theta)
                min_theta = min(bin_ang)
                closest_bin = bin_ang.index(min_theta)
                Angle_bins[closest_bin].append(i)
            CoM_Eigen_Angle_bins = CoM_Eigen_Angle_bins + Angle_bins

##### Final low-score diverse tags
        final_bins = [x for x in CoM_Eigen_Angle_bins if x != []]
        #print("# Final bins: ", final_bins)
        t  = np.loadtxt(args.file.name, dtype=str, usecols=19) # tag column of pep_stat.sc
        
        diverse = []
        for bin in final_bins:
            diverse.append(bin[0])  #first element is best by score due to initial ii sort by score

        N = len(diverse)
        print("# Binning Cycle %2i, bins = %6.3f, N = %6i"%(cycle+1, Nbin, N))
        #print("# ", Rmin, Rmax)
        print("# ")
        if 0 <= (N/Ng - 1) <= 0.1:
            break
        if Nbin == Rmax:
            break
        if Rmin == Rmax-1:
            Nbin = Rmax
        if N < Ng:
            Rmin = Nbin
            Nbin = int((Rmax+Nbin)/2)
        else:
            Rmax = Nbin
            Nbin = int((Rmin+Nbin)/2)


##### End binning cycle

    for n in range(len(final_bins)):
        print("# Bin ", n)
        for tag in final_bins[n]:
            print("# ", t[tag], s[tag])
            
    # Diversified tags
    t  = np.loadtxt(args.file.name, dtype=str, usecols=19) # tag column of pep_stat.sc
    print("# %d/%d"%(len(diverse), len(ii)))
    #print("#$ ospymol ", end = '')
    #for tag in t[diverse]:
    #    print(tag, end = '.pdb ')
    #print("")
    for tag in t[diverse]:
        print(tag)

    # Plot score v rmsd
    import matplotlib.pyplot as plt
    rms = f[..., 18]
    plt.plot(rms, s, "o", markersize=6, label="All")
    plt.plot(rms[ii], s[ii], "s", markersize=6, label="Filtered")
    #plt.plot(rms, s, "o", markersize=6)
    plt.plot(rms[diverse], s[diverse], "X", markersize=5, label="Diversified")
    mi = np.argmin(rms)
    plt.plot([rms[mi], rms[mi]], [-1e9, 1e9], "k:", label="Best rms")

    sdiff = s.max() - s.min()
    plt.xlim([0.95*rms.min(), 1.05*rms.max()])
    plt.ylim([s.min() - 0.05*sdiff, s.max() + 0.05*sdiff])

    #sdiff = max(0, ms, s[np.argmin(rms)]) - s.min()
    #plt.xlim([0, 5])
    #plt.ylim([s.min() - 0.05*sdiff, max(0, ms, s[np.argmin(rms)]) + 0.05*sdiff])

    plt.legend()

    #name = args.file.name.replace("output_", "").replace("/", "_").split(".")[0] + ".png"
    name = "rms_v_score.png"
    fig  = plt.gcf()

    fig.tight_layout()
    fig.set_size_inches(12.80, 7.20)
    #print("# Writing", name)
    fig.savefig(name, dpi=100)
    os.system("convert -trim -bordercolor White " + name + " " + name)

