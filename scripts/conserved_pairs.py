#!/usr/bin/env python3

# Copyright (C) 2019 Andreas Füglistaler <andreas.fueglistaler@epfl.ch>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import sys
import argparse
import util

groups = {'R': 0, 'H': 0, 'K': 0,
          'D': 1, 'E': 1,
          'S': 2, 'T': 2, 'N' : 2, 'Q' : 2,
          'C': 3, 'U': 3, 'G' : 3, 'P': 3,
          'A': 4, 'V':4, 'I': 4, 'L': 4, 'M' : 4, 'F' : 4, 'Y' : 4, 'W' : 4}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find pairs of conserved residues.')
    parser.add_argument("pdb", type=argparse.FileType('r'), nargs='?', default=sys.stdin,
                        help="Homologue pdb")
    parser.add_argument("--template", "-t", type=argparse.FileType('r'), required=True,
                        help="Template")
    parser.add_argument("--native", "-n", type=argparse.FileType('r'), required=True,
                        help="Native")
    parser.add_argument("--fasta", "-f", required=True, help="Alignement fasta-file")
    parser.add_argument("--cut_off", "-c", default=4., type=float, help="Distance constraints")
    #parser.add_argument("--chain_dist", "-d", default=10., type=float, help="Distance from chain residues")
    parser.add_argument("--print_residues", "-p",action='store_true', default=False,
                        help="Print all conserved residues")

    args  = parser.parse_args()
    fasta = util.parse_fasta(args.fasta)
    tem   = util.parse_pdb(args.template)
    hom   = util.parse_pdb(args.pdb)
    nat   = util.parse_pdb(args.native)[0]
    pep   = hom[1]
    ch    = tem[0]
    hi    = 1  # homologue index
    ti    = 1  # template index
    cons  = []
    grps  = []

    for f in fasta:
        if f[0] == '-':
            if f[1] != '-':
                ti += 1
        elif f[1] == '-':
            hi += 1

        else:
            if f[0] == f[1]:
                if f[1] != ch[ti][0].res:
                    print("Warning! Wrong residue-name in pdb. Expected %s, got %s"%(f[1], ch[ti][0].res))
                cons.append((hi, ti, f[1]))
            if groups[f[0]] == groups[f[1]]:
                grps.append((hi, ti, f[0]+f[1]))
            hi += 1
            ti += 1

    l = len(cons)
    if args.print_residues:
        print("# Concerved residues:")
        print("#  Res hom tem")
        for co in cons:
            print("# %4s%4d%4d"%(co[2], co[0], co[1]))

    pairs = []
    for i in range(l):
        ti_i  = cons[i][1]
        for j in range(i + 1, l):
            ti_j = cons[j][1]
            hi_i = cons[i][0]
            hi_j = cons[j][0]

            if (not hom[0][hi_i]) or (not hom[0][hi_j]):
                continue

            if ti_j - ti_i < 5: # No neighbors
                continue
            dmin = util.dist_min(ch[ti_i], ch[ti_j])
            d    = util.dist_CA(ch[ti_i], ch[ti_j])

            if dmin < args.cut_off:
                if (util.dist_chain(hom[0][hi_i], pep) > args.cut_off and
                    util.dist_chain(hom[0][hi_j], pep) > args.cut_off):
                    pairs.append([cons[i], cons[j], d, -1, util.dist_CA(nat[hi_i], nat[hi_j])])
                else:
                    pairs.append([cons[i], cons[j], d, -2, util.dist_CA(nat[hi_i], nat[hi_j])])

        for g in grps: #range(len(grps)):
            ti_j = g[1]
            hi_i = cons[i][0]
            hi_j = g[0]

            if (not hom[0][hi_i]) or (not hom[0][hi_j]):
                continue

            if abs(ti_j - ti_i) < 5: # No neighbors
                continue
            dmin = util.dist_min(ch[ti_i], ch[ti_j])
            d    = util.dist_CA(ch[ti_i], ch[ti_j])

            if dmin < args.cut_off:
                if (util.dist_chain(hom[0][hi_i], pep) > args.cut_off and
                    util.dist_chain(hom[0][hi_j], pep) > args.cut_off):
                    pairs.append([cons[i], g, d, -3, util.dist_CA(nat[hi_i], nat[hi_j])])
                else:
                    pairs.append([cons[i], g, d, -2, util.dist_CA(nat[hi_i], nat[hi_j])])


    for i, pi in enumerate(pairs):
        if pi[3] == -1 or pi[3] == -3:
            dm = pi[2]
            im = i
            for j in range(i+1, len(pairs)):
                pj = pairs[j]
                if abs(pi[0][0] - pj[0][0]) > 2:
                    break
                if pi[3] != -1 and pi[3] != -3:
                    break
                if abs(pi[1][0] - pj[1][0]) > 2:
                    continue

                if pj[3] > pi[3] or (pj[3] == pi[3] and pj[2] < dm):
                    pi[3] = j
                else:
                    pj[3] = im

        if pi[3] == -1 or pi[3] == -3:
            if abs(pi[4] - pi[2]) < 1:
                s = "AtomPair CA %4d CA %4d FLAT_HARMONIC %6.1f 1.0 0.5 # %s%d - %s%d native: %6.1f"%(
                pi[0][0], pi[1][0], pi[2], pi[0][2], pi[0][1], pi[1][2], pi[1][1], pi[4])
            else:
                s = "#AtomPair CA %4d CA %4d FLAT_HARMONIC %6.1f 1.0 0.5 # %s%d - %s%d native nonconfired: %6.1f"%(
                pi[0][0], pi[1][0], pi[2], pi[0][2], pi[0][1], pi[1][2], pi[1][1], pi[4])
        else:
            s = "#AtomPair CA %4d CA %4d FLAT_HARMONIC %6.1f 1.0 0.4 # %s%d - %s%d, redundant: %d native: %6.1f"%(
                pi[0][0], pi[1][0], pi[2], pi[0][2], pi[0][1], pi[1][2], pi[1][1], pi[3]+1, pi[4])

        print(s)
