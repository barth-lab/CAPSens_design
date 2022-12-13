#!/usr/bin/env python3

# Copyright (C) 2019 Andreas Füglistaler <andreas.fueglistaler@epfl.ch>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import sys
import argparse
import numpy as np
from matplotlib.cbook import flatten
import util

def parse_file(fn):
    f  = open(fn)
    ls = []
    for l in f:
        if l and not l[0] == '#':
            ll = l.split("#")[0].split()
            ls.append([int(c) for c in ll])
    f.close()
    return ls

donor_help = """
# Composition file. Examples:
# Use donor-residues 1-12, as aligned:
1 12

# Use donor-residues 1-12, 
# calculate offset at anchor-residue 1:
1  12 1

# Use donor-residues 1-12, 
# calculate offset at anchor-residue 12:
1  12 12 

# Use donor-residues 1-12, 
# calculate offset at anchor-residue 7:
1  12 7

"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compose with multiple templates",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("receiver", type=argparse.FileType('r'), nargs='?', default=sys.stdin,
                        help="Receiver PDB-file")
    parser.add_argument("--donor", "-d", type=argparse.FileType('r'), required=True,
                        help="Donor pdb-file")
    parser.add_argument("--file", "-f", required=True, help=donor_help)
    parser.add_argument("--rm", "-r", action='store_true', default=False,
                        help="Remove residues where the structure is unknown.")

    args = parser.parse_args()

    rs = util.parse_pdb(args.receiver)[0] # Receiver chain
    ds = util.parse_pdb(args.donor)[0]  # Donor chain
    ps = parse_file(args.file) # Parts

    for p in ps:
        start = min(p[0], p[1])
        stop  = max(p[0], p[1])

        off = np.zeros(3)
        ch  = rs[start][1].chain

        if len(p) > 2:
            ank   = p[2] # anker
            off   = rs[ank][1].pos - ds[ank][1].pos  # offset

        for i in range(start, stop+1):
            if rs[i] and ds[i] and not rs[i][0].res == ds[i][0].res:
                print("Warning! Not same residue-names", file=sys.stderr)

            rs[i] = ds[i]
            for ra in rs[i]:
                ra.pos  += off
                ra.chain = ch

    for a in flatten(rs):
        print(str(a))
    print("TER")
