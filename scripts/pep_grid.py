#!/usr/bin/env python3

# Copyright (C) 2019 Andreas Füglistaler <andreas.fueglistaler@epfl.ch>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import sys
import argparse
from copy import deepcopy
import itertools
import numpy as np
from numpy import r_
import util

def lattice(i, j, k):
    return r_[2*i + ((j + k)%2),
              np.sqrt(3)*(j + (1/3.)*(k%2)),
              2*np.sqrt(6.)/3.*k]

def neighbors(r):
    l = []

    #l.append(lattice(+0, +0, +0))

    l.append(lattice(+0, -1, -1))
    l.append(lattice(-1, +0, -1))
    l.append(lattice(+0, +0, -1))

    l.append(lattice(-1, -1, +0))
    l.append(lattice(+0, -1, +0))
    l.append(lattice(-1, +0, +0))
    l.append(lattice(+1, +0, +0))
    l.append(lattice(-1, +1, +0))
    l.append(lattice(+0, +1, +0))

    l.append(lattice(+0, -1, +1))
    l.append(lattice(-1, +0, +1))
    l.append(lattice(+0, +0, +1))

    return r/2*np.array(l)

def angs(a):
    l = []
    #l.append(r_[0, 0, 0])

    for i in itertools.product([-1, 1], repeat=3):
        l.append(r_[i])

    return a*np.array(l)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Puts peptide on grid around output')
    parser.add_argument("pdb", type=argparse.FileType('r'), nargs='?', default=sys.stdin,
                        help="PDB-file")
    parser.add_argument("--radius", "-r", type=float, default=1, help="Distance from input pdb")
    parser.add_argument("--angle", "-a", type=float, default=30, help="Angle from input pdb")
    parser.add_argument("--rms_max", "-m", type=float, default=5., help="Angle from input pdb")
    parser.add_argument("--prefix", "-p", default="",
                        help="Output prefix. Default is input-file-name")

    args = parser.parse_args()
    pdb  = util.parse_pdb(args.pdb)

    if not args.prefix:
        fn = args.pdb.name.split("/")[-1].split(".pdb")[0]
    else:
        fn = args.prefix

    #Find angles
    a0 = angs(args.angle/180.*np.pi)[0]
    v0 = neighbors(args.radius)[0]
    d  = 1.
    r  = 0.
    for i in range(20):
        p = deepcopy(pdb)
        util.translate(p[1], v0)
        util.rotate(p[1], a0[0], 0)
        util.rotate(p[1], a0[1]/d, 1)
        util.rotate(p[1], a0[2]/d, 2)
        r = util.rmsd(pdb[1], p[1])
        if r < args.rms_max:
            break
        else:
            d += 0.5

    print("Angles: (%6.3f, %6.3f, %6.3f)"%(args.angle, args.angle/d, args.angle/d))
    print("Translation: %6.3f"%args.radius)
    print("Max rmsd: %6.3f"%r)

    i   = 0
    #util.print_pdb(pdb, fn + "_%03d.pdb"%i)
    util.print_pdb(pdb, "%03d.pdb"%i)
    i += 1

    print("Input:%s%s_%03d.pdb"%(20*" ", fn, i-1))
    istart = i

    for v in neighbors(args.radius):
        p = deepcopy(pdb)

        util.translate(p[1], v)
        #util.print_pdb(p, fn + "_%03d.pdb"%i)
        util.print_pdb(p, "%03d.pdb"%i)

        i += 1

    print("Translations:%s%s_%03d.pdb - %s_%03d.pdb"%(13*" ", fn, istart, fn ,i-1))
    istart = i

    for a in angs(args.angle/180.*np.pi):
        p = deepcopy(pdb)

        #util.scaled_rotate(p[1], a)
        util.rotate(p[1], a[0], 0)
        util.rotate(p[1], a[1]/d, 1)
        util.rotate(p[1], a[2]/d, 2)
        #util.print_pdb(p, fn + "_%03d.pdb"%i)
        util.print_pdb(p, "%03d.pdb"%i)

        i += 1
    print("Rotations:%s%s_%03d.pdb - %s_%03d.pdb"%(16*" ", fn, istart, fn, i-1))
    istart = i

    t = 0
    for v in neighbors(args.radius):
        t += 1
        print("Translation %2d:%s%s_%03d.pdb - %s_%03d.pdb"%(t, 11*" ", fn, i, fn, i+7))
        for a in angs(args.angle/180.*np.pi):
            p = deepcopy(pdb)

            util.translate(p[1], v)
            #util.scaled_rotate(p[1], a)
            util.rotate(p[1], a[0], 0)
            util.rotate(p[1], a[1]/d, 1)
            util.rotate(p[1], a[2]/d, 2)
            #util.print_pdb(p, fn + "_%03d.pdb"%i)
            util.print_pdb(p, "%03d.pdb"%i)

            i += 1
    print("Translations + Rotations: %s_%03d.pdb - %s_%03d.pdb"%(fn, istart, fn, i-1))
