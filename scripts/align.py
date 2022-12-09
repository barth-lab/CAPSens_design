#!/usr/bin/env python3

# Copyright (C) 2019 Andreas Füglistaler <andreas.fueglistaler@epfl.ch>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import argparse
from pymol import cmd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Align two pdbs with respect with peptide-binding site')
    parser.add_argument("input", help="Input pdb to align")
    parser.add_argument("--to", "-t", required=True, help="Reference Pdb")
    parser.add_argument("--output", "-o", help="Output name. Default is overwriting input.")
    parser.add_argument("--distance", "-d", default=15, type=int,
                        help="Distance from peptide to consider for binding-site")
    parser.add_argument("--receptor", "-r", default="A",
                        help="Receptor chain name. Default is A")
    parser.add_argument("--peptide", "-p", default="B",
                        help="Peptide chain name. Default is B")

    args = parser.parse_args()

    if args.output:
        fn = args.output
    else:
        fn = args.input
    
    cmd.load(args.input, "I")
    cmd.load(args.to, "R")
    cmd.align("I", "R")
    if args.distance > 0:
        print("aligning pocket")
        for i in range(10):
            cmd.align("I and chain %s within %i of I and chain %s"%(args.receptor, args.distance,
                                                                    args.peptide),
                        "R and chain %s within %i of I and chain %s"%(args.receptor, args.distance,
                                                                    args.peptide))
    cmd.save(fn, "I")
