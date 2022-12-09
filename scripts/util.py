# Copyright (C) 2019 Andreas Füglistaler <andreas.fueglistaler@epfl.ch>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import sys
import numpy as np
from matplotlib.cbook import flatten

BB   = ["N", "CA", "C", "O"]
BBCB = ["N", "CA", "C", "O", "CB"]

class AtomLine:
    s0 = "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N"
    #sf = "ATOM  %5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  "
    sf = "ATOM  %5i %-4s %3s %1s%4i    %8s%8s%8s%6.2f%6.2f          %2s  "

    three21 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
               'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    def __new__(cls, s=s0):
        if s[17:20] not in AtomLine.three21:
            return None

        return object.__new__(cls)

    def __init__(self, s=s0):
        self.one23 = {v: k for k, v in self.three21.items()}

        self.natom = int(s[6:11])
        self.atom  = s[12:16]
        self.res   = self.three21[s[17:20]]
        self.chain = s[21:22]
        self.nres  = int(s[22:26])
        self.pos   = np.r_[float(s[30:38]), float(s[38:46]), float(s[46:54])]
        self.occ   = float(s[54:60])
        if s[60:66].strip():
            self.temp  = float(s[60:66])
        else:
            self.temp  = 0.

        self.ele   = s[76:78]

    def __str__(self):
        return self.sf%(self.natom, self.atom, self.one23[self.res], self.chain, self.nres,
                        ("%8.3f"%self.pos[0])[0:8], ("%8.3f"%self.pos[1])[0:8],
                        ("%8.3f"%self.pos[2])[0:8], self.occ, self.temp, self.ele)

def c2i(c):
    return ord(c) - ord("A")

# pdb is a list of chains
# chain is a list of residues, starting at 1 (index 0 is empty residued)
# residue is a list of atomlines
def parse_pdb(f):
    pdb = []
    res = 0
    ch  = 0
    chs = ""

    for l in f:
        if not l.startswith("ATOM"):
            continue

        a = AtomLine(l)
        if not a:
            continue

        if a.chain not in chs:
            pdb.append([])
            res  = 0
            chs += a.chain

        ch = chs.find(a.chain)

        if res < a.nres and (len(pdb[ch]) - 1) < a.nres:
            res = len(pdb[ch]) - 1
            while res < a.nres: # add empty list for missing residue
                pdb[ch].append([])
                res += 1

        if res != a.nres:
            res = a.nres

        pdb[ch][res].append(a)
    return pdb

def parse_fasta(fn):
    ls  = []
    f   = open(fn)

    for l in f:
        if not l:
            pass
        elif l[0] == '>':
            ls.append("-") # Start at index 1
        else:
            if not ls:
                ls.append("-")
            ls[-1] += l.strip()
    f.close()
    return list(map(list, zip(*ls)))

def print_pdb(pdb, fn=""):
    if not fn:
        f = sys.stdout
    else:
        f = open(fn, "w")

    for ch in pdb:
        for a in flatten(ch):
            print(str(a), file=f)
        print("TER", file=f)

    if fn:
        f.close()

def ca(p):
    cs = []
    for a in flatten(p):
        if a.atom.strip() == "CA":
            cs.append(a.pos)
    return np.array(cs)

def rmsd(pdb1, pdb2):
    c1 = ca(pdb1)
    c2 = ca(pdb2)
    d  = (c1-c2).flatten()
    return np.sqrt(np.sum(d.dot(d))/len(c1))

def translate(pdb, vec):
    for a in flatten(pdb):
        a.pos += vec

def pca(pdb):
    crds  = ca(pdb)
    com   = np.mean(crds, 0)
    crds -= com
    w, v  = np.linalg.eig(np.dot(crds.T, crds))
    ii    = np.argsort(w)[::-1]
    w     = w[ii]
    v     = v[:, ii].transpose()
    return com, v, w

def rotate(pdb, angle, axis=0):
    crds  = ca(pdb)
    com   = np.mean(crds, 0)
    crds -= com
    w, v  = np.linalg.eig(np.dot(crds.T, crds))
    ii    = np.argsort(w)[::-1]
    w     = w[ii]
    v     = v[:, ii].transpose()
    vr    = v[axis]/np.linalg.norm(v[axis])
    q     = np.r_[np.cos(angle/2), np.sin(angle/2)*vr]
    q_1   = np.r_[np.cos(angle/2), -np.sin(angle/2)*vr]

    for a in flatten(pdb):
        p     = np.r_[0, a.pos - com]
        p_    = quat_mul(quat_mul(q, p), q_1)
        a.pos = p_[1:] + com

def scaled_rotate(pdb, angles):
    crds  = ca(pdb)
    com   = np.mean(crds, 0)
    crds -= com
    w, v  = np.linalg.eig(np.dot(crds.T, crds))
    ii    = np.argsort(w)[::-1]
    w     = w[ii]
    v     = v[:, ii].transpose()

    for a in flatten(pdb):
        a.pos -= com

    for i, wi in enumerate(w):
        ai  = angles[i]*(wi/w[0])**3
        vi  = v[i]/np.linalg.norm(v[i])
        q   = np.r_[np.cos(ai/2), np.sin(ai/2)*vi]
        q_1 = np.r_[np.cos(ai/2), -np.sin(ai/2)*vi]

        for a in flatten(pdb):
            p     = np.r_[0, a.pos]
            p_    = quat_mul(quat_mul(q, p), q_1)
            a.pos = p_[1:]

    for a in flatten(pdb):
        a.pos += com

def quat_mul(q, r):
    t    = np.zeros(4)
    t[0] = r[0]*q[0] - r[1]*q[1] - r[2]*q[2] - r[3]*q[3]
    t[1] = r[0]*q[1] + r[1]*q[0] - r[2]*q[3] + r[3]*q[2]
    t[2] = r[0]*q[2] + r[1]*q[3] + r[2]*q[0] - r[3]*q[1]
    t[3] = r[0]*q[3] - r[1]*q[2] + r[2]*q[1] + r[3]*q[0]
    return t

def dist_fast(r1, r2):
    a1 = r1[0]
    a2 = r2[0]
    return np.linalg.norm(a1.pos - a2.pos)

def dist_CA(r1, r2):
    ca1 = next(filter(lambda a: a.atom.strip() == 'CA', r1))
    ca2 = next(filter(lambda a: a.atom.strip() == 'CA', r2))
    return np.linalg.norm(ca1.pos - ca2.pos)

def dist_min(r1, r2, H=False):
    dmin = 1e12
    for a1 in r1:
        if not H and isH(a1):
            continue
        for a2 in r2:
            if not H and isH(a2):
                continue
            d = np.linalg.norm(a1.pos - a2.pos)
            dmin = min(dmin, d)
    return dmin

def dist_min_hbond(r1, r2):
    if is_non_polar(r1) or is_non_polar(r2):
        return (None, None, 1e12)

    dmin = 1e12
    a1m  = None
    a2m  = None

    for a1 in r1:
        for a2 in r2:
            #if not ((isH(a1) and isO(a2)) or (isH(a2) and isO(a1))):
            #    continue
            if not test_hbond(a1, a2):
                continue
            d = np.linalg.norm(a1.pos - a2.pos)
            if d < dmin:
                dmin = d
                a1m  = a1
                a2m  = a2
    return (a1m, a2m, dmin)

def test_hbond(a1, a2):

    aa = [a1, a2]

    r  = False
    i = 0
    while (not r) and (i < 2):
        ai = aa[i]
        aj = aa[(i+1)%2]
        r  = ((isH_polar(ai) and isO(aj))
              or (ai.res == 'H' and ai.atom.strip() == 'ND1' and isH_polar(aj)))

        i += 1

    return r

def is_non_polar(r):
    non_p = ['L', 'I', 'V', 'A', 'F', 'M']
    return r[0].res in non_p


def dist_chain(r, c):
    dmin = 1e12
    for ri in c:
        if ri:
            dmin = min(dmin, dist_CA(r, ri))
    return dmin

def isH(a):
    return a.atom[1] == 'H'

def isH_polar(a):
    return a.atom[1] == 'H' #and len(a.atom.strip()) > 2

def isO(a):
    return a.atom[1] == 'O'

def isBB(a):
    return a.atom in BB
