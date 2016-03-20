# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-17 17:06:04
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-03-20 15:04:45
import os
import json
import math
import numpy as np
import scipy.spatial as scsp
import networkx as nx

from .SecondaryStructure import SecondaryStructure as SS
from ..constraints import ConstraintSet


class Form(object):
    _AA_DIST = 3.2
    _EDGE    = "CC"

    def __init__(self, sslist):
        self.sslist = sslist
        self.sequence  = ""
        self.structure = ""
        self.loops     = []

    def make_structure_sequence(self):
        self.structure += self._EDGE
        self.sequence  += "G" * len(self._EDGE)
        self.loops.append(len(self._EDGE))
        up = 1
        for x in range(len(self.sslist)):
            self.structure += (self.sslist[x].type * self.sslist[x].length)
            self.sequence  += self.sslist[x].sequence
            if x != len(self.sslist) - 1:
                dist = self.sslist[x].distance(self.sslist[x + 1], up)
                self.structure += ("C" * int(math.ceil(dist / self._AA_DIST)))
                self.sequence += ("G" * int(math.ceil(dist / self._AA_DIST)))
                self.loops.append(int(math.ceil(dist / self._AA_DIST)))
            up *= -1
        self.structure += self._EDGE
        self.sequence  += "G" * len(self._EDGE)
        self.loops.append(len(self._EDGE))

    def matches_desc(self):
        if not self._expected_edges():        return False
        if not self._expected_directions():   return False
        if not self._expected_intersection(): return False
        return True

    def to_fasta(self, name):
        return ">{0}\n{1.sequence}".format(name, self)

    def to_psipred_ss(self):
        text = []
        text.append("# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n")
        sse = self.structure
        for i in range(len(sse)):
            pH, pE, pC = 0, 0, 0
            if sse[i] == 'C': pC = 1
            else:
                edge = False
                if sse[i] != sse[i-1] or sse[i] != sse[i-2]: edge = True
                if sse[i] != sse[i+1] or sse[i] != sse[i+2]: edge = True
                if edge:
                    pC = 0.3
                    if sse[i] == 'E': pE = 0.7
                    else:             pH = 0.7
                else:
                    if sse[i] == 'E': pE = 1
                    else:             pH = 1

            line = "{0:>4} {1} {2}   {3:0.3f}  {4:0.3f}  {5:0.3f}".format(i + 1, 'A', sse[i], pC, pH, pE)
            text.append(line)
        return '\n'.join(text)

    def to_file_constraint(self):
        constraints = ConstraintSet()
        prelength   = 0
        positions   = []
        up          = 1
        for x in range(len(self.sslist)):
            prelength += self.loops[x]
            ini = (prelength + 1, self.sslist[x].get_xyz(up))
            mid = (ini[0] - 1 + self.sslist[x].length / 2, self.sslist[x].Mxyz)
            end = (ini[0] - 1 + self.sslist[x].length, self.sslist[x].get_xyz(up))
            prelength = end[0]
            up *= -1
            positions.append([ini, mid, end])
        for x in range(len(positions) - 1):
            for xx in range(len(positions[x])):
                xxx = positions[x][xx]
                for y in range(x + 1, len(positions)):
                    for yy in range(len(positions[y])):
                        yyy = positions[y][yy]
                        d = scsp.distance.cdist([xxx[1]], [yyy[1]], 'euclidean')[0][0]
                        constraints.add_constraint(xxx[0], yyy[0], d)
        return str(constraints)

    def _expected_intersection(self):
        dU = []
        dD = []
        up = 1
        for x in range(len(self.sslist) - 1):
            if up == 1:  dU.append((self.sslist[x], self.sslist[x + 1]))
            if up == -1: dD.append((self.sslist[x], self.sslist[x + 1]))
            up *= -1

        if len(dU) > 1:
            for x in dU:
                for y in dU:
                    if x != y and self._intersection(x, y, "up"):  return False
        if len(dD) > 1:
            for x in dD:
                for y in dD:
                    if x != y and self._intersection(x, y, "down"): return False
        return True

    def _intersection(self, s1, s2, key):
        left   = max(
                    min(s1[0].get_x(key), s1[1].get_x(key)),
                    min(s2[0].get_x(key), s2[1].get_x(key))
                 )
        right  = min(
                    max(s1[0].get_x(key), s1[1].get_x(key)),
                    max(s2[0].get_x(key), s2[1].get_x(key))
                 )
        top    = max(
                    min(s1[0].get_z(key), s1[1].get_z(key)),
                    min(s2[0].get_z(key), s2[1].get_z(key))
                 )
        bottom = min(
                    max(s1[0].get_z(key), s1[1].get_z(key)),
                    max(s2[0].get_z(key), s2[1].get_z(key))
                 )
        if bottom < top or right > left: return True
        return False

    def _expected_directions(self):
        d0 = [x.direct for x in self.sslist]
        d1 = set(d0[::2])
        d1.discard(0)
        if 1 in d1 and -1 in d1: return False
        d2 = set(d0[1::2])
        d2.discard(0)
        if 1 in d2 and -1 in d2: return False
        return len(d1.intersection(d2)) == 0

    def _expected_edges(self):
        edges = 0
        if self.sslist[0].edge == -1 or self.sslist[-1].edge == -1: return False
        for x in self.sslist:
            if x.edge == 1: edges += 1
        if edges > 2: raise AttributeError("More than two edge structures is not possible")
        if edges != 0 and self.sslist[0].edge + self.sslist[-1].edge != edges: return False
        return True

    def __str__(self):
        return "-".join([x.identifier for x in self.sslist])


class FormFabric(object):
    _DEPTH     = 9.1
    _DEPTH_VAR = 1.1
    _WIDTH     = {"H": 10.0, "E": 4.8, "C": 4.8, "X": 4.8}
    _WIDTH_VAR = {"H": 0.7,  "E": 0.2, "C": 0.4, "X": 0.3}

    def __init__(self):
        self._id     = None
        self._desc   = None
        self._layers = []

        self.forms     = []
        self.Estandard = 8

    def build(self, identifier, filename):
        self._id = identifier
        with open(filename) as fd:
            self._process(json.loads(''.join([l.strip() for l in fd])))

    def dump(self, outdir = os.path.join(os.getcwd(), 'forms')):
        outdir = os.path.join(outdir, self._id)
        if not os.path.isdir(outdir): os.makedirs(outdir)
        for x in range(len(self.forms)):
            ident = "{0}_{1:06d}".format(self._id, x + 1)
            finaldir = os.path.join(outdir, ident)
            if not os.path.isdir(finaldir): os.mkdir(finaldir)
            fasta   = os.path.join(finaldir, "fasta.fa")
            psipred = os.path.join(finaldir, "psipred.ss2")
            constrs = os.path.join(finaldir, "constraints.cs")
            with open(fasta, "w") as fd:   fd.write(self.forms[x].to_fasta(self._id))
            with open(psipred, "w") as fd: fd.write(self.forms[x].to_psipred_ss())
            with open(constrs, "w") as fd: fd.write(self.forms[x].to_file_constraint())

    def print_structures(self):
        for x in self._layers:
            for y in x:
                print y

    def reset(self):
        self._id     = None
        self._desc   = None
        self._layers = None
        self.forms   = []

    def _process(self, description):
        self._desc = description
        maxL = {"H": 0, "E": 0, "C": 0, "X": 0}
        for x in range(len(self._desc["layers"])):
            self._layers.append([])
            for y in range(len(self._desc["layers"][x])):
                ss = SS(self._desc["layers"][x][y], x+1, y+1)
                self._layers[-1].append(ss)
                if ss.length > maxL[ss.type]: maxL[ss.type] = ss.length

        self._apply_lengths(maxL)
        self._place_xz()
        self._create_forms()

    def _create_forms( self ):
        i = 1
        G = self._create_graph()
        path_length = len( G.nodes() )-1
        self.forms = []
        for node in G.nodes():
            for path in self._find_paths(G, node, path_length):
                f = Form(path)
                if not f.matches_desc(): continue
                f.make_structure_sequence()
                self.forms.append(f)
                i += 1

    def _create_graph( self ):
        G = nx.Graph()
        for x in self._layers:
            for sse1 in x:
                for sse2 in x:
                    if sse1 < sse2:
                        G.add_edge( sse1 , sse2, object=SS )
        for lyr1 in range(len(self._layers)):
            for lyr2 in range(len(self._layers)):
                if abs(lyr1 - lyr2) == 1:  # Only consecutive layers
                    for sse1 in self._layers[lyr1]:
                        for sse2 in self._layers[lyr2]:
                            G.add_edge( sse1 , sse2, object=SS )
        return G

    def _find_paths( self, G, u, n ):
        if n == 0: return [[u]]
        paths = [[u]+path for neighbor in G.neighbors(u) for path in self._find_paths(G, neighbor, n-1) if u not in path]
        return paths

    def _place_xz(self):
        r = np.random.random_sample()

        for x in range(len(self._layers)):
            dvar = float(self._DEPTH_VAR * r) - (self._DEPTH_VAR / 2.0)
            z    = float((self._DEPTH + dvar) * x)
            last = 0
            for y in range(len(self._layers[x])):
                self._layers[x][y].set_z(z)
                tp   = self._layers[x][y].type
                wvar = (self._WIDTH_VAR[tp] * r) - (self._WIDTH_VAR[tp] / 2.0)
                xp   = last + self._WIDTH[tp] + wvar
                last = xp
                self._layers[x][y].set_x(xp)

    def _apply_lengths(self, lengths):
        if lengths["H"] == 0 and lengths["E"] == 0 and lengths["C"] == 0:
            lengths["E"] = self.Estandard
        if lengths["H"] != 0:
            if lengths["E"] == 0: lengths["E"] = int(lengths["H"]/2.0)
            if lengths["C"] == 0: lengths["C"] = int(lengths["H"]/2.0)
        if lengths["E"] != 0:
            if lengths["H"] == 0: lengths["H"] = int(lengths["E"] * 2.0)
            if lengths["C"] == 0: lengths["C"] = int(lengths["E"] * 2.0)
        if lengths["C"] != 0:
            if lengths["H"] == 0: lengths["H"] = int(lengths["C"] * 2.0)
            if lengths["E"] == 0: lengths["E"] = int(lengths["C"] * 2.0)
        for x in self._layers:
            for y in x:
                if y.length == 0: y.length = lengths[y.type]

