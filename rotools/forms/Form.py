# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-23 15:01:29
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-03-23 16:50:02
import math
import numpy as np
import scipy.spatial as scsp

from ..constraints import ConstraintSet


class Form(object):
    _AA_DIST = 3.2
    _EDGE    = "CC"

    def __init__(self, sslist):
        self.sslist    = sslist
        self.sequence  = ""
        self.structure = ""
        self.loops     = []
        self.edges     = []

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
            self.edges.append((ini[0], end[0]))

        positions = self._expand_points(positions)  # EXPAND POINT... to.... manual...

        for x in range(len(positions) - 1):
            for xx in range(len(positions[x])):
                xxx = positions[x][xx]
                for y in range(x + 1, len(positions)):
                    for yy in range(len(positions[y])):
                        yyy = positions[y][yy]
                        d = scsp.distance.cdist([xxx[1]], [yyy[1]], 'euclidean')[0][0]
                        constraints.add_constraint(xxx[0], yyy[0], d)
        return str(constraints)

    def _expand_points(self, positions):
        new_p = []
        for x in range(len(positions)):
            new_pos = []
            for y in range(len(positions[x])):
                new_pos.append(positions[x][y])
                if y < len(positions[x]) - 1:
                    p = positions[x][y][0] + (positions[x][y + 1][0] - positions[x][y][0]) / 2
                    xx = (positions[x][y + 1][1][0] + positions[x][y][1][0]) / 2.0
                    yy = (positions[x][y + 1][1][1] + positions[x][y][1][1]) / 2.0
                    zz = (positions[x][y + 1][1][2] + positions[x][y][1][2]) / 2.0
                    new_pos.append((p, np.array([xx, yy, zz])))
            new_p.append(new_pos)
        return new_p

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
