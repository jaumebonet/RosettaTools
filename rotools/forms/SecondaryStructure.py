# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-18 11:53:53
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-03-23 16:47:44

from collections import OrderedDict
import json
from random import random
from bisect import bisect
import numpy as np
import scipy.spatial as scsp

# CHOP780201 alpha-helix propensity AAindex (Chou-Fasman, 1978b)
helix_aa = [("A", 1.42), ("L", 1.21), ("R", 0.98), ("K", 1.16), ("N", 0.67),
            ("M", 1.45), ("D", 1.01), ("F", 1.13), ("C", 0.70), ("P", 0.57),
            ("Q", 1.11), ("S", 0.77), ("E", 1.51), ("T", 0.83), ("G", 0.57),
            ("W", 1.08), ("H", 1.00), ("Y", 0.69), ("I", 1.08), ("V", 1.06)]
# CHOP780202 beta-sheet propensity AAindex (Chou-Fasman, 1978b)
betas_aa = [("A", 0.83), ("L", 1.30), ("R", 0.93), ("K", 0.74), ("N", 0.89),
            ("M", 1.05), ("D", 0.54), ("F", 1.38), ("C", 1.19), ("P", 0.55),
            ("Q", 1.10), ("S", 0.75), ("E", 0.37), ("T", 1.19), ("G", 0.75),
            ("W", 1.37), ("H", 0.87), ("Y", 1.47), ("I", 1.60), ("V", 1.70)]
rand_seq = {"H": helix_aa, "E": betas_aa, "C": [("G", 1.00)]}


def weighted_choice(choices):
    values, weights = zip(*choices)
    total = 0
    cum_weights = []
    for w in weights:
        total += w
        cum_weights.append(total)
    x = random() * total
    i = bisect(cum_weights, x)
    return values[i]


class SecondaryStructure(object):
    _PER_RES = {"H": 1.5, "E": 3.5, "C": 3.5, "X": 0}

    def __init__(self, desc, layer, lcount):
        self.desc    = desc
        self.layer   = int(layer)
        self.lcount  = int(lcount)
        self.type    = str(desc["type"])
        self._length = 0
        self.direct  = int(desc["direction"]) if "direction" in desc else 0
        self.edge    = int(desc["edge"])      if "edge"      in desc else 0
        self.Uxyz    = np.array([0.0, 0.0, 0.0])
        self.Mxyz    = np.array([0.0, 0.0, 0.0])
        self.Dxyz    = np.array([0.0, 0.0, 0.0])
        self.ref     = desc["ref"]            if "ref"       in desc else None
        self._seq    = int(desc["sequence"])  if "sequence"  in desc else None

        self.length  = desc["length"] if "length" in desc else 0

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, value):
        self._length = value
        self.Uxyz[1] = float(float(self._length) * self._PER_RES[self.type])/2.0
        self.Dxyz[1] = -self.Uxyz[1]

    @property
    def identifier(self):
        return chr(self.layer + 65 - 1) + str(self.lcount) + self.type

    @property
    def sequence(self):
        if self._seq is None:
            self._seq = ""
            for x in range(self.length):
                self._seq += weighted_choice(rand_seq[self.type])
        return self._seq

    def has_ref(self):
        return self.ref is not None

    def get_xyz(self, key):
        if key == "up" or key == 1:     return self.Uxyz
        if key == "middle" or key == 0: return self.Mxyz
        if key == "down" or key == -1:  return self.Dxyz

    def get_x(self, key):
        return self.get_xyz(key)[0]

    def get_y(self, key):
        return self.get_xyz(key)[1]

    def get_z(self, key):
        return self.get_xyz(key)[2]

    def set_x(self, value):
        self.Uxyz[0] = value
        self.Mxyz[0] = value
        self.Dxyz[0] = value

    def set_z(self, value):
        self.Uxyz[2] = value
        self.Mxyz[2] = value
        self.Dxyz[2] = value

    def Udistance(self, other):
        return scsp.distance.cdist([self.Uxyz], [other.Uxyz], 'euclidean')[0][0]

    def Mdistance(self, other):
        return scsp.distance.cdist([self.Mxyz], [other.Mxyz], 'euclidean')[0][0]

    def Ddistance(self, other):
        return scsp.distance.cdist([self.Dxyz], [other.Dxyz], 'euclidean')[0][0]

    def distance(self, other, key):
        if key == "up" or key == 1:     return self.Udistance(other)
        if key == "middle" or key == 0: return self.Mdistance(other)
        if key == "down" or key == -1:  return self.Ddistance(other)

    def to_dict(self, full = False):
        info = OrderedDict()
        if full:
            info["layer"] = self.layer
            info["count"] = self.lcount
        info["type"]      = self.type
        info["length"]    = self.length
        info["direction"] = self.direct
        info["edge"]      = self.edge
        if self.ref is not None:
            info["ref"]   = self.ref
        if full:
            info["Uxyz"]  = [round(x, 2) for x in self.Uxyz]
            info["Mxyz"]  = [round(x, 2) for x in self.Mxyz]
            info["Dxyz"]  = [round(x, 2) for x in self.Dxyz]
        return info

    def __eq__( self, other ):
        return self.__hash__() == other.__hash__()

    def __ne__( self, other ):
        return not self.__eq__( other )

    def __lt__( self, other ):
        return self.__hash__() < other.__hash__()

    def __hash__( self ):
        return self.layer * 100 + self.lcount

    def __str__( self ):
        return json.dumps(self.to_dict(True))

    def __repr__( self ):
        return self.__str__()
