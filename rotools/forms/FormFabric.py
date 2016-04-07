# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-17 17:06:04
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-04-07 07:21:56
import os
import json
import numpy as np
import networkx as nx

from .SecondaryStructure import SecondaryStructure as SS
from .Form import Form
from .FormMotif import FormMotif


class FormFabric(object):
    # _DEPTH     = 7.1
    _DEPTH     = 9.1
    _DEPTH_VAR = 1.1
    # _WIDTH     = {"H": 9.0, "E": 4.3, "C": 4.3, "X": 4.3}
    _WIDTH     = {"H": 10.0, "E": 4.8, "C": 4.8, "X": 4.8}
    _WIDTH_VAR = {"H": 0.7,  "E": 0.2, "C": 0.4, "X": 0.3}

    def __init__(self):
        self._id     = None
        self._desc   = None
        self._layers = []
        self._motif  = None

        self.forms     = []
        self.Estandard = 8

    def build(self, identifier, filename):
        self._id = identifier
        data     = {}
        with open(filename) as fd:
            data = json.loads(''.join([l.strip() for l in fd]))
        if "motif" in data:
            p1 = os.path.dirname(os.path.abspath(filename))
            p2 = os.path.abspath(data["motif"]["pdb_file"])
            data["motif"]["pdb_file"] = os.path.relpath(p2, p1)
        self._process(data)

    def dump(self, outdir = os.path.join(os.getcwd(), 'forms')):
        outdir = os.path.join(outdir, self._id)
        if not os.path.isdir(outdir): os.makedirs(outdir)
        for x in range(len(self.forms)):
            ident = "{0}_{1:06d}".format(self._id, x + 1)
            finaldir = os.path.join(outdir, ident)
            if not os.path.isdir(finaldir): os.mkdir(finaldir)
            identf  = os.path.join(finaldir, "info.md")
            fasta   = os.path.join(finaldir, "fasta.fa")
            psipred = os.path.join(finaldir, "psipred.ss2")
            constrs = os.path.join(finaldir, "constraints.cst")
            targetl = os.path.join(finaldir, "target.loop")
            tmplatl = os.path.join(finaldir, "template.loop")
            comndfl = os.path.join(finaldir, "run.command")
            with open(identf, "w")  as fd: fd.write(str(self.forms[x]))
            with open(fasta, "w")   as fd: fd.write(self.forms[x].to_fasta(self._id))
            with open(psipred, "w") as fd: fd.write(self.forms[x].to_psipred_ss())
            with open(constrs, "w") as fd: fd.write(self.forms[x].to_file_constraint())
            with open(targetl, "w") as fd: fd.write(self._motif.to_target_loops())
            with open(tmplatl, "w") as fd: fd.write(self._motif.to_tmpl_loops(self.forms[x]))
            if not os.path.isfile(comndfl):
                with open(comndfl, "w") as fd: fd.write(self._motif.to_command(targetl, tmplatl, fasta, psipred, constrs, ident, self.forms[x], comndfl))

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
        self._desc  = description
        self._motif = FormMotif(description["motif"])
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

