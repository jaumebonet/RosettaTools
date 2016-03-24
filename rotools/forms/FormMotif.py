# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-23 14:48:09
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-03-24 10:09:08
import os
from ..loops import Loops


class FormMotif(object):

    def __init__(self, desc):
        self.file     = desc["pdb_file"]
        self.chain    = desc["chain"] if "chain" in desc else "A"
        self.segments = desc["segments"]
        self.order    = {}
        for x in range(len(self.segments)):
            self.order[self.segments[x]["name"]] = x + 1

    def to_target_loops(self):
        l = Loops()
        for x in self.segments:
            l.add_loop(x["ini"], x["end"])
        return str(l)

    def to_tmpl_loops(self, form):
        l = Loops()
        for x in range(len(form.sslist)):
            if form.sslist[x].has_ref():
                print "---", form.sslist[x]
                l.add_loop(form.edges[x][0], form.edges[x][1])
        return str(l)

    def _set_order(self, form):
        refs  = filter(None, [x.ref for x in form.sslist])
        order = [self.order[x] for x in refs]
        return order

    def to_command(self, targetl, templtl, fastaf, ssefil, constrfl, ident, form, commfl):

        text = []
        text.append("-fold_from_loops:target:pose {0}".format(os.path.relpath(self.file, os.path.dirname(commfl))))
        text.append("-fold_from_loops:target:chain {0}".format(self.chain))
        text.append("-fold_from_loops:target:pdb_count")
        text.append("-fold_from_loops:target:loops {0}".format(os.path.relpath(targetl, os.path.dirname(commfl))))
        text.append("-fold_from_loops:target:order {0}".format(" ".join([str(x) for x in self._set_order(form)])))

        text.append("-in:file:s {0}".format(os.path.relpath(self.file, os.path.dirname(commfl))))
        text.append("-in:file:fasta {0}".format(os.path.relpath(fastaf, os.path.dirname(commfl))))
        text.append("-in:file:psipred_ss2 {0}".format(os.path.relpath(ssefil, os.path.dirname(commfl))))
        text.append("-loops:loop_file {0}".format(os.path.relpath(templtl, os.path.dirname(commfl))))
        text.append("-fold_from_loops:scaffold:nopose")
        text.append("-fold_from_loops:scaffold:cst_file {0}".format(os.path.relpath(constrfl, os.path.dirname(commfl))))

        text.append("-fold_from_loops:loop_mov_nterm 2")
        text.append("-fold_from_loops:loop_mov_cterm 2")

        text.append("-fold_from_loops:native_ca_cst")

        text.append("-in:ignore_unrecognized_res")
        text.append("-run:intermediate_structures")
        text.append("-out:nstruct 2")
        text.append("-out:file:silent_struct_type binary")
        text.append("-out:prefix {0}".format(ident))
        text.append("-out:file:silent {0}".format(ident))

        text.append("-out:overwrite")

        return "\n".join(text)
