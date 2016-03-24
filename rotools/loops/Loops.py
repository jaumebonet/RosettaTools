# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-23 16:33:56
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-03-24 10:56:59


class Loops(object):
    def __init__(self):
        self.loops = []

    def add_loop(self, ini, end):
        self.loops.append((ini, end))

    def __str__(self):
        text = []
        text.append("#LOOP start end cutpoint skip-rate extend")
        for l in self.loops:
            text.append("LOOP {0[0]} {0[1]} 0 0.0 1".format(l))
        return "\n".join(text) + "\n"
