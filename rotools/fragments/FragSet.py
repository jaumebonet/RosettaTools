# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-15 15:54:35
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-03-17 15:19:15
import os
import re
from collections import Counter

from .Fragment import Fragment, FragmentPosition


class FragSet(object):
    """Set of Fragments for multiple positions"""
    def __init__(self):
        self.title     = "Fragments"
        self.fragments = {}

    @staticmethod
    def parse(fragfile):
        fset    = FragSet()
        ppos    = FragmentPosition()
        re_pnum = re.compile("position:\s+(\d+)")
        posenum = 0
        fragnum = 0
        with open(fragfile) as fd:
            for line in fd:
                g = re.match(re_pnum, line)
                if g:
                    posenum = int(g.group(1))
                    fragnum = 0
                    continue
                if len(line.strip()) == 0:
                    fragnum += 1
                    continue

                fpos = FragmentPosition.parse(line)
                if fpos.posenum == 0: fpos.posenum = posenum
                if fpos.fragnum == 0: fpos.fragnum = fragnum
                fset.fragments.setdefault(fpos.posenum, [])
                if fpos.fragnum != ppos.fragnum:
                    fset.fragments[fpos.posenum].append(Fragment())
                fset.fragments[fpos.posenum][-1].append(fpos)
                ppos = fpos
        fset.title = os.path.split(fragfile)[-1]
        return fset

    def candidates(self):
        return len(self.fragments[1])

    def min_position(self):
        return min(self.fragments.keys())

    def max_position(self):
        return max(self.fragments.keys())

    def plot(self, fileprefix = None, format = "svg", to_file = True, show = True):
        try:
            import seaborn as sns
            import matplotlib.pyplot as plt
        except:
            raise ImportError("Plotting fragments requires seaborn and matplotlib.")

        phi = [[] for i in range(self.max_position() + 1)]
        psi = [[] for i in range(self.max_position() + 1)]
        sse = [[] for i in range(self.max_position() + 1)]
        seq = [[] for i in range(self.max_position() + 1)]
        for x in range(self.min_position(), self.max_position() + 1):
            for y in range(len(self[x])):
                for z in range(len(self[x][y])):
                    if (x + z) >= len(phi):
                        phi.append([])
                        psi.append([])
                        sse.append([])
                        seq.append([])
                    phi[x + z].append(self[x][y][z].phi)
                    psi[x + z].append(self[x][y][z].psi)
                    sse[x + z].append(self[x][y][z].secstr)
                    seq[x + z].append(self[x][y][z].aatype)
        for x in range(len(sse)):
            if len(sse[x]) == 0: sse[x] = ""
            else:
                data   = Counter(sse[x])
                sse[x] = data.most_common()[0][0]
            if len(seq[x]) == 0: seq[x] = ""
            else:
                data   = Counter(seq[x])
                seq[x] = data.most_common()[0][0]

        fig = plt.figure()
        fig.suptitle(self.title)
        fig.set_tight_layout(True)
        fig.set_figheight(8)
        fig.set_figwidth(1  + 0.2 * len(sse))

        ax1 = fig.add_subplot(2, 1, 1)
        ax1.set_xlim([-360, 360])
        ax1.set_ylabel('phi')
        ax1.set_ylim([-360, 360])
        axs1 = sns.boxplot(data=phi, ax=ax1)
        ax1.set_xticklabels(sse)

        boxes = axs1.artists
        for i, box in enumerate(boxes):
            if 'H' in sse[i+1]:   box.set_facecolor('b')
            elif 'E' in sse[i+1]: box.set_facecolor('r')
            else:                 box.set_facecolor('#98FB98')

        ax2 = fig.add_subplot(2, 1, 2)
        ax2.xaxis.set_label_position('top')
        ax2.set_ylabel('psi')
        ax2.set_ylim([-360, 360])
        axs2 = sns.boxplot(data=psi, ax=ax2)
        ax2.xaxis.set_ticks_position('top')
        ax2.set_xticklabels(seq)

        boxes = axs2.artists
        for i, box in enumerate(boxes):
            if 'H' in sse[i+1]:   box.set_facecolor('b')
            elif 'E' in sse[i+1]: box.set_facecolor('r')
            else:                 box.set_facecolor('#98FB98')

        if to_file:
            filename = ""
            if fileprefix is None: filename = self.title + '.' + format
            else:                  filename = fileprefix + '.' + format
            fig.savefig(filename, dpi=300)

        if show: plt.show()

    def __getitem__(self, key):
        return self.fragments[key]

    def __iter__(self):
        for x in range(self.min_position(), self.max_position() + 1):
            yield self.fragments[x]

    def __len__(self):
        return len(self.fragments)

    def __str__(self):
        text = []
        for x in range(self.min_position(), self.max_position() + 1):
            text.append(" position:{0:>13} neighbors:{1:>13}\n".format(x, len(self[x])))
            for y in self[x]:
                text.append(str(y))
                text.append("")
        text.pop()
        return "\n".join(text)
