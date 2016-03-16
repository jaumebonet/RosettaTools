# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-15 14:47:41
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-03-16 18:49:08
"""
According to RosettaWiki:

Column -- Meaning
1      -- blank
2-5    -- PDB code for the fragment origin
7      -- chain ID for the origin PDB
9-13   -- PDB residue number for the origin PDB
15     -- amino acid identity in the origin PDB
17     -- secondary structure for the origin PDB (Helix, Loop, Extended/beta)
19-27  -- phi
28-36  -- psi
37-45  -- omega
46-54  -- C-alpha x coordinate for origin PDB (optional)
55-63  -- C-alpha y coordinate for origin PDB (optional)
65-73  -- C-alpha z coordinate for origin PDB (optional)
74-79  -- unknown (unused)
80-85  -- unknown (unused)
86     -- Literal "P" (unused)
87-89  -- fragment position number, pose numbered (unused)
91     -- Literal "F"(unused)
92-94  -- fragment number (unused)

0         1         2         3         4         5         6         7         8         9
 2od4 A     5 G L   59.826 -156.477 -179.462   -0.936    4.890     8.153 3     0.000 P  1 F  1
"""


class FragmentPosition( object ):
    """Contains one position of Rosetta Fragment data"""
    def __init__(self):
        self.pdb     = None
        self.chain   = None
        self.resnum  = None
        self.aatype  = None
        self.secstr  = None
        self.phi     = None
        self.psi     = None
        self.omega   = None
        self.coord   = [None, None, None]
        self.ukn1    = None
        self.ukn2    = None
        self.p       = "P"
        self.posenum = 0
        self.f       = "F"
        self.fragnum = 0

    @staticmethod
    def parse(line):
        line = line.rstrip()
        frag = FragmentPosition()
        frag.pdb     = line[1:5].strip()
        frag.chain   = line[6].strip()
        frag.resnum  = int(line[8:13].strip())
        frag.aatype  = line[14].strip()
        frag.secstr  = line[16].strip()
        frag.phi     = float(line[17:26].strip())
        frag.psi     = float(line[26:35].strip())
        frag.omega   = float(line[35:44].strip())
        if len(line) > 44:
            if len(line[44:53].strip()) > 0:
                frag.coord   = [float(line[44:53].strip()),
                                float(line[53:62].strip()),
                                float(line[62:72].strip())]
            frag.ukn1    = line[72:78]
            frag.ukn2    = line[79:85]
            frag.p       = line[85].strip()
            frag.posenum = int(line[86:89].strip())
            frag.f       = line[90].strip()
            frag.fragnum = int(line[91:].strip())

        return frag

    def __str__(self):
        text  = " {0.pdb} {0.chain} {0.resnum:>5} ".format(self)
        text += "{0.aatype} {0.secstr} {0.phi:>8.3f} ".format(self)
        text += "{0.psi:>8.3f} {0.omega:>8.3f}".format(self)
        if self.coord[0] is not None:
            text += "{0.coord[0]:>9.3f}{0.coord[1]:>9.3f}".format(self)
            text += " {0.coord[2]:>9.3f}".format(self)
            text += "{0.ukn1} {0.ukn2}".format(self)
        else: text += "{0:>41}".format("")
        text += "{0.p}{0.posenum:>3} ".format(self)
        text += "{0.f}{0.fragnum:>3}".format(self)

        return text


class Fragment(object):
    """Containt a Rosetta Fragment data"""
    def __init__(self):
        self.positions = []

    def append(self, position):
        self.positions.append(position)

    def __iter__(self):
        return self.positions.__iter__()

    def __getitem__(self, key):
        return self.positions[key]

    def __len__(self):
        return len(self.positions)

    def __str__(self):
        return "\n".join([str(x) for x in self.positions])
