# -*- coding: utf-8 -*-
# @Author: Jaume Bonet
# @Date:   2016-03-24 13:17:30
# @Last Modified by:   Jaume Bonet
# @Last Modified time: 2016-03-24 14:18:05
import os
import sys
from collections import Counter
from SBI.structure import PDB, InnerContacts

# Including rotools without adding them to the PYTHONPATH
scrdir  = os.path.dirname(os.path.realpath(__file__))
rotools = os.path.join(scrdir, "../..")
sys.path.append(rotools)
from rotools.constraints import ConstraintSet


def similar_to(d, mind, ca, cb, geo, bck):
    if mind >= float(d) - 3 and mind <= float(d) + 3: return "MIN"
    if ca >= float(d) - 3 and ca <= float(d) + 3: return "CA"
    if cb >= float(d) - 3 and cb <= float(d) + 3: return "CB"
    if geo >= float(d) - 3 and geo <= float(d) + 3: return "GEO"
    if bck >= float(d) - 3 and bck <= float(d) + 3: return "BCK"
    return "NONE"


cs = ConstraintSet.parse(sys.argv[1])
for f in os.listdir(os.getcwd()):
    if f.endswith('pdb'):
        data = []
        inc = InnerContacts(PDB(f), AA=True, AA_distance=35, HT=False)
        for cont in inc.AAcontacts[0].contacts:
            if cs.has_contact(cont.aminoacid1.number, cont.aminoacid2.number):
                d = cs.get_contact(cont.aminoacid1.number, cont.aminoacid2.number).value
                s = similar_to(d, cont.min_distance, cont.ca_distance, cont.cb_distance, cont.geometric_distance, cont.backbone_distance)
                data.append(s)
                print f, cont.aminoacid1.number, cont.aminoacid2.number, d, cont.min_distance, cont.ca_distance, cont.cb_distance, cont.geometric_distance, cont.backbone_distance, s
        cc = Counter(data)
        print "SUMMARY", f, cc, len(data)
