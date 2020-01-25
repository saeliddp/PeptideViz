# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 11:18:39 2019

@author: saeli
"""

import pickle
from sequence_tree import *
import statistics

with open('delta_aa_lists.pickle', 'rb') as fr:
    das = pickle.load(fr)

aas = "ACDEFGHIKLMNPQRSTVWY"
means = {}
stdevs = {}
for aa in das:
    means[aa] = {}
    stdevs[aa] = {}
    for aa2 in das[aa]:
        if len(das[aa][aa2]) > 0:
            means[aa][aa2] = statistics.mean(das[aa][aa2])
        else:
            means[aa][aa2] = None
            
        if len(das[aa][aa2]) > 1:
            stdevs[aa][aa2] = statistics.stdev(das[aa][aa2])
        else:
            stdevs[aa][aa2] = None


with open('means.pickle', 'wb') as fw:
    pickle.dump(means, fw)

with open('stdevs.pickle', 'wb') as fw:
    pickle.dump(stdevs, fw)
