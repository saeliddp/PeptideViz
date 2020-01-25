# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 19:33:46 2019

@author: saeli
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 11:18:39 2019

@author: saeli
"""

import pickle
from sequence_tree import *
import statistics

with open('delta_positional_lists.pickle', 'rb') as fr:
    das = pickle.load(fr)

aas = "ACDEFGHIKLMNPQRSTVWY"
means = {}
stdevs = {}
for i in range(12):
    means[i + 1] = {}
    stdevs[i + 1] = {}
    for aa2 in das[i + 1]:
        if len(das[i + 1][aa2]) > 0:
            means[i + 1][aa2] = statistics.mean(das[i + 1][aa2])
        else:
            means[i + 1][aa2] = None
            
        if len(das[i + 1][aa2]) > 1:
            stdevs[i + 1][aa2] = statistics.stdev(das[i + 1][aa2])
        else:
            stdevs[i + 1][aa2] = None


with open('means_pos.pickle', 'wb') as fw:
    pickle.dump(means, fw)

with open('stdevs_pos.pickle', 'wb') as fw:
    pickle.dump(stdevs, fw)
