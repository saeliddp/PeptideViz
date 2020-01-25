# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 21:18:39 2019

@author: saeli
"""

import pickle
from sequence_tree import *
import statistics
import pandas

aas = "ACDEFGHIKLMNPQRSTVWY"
data = []
for i in range(len(aas)):
    data.append([])
    for j in range(len(aas)):
        data[i].append([])
        for k in range(12):
            data[i][j].append([])

def getIndices(aa1, aa2):
    for i in range(len(aa1)):
        if aa1[i] != aa2[i]:
            return [aas.index(aa1[i]), aas.index(aa2[i]), i]

# this set represents the peptide pairs that have already been logged
# e.g. if A -> B has been logged already, it will not be recorded again
# for the tree that is being explored
deltas_explored = set()
def analyze_tree(node):
    if node.hasChildren():
        for child in node.children:
            delta = child.binding_affinity - node.binding_affinity
            inds = getIndices(node.aa_sequence, child.aa_sequence)
            key = node.aa_sequence + child.aa_sequence
            if key not in deltas_explored:
                data[inds[0]][inds[1]][inds[2]].append(delta)
                deltas_explored.add(key)
            analyze_tree(child)
            
def analyzeTree(tree):
    analyze_tree(tree.root)
    
with open('tree_dict.pickle', 'rb') as fr:
    tree_dict = pickle.load(fr)
    
for og_aa in tree_dict:
    for tree in tree_dict[og_aa]:
        analyzeTree(tree)
        deltas_explored = set()
means = []
stdevs = []
#100 is an arbitrarily high number

# for plotting purposes
meansxyzm = []
stdevsxyzm = []
for i in range(4):
    meansxyzm.append([])
    stdevsxyzm.append([])
    
for i in range(len(data)):
    means.append([])
    stdevs.append([])
    for j in range(len(data[i])):
        means[i].append([])
        stdevs[i].append([])
        for k in range(len(data[i][j])):
            deltas = data[i][j][k]
            if (len(deltas) > 0):
                mean = statistics.mean(deltas)
                means[i][j].append(mean)
                meansxyzm[0].append(i)
                meansxyzm[1].append(j)
                meansxyzm[2].append(k)
                meansxyzm[3].append(mean)
            else:
                means[i][j].append(0)
            if (len(deltas) > 1):
                stdev = statistics.stdev(deltas)
                stdevs[i][j].append(stdev)
                stdevsxyzm[0].append(i)
                stdevsxyzm[1].append(j)
                stdevsxyzm[2].append(k)
                stdevsxyzm[3].append(stdev)
            else:
                stdevs[i][j].append(-1)

with open('meansxyzm.pickle', 'wb') as fw:
    pickle.dump(meansxyzm, fw)
with open('stdevsxyzm.pickle', 'wb') as fw:
    pickle.dump(stdevsxyzm, fw)
    
