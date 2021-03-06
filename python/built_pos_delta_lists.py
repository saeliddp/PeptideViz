# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 19:28:12 2019

@author: saeli
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 09:56:19 2019

@author: saeli
"""
import pickle
from sequence_tree import *
# create a 12X20 matrix of lists. each list will store every occurence of the
# change in binding affinity to a specific amino acid in a specific position

delta_aa_lists = {}
aas = "ACDEFGHIKLMNPQRSTVWY"

for i in range(12):
    delta_aa_lists[i + 1] = {}
    for aa2 in list(aas):
        delta_aa_lists[i + 1][aa2] = []

# returns what keys represent aa1 -> aa2 in delta_aa_lists
def getKeys(aa1, aa2):
    for i in range(len(aa1)):
        if aa1[i] != aa2[i]:
            return [i + 1, aa2[i]]
        
deltas_explored = set()      
def analyze_tree(node):
    if node.hasChildren():
        for child in node.children:
            delta = child.binding_affinity - node.binding_affinity
            keys = getKeys(node.aa_sequence, child.aa_sequence)
            key = node.aa_sequence + child.aa_sequence
            if key not in deltas_explored:
                delta_aa_lists[keys[0]][keys[1]].append(delta)
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

with open('delta_positional_lists.pickle', 'wb') as fw:
    pickle.dump(delta_aa_lists, fw)