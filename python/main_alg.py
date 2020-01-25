# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 12:52:35 2019

@author: saeli
"""
import transfer # will build the merged csv
print("Merged CSV built")
import pandas
import numpy
from sequence_tree import *
import pickle
initial_list = pandas.read_csv("merged.csv")

# gets a list of the original 96 peptides
og = open('og_peptides.txt', 'r').readlines()

# maps original peptides to all of their 7-similar peptides (in the form of
# SequenceNodes)
og_dict = {}

# maps original peptides to all of the trees built from their 7-similar peptides
tree_dict = {}


# remove trailing newlines and sets dictionary keys
for pep in og:
    og_dict[pep[:12]] = []
    tree_dict[pep[:12]] = []

# returns true if pep1 is similar to pep2 by at least 7/12 amino acids
# with 'similar' meaning the same amino acid in the same position
def isSimilarByAtLeast(num_aa, pep1, pep2):
    sim = 0
    for i in range(len(pep1)):
        if pep1[i] == pep2[i]:
            sim += 1
    return sim >= num_aa

# for each of the original 96, get all similar by 7
numparsed = 1
for og_seq in og_dict:
    for row in initial_list.itertuples():
        if isSimilarByAtLeast(7, row.AA_seq, og_seq):
            og_dict[og_seq].append(SequenceNode(row=row))
    print(og_seq + " " + str(len(og_dict[og_seq])) + " " + str(numparsed))
    numparsed += 1
    
    
# build trees for each of the nodes in og_dict and place those trees into tree_dict
for aa in og_dict:
    #print('original: ' + aa)
    #print('num 7 similar: ' +str(len(og_dict[aa])))
    for node in og_dict[aa]:
        #print(node.aa_sequence + "########")
        tree = SequenceTree(node)
        tree.buildTree(initial_list)
        tree.printToFile('text_trees/' + aa + '_' + node.aa_sequence + '.txt')
        tree_dict[aa].append(tree)

with open('tree_dict.pickle', 'wb') as fw:
    pickle.dump(tree_dict, fw)
