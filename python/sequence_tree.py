# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 16:19:24 2019

@author: saeli
"""
import pandas
import numpy

# a sequence node represents a specific peptide in the set
class SequenceNode(object):
    #@param row contains AA_seq, binding_affinity, 
    def __init__(self, row=None, aa=None, ba=None):
        if row is not None:
            self.aa_sequence = row.AA_seq
            self.binding_affinity = row.binding_affinity
        else:
            self.aa_sequence = aa
            self.binding_affinity = ba
            
        self.children = []
    
    def __eq__(self, other):
        return self.aa_sequence == other.aa_sequence
    
    def addChild(self, row):
        self.children.append(SequenceNode(row=row))
    
    def hasChildren(self):
        return len(self.children) > 0
    
    def offByOne(self, aa):
        num = 0
        for i in range(len(aa)):
            if aa[i] == self.aa_sequence[i]:
                num += 1
        return num == 11

    # adds all 'child nodes' (peptides off by one) of this node and returns
    # a new dataframe excluding rows of the child nodes and this node
    def addAllChildren(self, df):
        for row in df.itertuples():
            if self.offByOne(row.AA_seq):
                self.addChild(row)
                
# represents a tree of 'off by one' peptides
class SequenceTree(object):
    def __init__(self, node):
        self.root = node
        self.roots_built = [node]
    
    # builds a tree of depth 4
    def buildTree(self, df):
        self.build_tree(df, 4, self.root)
    
    def build_tree(self, df, layers, rt):
        if layers != 0:
            rt.addAllChildren(df)
            for child in rt.children:
                if child not in self.roots_built:
                    self.roots_built.append(child)
                    self.build_tree(df, layers - 1, child)
                    self.roots_built.pop()
                
    # prints tree to console
    def printTree(self):
        self.print_tree(self.root, 0)
    
    def print_tree(self, rt, num_ind):
        st = ""
        for i in range(num_ind):
            st = st + "\t"
        print(st + rt.aa_sequence)
        if len(rt.children) > 0:
            for child in rt.children:
                self.print_tree(child, num_ind + 1)
    
    # prints tree to specified file
    def printToFile(self, filename):
        file = open(filename, 'w')
        self.print_to_file(self.root, file, 0)
        file.close()
        
    def print_to_file(self, rt, file, num_ind):
        st = ""
        for i in range(num_ind):
            st = st + "\t"
        file.write(st + rt.aa_sequence + '\n')
        if len(rt.children) > 0:
            for child in rt.children:
                self.print_to_file(child, file, num_ind + 1)