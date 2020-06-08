# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 17:37:26 2017
@author: Ivan Alisson Cavalcante Nunes de Lima
"""

class node:
    """
        Class tha stores each string in a node, with a character, a index before alignment and a index after the alignment
    """
    amino = None
    degree = None
    aminoId = None
    thisNode = None     #pdbPos = None
    targetNode = None   #seqPos = None
    ClusteringCoef = None
    Betweennessweighted = None
    
    def __init__(self, amino, degree, aminoId, thisNode, targetNode, clustCoef, betWeigh):
        self.amino = amino
        self.degree = degree
        self.aminoId = aminoId
        self.thisNode = thisNode
        self.targetNode = targetNode
        self.ClusteringCoef = clustCoef
        self.Betweennessweighted = betWeigh

    def getAll(self):
        return self.amino, self.degree, self.aminoId, self.thisNode, self.targetNode, self.ClusteringCoef, self.Betweennessweighted
    
    def setAll(self, amino, degree, aminoId, thisNode, targetNode, clustCoef, betWeigh):
        self.amino = amino
        self.degree = degree
        self.aminoId = aminoId
        self.thisNode = thisNode
        self.targetNode = targetNode
        self.ClusteringCoef = clustCoef
        self.Betweennessweighted = betWeigh