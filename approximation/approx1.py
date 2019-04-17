#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
---------------
Algorithme naïf
---------------
Ce fichier contient l'algorithme d'approximation "naïf" : On parcourt la matrice de dissimilarité, et on l'ajuste, si une valeur est trop grande par rapport à ses voisines, on la diminue au maximum de ses voisines.
"""
from dissimilarity import *

def approximate(d):
    """ Modifie "naïvement" la dissimilarité d pour la rendre Robinsonienne."""
    while not d.isCompatibleOrder():
        for x in range(0,d.n-1):
            for y in range(x+1,d.n):
                if x < y:
                    m = max(d[x,y], d[x+1,y],d[x,y-1])
                    if m > d[x,y]: d[x,y] = m

 

