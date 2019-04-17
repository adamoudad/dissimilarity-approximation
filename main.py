#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Programme principal
===================
Il s'agit du programme principal, servant à lancer une simulation.

Une fonction d'affichage de titre pour mieux organiser les différentes parties des tests est fournie.
"""
from dissimilarity import *
from sort.sortElements import *
from approximation  import approx1,approx2
from generation.generateRandom import *
import numpy as np
import time

def printTitle(s, level=1):
    """ Affiche un titre mis en forme dans la console """
    if level == 0:
        max_length = 30
        sep = " "
        line = sep*((max_length - len(s)//2 - 2))
        if not len(s) % 2: s += " "
        print("\n\n",line + " ",s.upper()," " + line, "\n\n")
    else:
        if level == 1:
            max_length = 30
            sep = "-"
        elif level == 2:
            max_length = 10
            sep = "*"
        line = sep*((max_length - len(s)//2 - 2))
        if not len(s) % 2: s += " "
        print(line + " ",s," " + line)

printTitle("Approximation par dissimilarités Robinsoniennes",0)

printTitle("Dissimlarité Utilisée")
d = generateRobinsonian(20)
m0 = np.array(d)                # Matrice de départ, robinsonienne
# print(d)
randomize(d)
m = np.array(d)
distance = np.linalg.norm(m-m0)
print("Mesure de la perturbation introduite : |d0 - d| = ", np.linalg.norm(m-m0))
print(d)

t0 = time.clock()

printTitle("Tri")
sortElements(d)
print(d)
sortElements(d)
sortElements(d)
sortElements(d)
sortElements(d)
printTitle("Résultats", 2)
print(d)
printTitle("Performances", 2)
print("Temps : ",time.clock() - t0)
print("Ordre compatible ?", d.isCompatibleOrder())

old_order = list(d.order)
dR = Dissimilarity(list(d.matrix))
t0 = time.clock()

printTitle("Approximation 1 : méthode naïve")
approx1.approximate(dR)
printTitle("Résultats", 2)
# print(dR)
printTitle("Performances", 2)
print("Temps : ",time.clock() - t0)
print("Ordre compatible ?", dR.isCompatibleOrder())
dR.update(old_order)
mR = np.array(dR)
print("Distance L² : ", np.linalg.norm(m - mR))

dR = Dissimilarity(list(d.matrix))
t0 = time.clock()

printTitle("Approximation 2")
dR = approx2.approximate(dR)
printTitle("Résultats", 2)
# print(dR)
printTitle("Performances", 2)
print("Temps : ",time.clock() - t0)
print("Ordre compatible ?", dR.isCompatibleOrder())
dR.update(old_order)
mR = np.array(dR)
print("Distance L² : ", np.linalg.norm(m - mR))


