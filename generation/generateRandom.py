#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Pour réaliser des tests de performances des algorithmes, on utilise des matrices, d'abord Robinsoniennes, mais qui sont perturbées (ordre modifié, introduction de bruit blanc sur les valeurs.

Ce fichier contient les fonctions de génération de matrices aléatoires.
"""

from dissimilarity import *
import random

def generateRobinsonian(length, max_variation=5):
    """ Genere une matrice robinsonienne """
    # Genere une matrice remplie de 0 et taille "length"
    d = Dissimilarity([0]*(length*(length-1)//2))
    for j in range(1,d.n):
        for i in range(0,d.n - j):
            m = max(d[i,i+j-1],d[i+1,i+j])
            d[i,i+j] = m + random.randint(1,max_variation)
    return d
def permuteRandom(d):
    """Permute aléatoirement deux éléments d'une matrice"""
    i = random.randint(0, d.n-1)
    j = random.randint(0, d.n-1)
    d.permute(i,j)
def modifyRandom(d, min_noise=-5, max_noise=5):
    """Modifie une valeur choisie au hasard, en introduisant un bruit uniforme compris entre *min_noise* et *max_noise*\ ."""
    x = random.randint(0, d.n-1)
    y = random.randint(0, d.n-1)
    d[x,y] = max(0, d[x,y] + random.randint(min_noise, max_noise))
def randomize(d, repetition=50, min_noise=-5, max_noise=5):
    """Introduit de l'alétoire dans une matrice, en permutant pui modifiant ses valeurs au hasard (le processus est répété 50 fois)"""
    for i in range(repetition):
        permuteRandom(d)
        modifyRandom(d, min_noise, max_noise)
    
