#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
--------------
Tri par pivots
--------------
Ce fichier contient l'algorithme de Tri des éléments par pivots
"""
import random
from dissimilarity import *

def sortElements(d, start=0, end=-1, cur_class=-1, last_pivot=-1, verbose=False):
    """ Algorithme de réarrangement des éléments
    Le but est de trouver un arrangement des éléments qui pourrait se rapporcher au maximum d'un ordre compatible

    Paramètres : 
    d : la dissimilarité
    start : l'indice de début de la sous-matrice considérée
    end : l'indice de fin de la sous-matrice considérée (non inclus)
    verbose : affiche toutes les permutations et cas terminaux
    Renvoi: rien
    """
    # Initialisation
    if end == -1: end = d.n
    # Choix des pivots
    if (end - start) > 1:
        x,y = choosePivots(d.order[start:end])
        # On organise les deux pivots, par rapport au pivot précédent
        if cur_class == 0:
            if getClass(d,x,y,last_pivot) == 0:
                x,y = y,x
                d.permute(x,y, by="element")
                if verbose: print(y,",",x," -> ",x,",",y)
        elif cur_class == 1:
            if getClass(d,x,y,last_pivot) == 2:
                x,y = y,x
                d.permute(x,y, by="element")
                if verbose: print(y,",",x," -> ",x,",",y)
        elif cur_class == 2:
            if getClass(d,x,y,last_pivot) == 2:
                x,y = y,x
                d.permute(x,y, by="element")
                if verbose: print(y,",",x," -> ",x,",",y)
    # On vérifie la taille de la sous-matrice
    if (end - start) > 2:
        # On prépare la liste qui stockera les éléments dans leurs classes respectives
        C = ([],[],[])
        # On parcourt les éléments à arranger, et on les classe
        for z in d.order[start:end]:
            # Pour les éléments qui ne sont pas les pivots
            if z not in (x,y):
                # On récupère la classe à laquelle  k appartient
                c = getClass(d,x,y,z)
                C[c].append(z)
        # Affichage des infos sur cet appel récursif
        if verbose:
            print("Pivots : ",x,",",y, " >> ",d.order[start:end], " --> ", C[0]+[x]+C[1]+[y]+C[2])
        # Les C[i] contiennent maintenant les indices des éléments classés dans la classe i, on applique les modifications à la dissimilarité d comme suit :
        # Réarrangement de l'ordre dans la dissimilarité
        for index,element in enumerate(C[0]+[x]+C[1]+[y]+C[2]):
            d.permute(index+start, d.index(element), by="index")
        # Appels récursifs sur chaque classe
        sortElements(d, start, len(C[0]), 0, x, verbose)
        sortElements(d, start + len(C[0]) + 1, start + len(C[0]) + 1 + len(C[1]), 1, x, verbose)
        sortElements(d, start + len(C[0]) + 2 + len(C[1]), end, 2, y, verbose)
    else:
        if verbose:
             print("Cas terminal !")

def choosePivots(order, mode="random"):
    """ Choisit deux indices, pour les pivots dans l'algorithme de tri
    Params :
    order : liste des éléments dans laquelle piocher les deux pivots
    mode : indique la façon de piocher les pivots ; random(aléatoirement), middle(couper en trois parties équilibrées)
    Renvoi: Les deux pivots
    """
    if mode == "middle":
        i = len(order)//3
        j = 2*i
        return (order[i],order[j])
    elif mode == "random":
        x,y = random.sample(order,2)
        if order.index(x) < order.index(y): return x,y
        else: return y,x

def getClass(d,x,y,z):
    """ Renvoie la classe de l'élément z pour les pivots x et y """
    m = max(d[x,y],d[x,z],d[y,z])
    if d[y,z] == m:
        return 0
    elif d[x,y] == m:
        return 1
    else:
        return 2

if __name__ == "__main__":
    print("Test de tri")
    # d = Dissimilarity([ 1,2,3,4,5,1,2,3,4,1,2,3,1,2,1 ])
    # d = Dissimilarity([ 5,4,5,10,3,6,1,3,9,8,2,5,4,7,1 ])
    # d = Dissimilarity([ 4, 3, 6, 1, 8, 5  ])
    # d = Dissimilarity([1, 5, 7, 7, 3, 4, 6,3, 3,2])
    # d = Dissimilarity([10, 5, 2, 3, 3, 4, 10,3, 3,2])
    d = Dissimilarity([1,2,3,1,2,1])
    print(d)
    d.permute(2,3)
    print(d)
    sortElements(d)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    # sortElements(d, 0, d.n)
    print(d)
    print(d.order)
