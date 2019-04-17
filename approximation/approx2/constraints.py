#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
***********************
Gestion des contraintes
***********************

Fonctions de manipulation des contraintes

'''

from dissimilarity import *
from approximation.approx2.functions import *

import numpy as np


def ini_contraints(D0):
        """
        Fonction: ini_contraints
        Description: Parcourir D0 et trouver les contraintes non satisfaits.
        Entree: D0 (matrice initiale) est un objet de la classe Dissimularité
        Sortie: Ca. Il est une liste des listes dont chacune a la forme (i, j, i',j'), soit (i=i' et j=j'-1) soit (i=i'+1 et j=j'),avec la condition D0[i,j]<D0[i',j']
        Complexité:	O(n^2)
        """
        n=D0.n
        Ca=[]
        for i in range(0, n-2):
                for j in range(i+2, n):
                        if D0[i, j]<D0[i, j-1]:
                                Ca.append( [i, j, i, j-1] )
                        if D0[i, j]<D0[i+1, j]:
                                Ca.append( [i, j, i+1, j] )
        return Ca

def regroupe_contraints(Ca, n):
        """
        regrouper les containtes de Ca en sous_groupe: soitent [i, j, k, l] et [i', j', k', l'] deux listes de Ca, 
        
        * si l'une des conditions suivantes est satisfaite: (i=i', j=j') ou (i=k', j=l') ou (k=i', l=j') ou (k=k', l=l'),	
        * on ajoute alors les deux contraintes [i, j, k, l] et [i', j', k', l'] dans la même sous-liste.  


        Entree: Ca. C'est une liste dont chaque élément a la forme [i, j, k, l], représentant la contrainte D(i, j)>D(k,l), n. Dimension de la matrice
        Sortie: Ca_sub. C'est une liste des listes, dont chaque liste est une groupe de contrainte connectée.
        
        Complexité:	O(n^2)
        """
	# Représenter Ca par un graphe, chaque contrainte est représentée par un arrêt
	# Les éléments(points) de ce graphe est les points d'une matrice triangulaire de dimension n*n 
        C_graphe=np.zeros((n, n, 4), dtype=bool)
        for c in Ca:
		#prendre les quatres indices de chaque contrainte: (row1, col1), (row2, col2)
                row1=c[0]
                col1=c[1]
                row2=c[2]
                col2=c[3]
		# Au chaque point de C_graphe, quatre valeurs de type bool sont attribuées aux quatre direction [up, right, down, left], pour indiquer s'il existe un arrêt dans cette direction.
		# C'est un arrêt horizontal:
                if row1==row2 and col1==col2+1 :
                        C_graphe[row1][col1][3]=True
                        C_graphe[row2][col2][1]=True
                # C'est un arrêt vertical:
                if row1==row2-1 and col1==col2 :
                        C_graphe[row1][col1][2]=True
                        C_graphe[row2][col2][0]=True

	# Séparer cette graphe en sous-graphes qui ne sont pas connectées entre eux. Mais les points dans une graphe sont connectés.
        Ca_sub=[]
        # Parcourir ce graphe.
        for k in range(0, n-2):
                # i, j désignent le position dans le matrice
                for i in range( 0, k+1 ):
                        j=i+n-1-k
                        sub=[]
                        take_connected_graphe(C_graphe, i,j, sub) 
                        if len(sub)>0:
                                Ca_sub.append(sub)
        return Ca_sub


# complexité: O(m^2)
def remove_contraints(Ca, Lambda):
    """ Retire de Ca les contraintes de Robinson dont le Lambda est négatif """
    removed=False
    # print 'contraintes_removed'
    nb=0
    for lambda_elements in Lambda:
        if lambda_elements[0]<0:
            #print lambda_elements[1]
            Ca.remove(lambda_elements[1])
            removed=True
            nb+=1
    # print  nb
    return removed

# complexité: O(n^2*m)
def addContraints(d, Ca):
    """ Ajoute les conditions de Robinson qui ne sont pas encore satisfaites à Ca """
    addedConstraints = False
    precision=0.0001
    nb=0
    # print 'Contraintes ajoutees'
    for i in range(d.n-2):
        for j in range(i+2, d.n):
            if (d[i,j] < d[i,j-1] - precision) and ([i,j,i,j-1] not in Ca):
                # print [i,j,i,j-1]
                Ca.append([i,j,i,j-1]) 
                addedConstraints = True
                nb+=1
            if (d[i,j] < d[i+1,j] - precision) and ([i,j,i+1,j] not in Ca):
                # print [i,j,i+1,j]
                Ca.append([i,j,i+1,j])
                addedConstraints = True
                nb+=1
                # print 'Ca in addConstraints\n', Ca
                # print 'addedConstraints\n', addedConstraints
    # print nb
    return addedConstraints


def remove_one_contraints(Ca, Lambda):
    """ Retire de Ca une contraintes de Robinson dont le Lambda est négatif """
    removed=False
    # print 'contraintes_removed'
    for lambda_elements in Lambda:
        if lambda_elements[0]<0:
            # print lambda_elements[1]
            Ca.remove(lambda_elements[1])
            removed=True
            break
    return removed

def addOneContraints(d, Ca):
    """ Ajoute une conditions de Robinson qui ne sont pas encore satisfaites à Ca """
    addedConstraints = False
    precision=0.0001
    # print 'Contraintes ajoutees'
    for i in range(d.n-2):

        for j in range(i+2, d.n):
            if (d[i,j] < d[i,j-1] - precision) and ([i,j,i,j-1] not in Ca):
                # print [i,j,i,j-1]
                Ca.append([i,j,i,j-1]) 
                addedConstraints = True
                return addedConstraints
            if (d[i,j] < d[i+1,j] - precision) and ([i,j,i+1,j] not in Ca):
                # print [i,j,i+1,j]
                Ca.append([i,j,i+1,j])
                addedConstraints = True
                # print 'Ca in addConstraints\n', Ca
                # print 'addedConstraints\n', addedConstraints
                return addedConstraints
    return addedConstraints

def remove_part_contraints(Ca, Lambda, percent):
    """ Retire de Ca un certain pourcentage des contraintes de Robinson dont le Lambda est négatif """
    import math
    removed=False
    # print 'contraintes_removed'
    nb=0
    nb_negative= len([ x[0]  for x in Lambda if x[0] < 0 ])
    nb_should_remove=math.ceil(nb_negative*percent)
    for lambda_elements in Lambda:
        if lambda_elements[0]<0:
            #print lambda_elements[1]
            Ca.remove(lambda_elements[1])
            removed=True
            nb+=1
        if nb > nb_should_remove+1:
            break
    # print  nb
    return removed

def addPartContraints(d, Ca, percent):
    """ Ajoute un certain pourcentage des conditions de Robinson qui ne sont pas encore satisfaites"""
    import math
    addedConstraints = False
    precision=0.0001
    nb=0
    nb_add_total=0
    for i in range(d.n-2):
        for j in range(i+2, d.n):
            if ((d[i,j] < d[i,j-1] - precision) and ([i,j,i,j-1] not in Ca)) or ( (d[i,j] < d[i+1,j] - precision) and ([i,j,i+1,j] not in Ca)): nb_add_total+=1
    nb_should_add= math.ceil(nb_add_total*percent)

    for i in range(d.n-2):
        for j in range(i+2, d.n):
            if (d[i,j] < d[i,j-1] - precision) and ([i,j,i,j-1] not in Ca):
                Ca.append([i,j,i,j-1]) 
                addedConstraints = True
                nb+=1
            if (d[i,j] < d[i+1,j] - precision) and ([i,j,i+1,j] not in Ca):
                Ca.append([i,j,i+1,j])
                addedConstraints = True
                nb+=1
            if nb > nb_should_add+1.0:
                return addedConstraints
    return addedConstraints


