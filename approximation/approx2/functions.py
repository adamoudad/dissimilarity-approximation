#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
****************
Autres fonctions
****************
"""
from dissimilarity import *
import numpy as np
import math
import copy
# garbage collector
import gc



def take_connected_graphe(C_graphe, i, j, sub):
        '''Trouver le sous graphe passé par C_graphe(i,j) par reccurence'''
        # Les 4 directions de cette position
        jump=[[i-1, j], [i, j+1], [i+1, j], [i, j-1]]
        for direction in range(0, 4):
                if C_graphe[i][j][direction]:
                        next_x=jump[direction][0]
                        next_y=jump[direction][1]
                        # prendres cette contrainte
                        sub.append([min(i, next_x), max(j, next_y), max(i, next_x), min(j, next_y)])
                        # supprimer l'arrêt
                        C_graphe[i][j][direction]=False
                        #supprimer dans l'arrêt dans l'autre côté
                        next_direction=(direction+2)%4
                        C_graphe[next_x][next_y][next_direction]=False
                        # sauter au prochaint point
                        take_connected_graphe(C_graphe, next_x, next_y, sub)
        return True


def base_orthonorme(Ca_sub, n):
        '''
        Générer une suite des bases orthonormées de l'hyperplane formée par les contraintes.
		
        * Entree:Ca_sub. C'est une liste des listes, dont chaque liste contient des listes de forme [i, j, k, l].
        * Sortie: Une liste des bases orthonormées. 
                
        Chaque base est représentée par l'example: [valeur, [[i, j], [k, l], ...)]] avec : 
        
        * [[i, j], [k, l], ...] : indices des composantes non nulles
        * valeur: valeur des composantes non nulles
        
        Complexité:	O(n^2+m*g)
        '''
        Bases=[ ]
        # Chaque sous groupe de Ca_sub corresponde à une base de l'hyperplanne 
        # d_contraint enregistre les contraints qui s'agit de plusieur dimension.
        d_contraint=[]
        for sub in Ca_sub:
                b=[]
                # pour chaque base, trouver des dimension non nulle
                for c in sub:
                        if [ c[0], c[1] ] not in b:
                                b.append([ c[0], c[1] ])
                        if [ c[2], c[3] ] not in b:
                                b.append([ c[2], c[3] ])
                d_contraint.extend(b)
                Bases.append( [ 1/ math.sqrt(len(b)), b ] )

        # Les autres bases dont une seul dimension prend la valeur non nulle.

        for i in range(0, n-1):
                for j in range(i+1, n):
                        if [i, j] not in d_contraint:
                                b=[]
                                b.append([i,j])
                                Bases.append( [1.0, b] )
        return Bases


def projection(D0, Bases):
        '''
        Projecter D0 sur une suite de bases orthonormées.

        Entree: 		D0 est une matrices triangulaire.
                        Bases. Une liste des bases orthonormées de l'hyperplane dans l'espace où se trouve D0. 
                        Chaque base est représentée par l'example: [valeur, [[i, j], [k, l], ...] ], 
                        avec : 
                        [[i, j], [k, l], ...] : indices des composantes non nulles
                        valeur: valeur des composantes non nulles

        Sortie: 		Dr. C'est le point de projection
        Complexité:	O(n^2*g)
        '''
	# initialiser Dr
        Dr=copy.deepcopy(D0)
	# Projection= sum((D0 . b) b);  b est un vecteur dans Bases
	# On va projecter D0 sur chaque base b:
	# On remarque que toutes les base b ont des composantes non nulles totalement différentes( n'ont pas de composantes non nulles en commun)
        for b in Bases:
		#b_value: valeur des composantes non nulles de base b
                b_value=b[0]
		# b_points est les composantes non nulles de base b
                b_points=b[1]
		# D0_useful enregistre les valeurs de D0 dont l'indice est dans b_points; ce sont des indices utiles pour la projection sur b
                D0_useful=[]
                for p in b_points:
                        row=p[0]
                        col=p[1]
                        D0_useful.append(D0[row, col])
		# projection de D0 sur b
                proj=sum(D0_useful)*b_value*b_value
                for p in b_points:
                        row=p[0]
                        col=p[1]
                        Dr[row, col]=proj
        return Dr



def multiplicateurs(Dr, D0, Ca, Ca_sub, Bases):
        """
        Résoudre les multiplicateur de lagrange
        
        * Entree: D0 et Dr (matrice initiale) sont des objets de la classe Dissimularité 
                ** Ca_sub est des contraintes regroupées
                ** Bases est une liste de bases de l'hyperplane des contraintes
        * Sortie: Une liste de multiplicateurs Lambda sous forme: [[lambda, [i, j, k, l]], ...], avec:  [i, j, k, l], contrainte; lambda, multiplicateur corresponde à cette contrainte.
 
        Complexité:	O(m*g^2) ??
        """
        Lambda=[]
        #on va résoudre les mulitiplicateurs dans des sous espaces, chaque sous espaces corresponde à une sous-groupe de contraintes dans Ca_sub 
        Len=len(Ca_sub)
        for i in range(0,Len):
		#Dim: les indices de cette sous espaces 
                Dim=Bases[i][1][:]
		# supprimer les contraintes inutiles pour chaque Ca_sub
                Contraints=delete_ctr_useless(Ca_sub[i])

                n=len(Dim)
                m=len(Contraints)

		# preparation de la matrice de contrainte M_contraints
                M_contraints=np.zeros((n, m), dtype='float') #M_contraints est de type array
                for j in range(0, m):
                        id1=Dim.index([ Contraints[j][0], Contraints[j][1] ])
                        id2=Dim.index([ Contraints[j][2], Contraints[j][3] ])
                        M_contraints[id1, j]=1.0
                        M_contraints[id2, j]= -1.0
		# Supprimer le dernière ligne de M_contraintes pour le rendre de dimension carrée.
                M_contraints_sq=np.delete(M_contraints, n-1, 0)
		# relaease l'espace de stockage
                del M_contraints
                gc.collect()
		# préparation de la vecteur Valeurs en type numpy array.
                Valeurs=np.zeros(n-1, dtype=float)
                for j in range(0, n-1):
                        position=Dim[j]
                        Valeurs[j]= 2.0*(Dr[position[0], position[1]] - D0[position[0], position[1]])

		# Résoudre les multiplicateur lambda de cette sous-groupe de contraint.
                sol=np.linalg.solve(M_contraints_sq, Valeurs)

		# Enregistrer les multiplicateur dans son contrainte correspondant.
                for j in range(0, m):
                        Lambda.append([ sol[j], Contraints[j] ])
        return Lambda



# supprimer les contraintes inutiles
def delete_ctr_useless(Ctr_ini):
	Contraints=[]
	dim_exist=[]
	for c in Ctr_ini:
		if ([c[0], c[1]] not in dim_exist) or  ( [c[2], c[3]] not in dim_exist):
			Contraints.append(c)

			if [c[0], c[1]] not in dim_exist:
				dim_exist.append( [c[0], c[1]] )
			if [c[2], c[3]] not in dim_exist:
				dim_exist.append( [c[2], c[3]] )

	return Contraints













