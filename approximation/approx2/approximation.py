#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dissimilarity import *
from approximation.approx2.constraints import *
from approximation.approx2.functions import *
# garbage collector
import gc

def approximate(D0):
        """
        Approxime D0 par une matrice Robisonienne Dr en Cardant l'ordre des éléments de D0
        
        * Entree: D0 (matrice initiale) est un objet de la classe Dissimularité 
        * Sortie: Dr (matrice Robisonnienne)
        
        Complexité: O(n^2*m+g^2*m)* nb_boucle,
        
        * n: nombre d'élément
        * m: nombre de contraintes non satisfaites
        * g: taille de graphe séparée
        * nb_boucle: nombre de boucles pour trouver Dr
        
        """
	# Ca enregistre les contraintes atteintes sous une liste des listes. Chaque liste a la forme [i, j, i',j'],  qui désigne la contrainte D0[i,j]>=D0[i',j']. ici D0[i, j] et D0[i, j] sont deux valeurs voisines: soit (i=i' et j=j'-1),  ou (i=i'+1 et j=j').

        # Initialisation de Ca
        Ca=ini_contraints(D0) #complexity: O(n^2)
        # n: dimension de matrice, 
        n=D0.n
        first_bool=True
        boocle=0
        # while True:
        while boocle<50:
        #gestion de stockage RAM
                boocle+=1
                if first_bool==False:
                        del Ca_sub, Bases, Dr, Lambdas
                        gc.collect()
                first_bool=False

                #Regrouper les contraintes Ca sous des groupes none connectés. Ca_sub est une liste des listes, et a la forme [[[i,j,i',j'],...],...]. 
                #complexite: O(n^2)
                Ca_sub = regroupe_contraints(Ca, n)

                # Générer une suite de bases orthonormées de l'hyperplane de Ca(toutes les contraintes atteintes).
                # complexite: O(n^2+m*g)
                Bases = base_orthonorme(Ca_sub, n)

                # Projecter D0 sur l'hyperplanes formé par les contraintes de Ca.
                # Complexité:	O(n^2*g)
                Dr = projection(D0, Bases)

                
                #résoudre les multiplicateurs de Ca
                # Complexité:	O(m*g^2) ??
                Lambdas = multiplicateurs(Dr, D0, Ca, Ca_sub, Bases)

                ##--------------------------------------------------------------------------------------------------------------------
                ##_________________Les méthodes pour ajouter ou supprimer des contraintes:________
                ##---------------------------------------------------------------------------------------------------------------------
        
                ##---------------------------------methode 1-------------------------------------------------
                # if (not remove_contraints(Ca, Lambdas)) and (not addContraints(Dr, Ca)):
                # 	break
		##test resultat: possible de rentre dans boucle infinie

		##---------------------------------methode 2-------------------------------------------------
		# complexité: O(m^2 )
                removed=remove_contraints(Ca, Lambdas)
                # complexité: O(n^2*m)
                added=addContraints(Dr, Ca)
                if (not removed ) and (not added ):
                        break
                #test resultat: possible de rentre dans boucle infinie, un peu plus rapide que methode 1

                ##---------------------------------methode 3-------------------------------------------------
                # removed=remove_one_contraints(Ca, Lambdas)
                # added=addOneContraints(Dr, Ca)
                # if (not removed ) and (not added ):
                # 	break		

		##---------------------------------methode 4-------------------------------------------------
		# if boocle < 100:
		# 	removed=remove_contraints(Ca, Lambdas)
		# 	added=addContraints(Dr, Ca)
		# 	if (not removed ) and (not added ):
		# 		break	
		# else:
		# 	removed=remove_one_contraints(Ca, Lambdas)
		# 	added=addOneContraints(Dr, Ca)
		# 	if (not removed ) and (not added ):
		# 		break		

		##---------------------------------methode 5-------------------------------------------------
		# percent=0.1
		# removed=remove_part_contraints(Ca, Lambdas, percent)
		# added=addPartContraints(Dr, Ca, percent)
		# if (not removed ) and (not added ):
		# 	break	

        return Dr
