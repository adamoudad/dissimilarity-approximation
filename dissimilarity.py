#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Une dissimilarité est représentée par un ensemble de vecteurs et une matrice, qui associe à chaque pair de vecteurs, leur dissimilarité.

Une matrice de dissimilarité, en plus d'être à valeur positive et diagonale nulle, est symétrique.

Ainsi, la classe Dissimilarity implémente une représentation optimisée pour la manipulation de Dissimilarités.


"""
from array import *
import random


def triangularNumber(n):
    """ Renvoie i si n est le ième nombre triangulaire, renvoie -1 si n n'est pas triangulaire.(cf : https://fr.wikipedia.org/wiki/Nombre_triangulaire ) """
    m = n
    i = 0                   # Compteur
    while m > 0:
        # On retire à chaque étape i, qui vaut : 1, puis 2, puis 3, ...
        i += 1
        m -= i
    # n est triangulaire ssi m est nul à la fin de l'opération précédente
    # Seulement si n est triangulaire, on renvoie i, le dernier nombre de la suite.
    return i if (m == 0) else -1


class Dissimilarity:
    def __init__(self,matrix, order=[]):
        """ Initialise une dissimilarité
        matrix doit décrit le triangle supérieur de la matrice de dissimilarité.
        order peut être renseigné si on veut créer une dissimilarité avec un ordre précis. Si order n'a pas les bonnes dimensions, on initialisera la matrice avec l'ordre canonique."""
        # Test si la tableau de valeur "matrix" fournie définit bien une matrice triangulaire <=> la taille de "matrix" est un nombre triangulaire
        n = triangularNumber(len(matrix))
        if n != -1 : self.n = n + 1 # n+1 car la diagonale est nulle, et dans "matrix", on a seulement mis le triangle supérieur (ou inférieur)
        else: raise IndexError("Matrice de mauvaises dimensions fournie")
        self.matrix = array('d', matrix)
        if order and (len(order) == self.n): self.order = list(order)
        else: self.order = list(range(self.n))
    def __len__(self):
        return self.n
    def __getitem__(self, elements):
        """ Renvoie la dissimilarité entre les éléments x et y"""
        x,y = elements          # Il s'agit des éléments ! pas des positions
        # i,j = self.order[i], self.order[j]
        if x >= self.n or x < 0 or y >= self.n or y < 0:
            raise IndexError("L'indice sort des limites du tableau")
        if x == y: return 0
        # Version lecture du triangle haut
        elif x > y: return self.matrix[x-1 + (self.n-1)*(self.n-2)//2 - (self.n-y-1)*(self.n-y-2)//2]
        else: return self.matrix[y-1 + (self.n-1)*(self.n-2)//2 - (self.n-x-1)*(self.n-x-2)//2]
    def __iter__(self):
        """ Permet de rendre l'objet "iterable", il peut donc être parcouru par une boucle <<for i in d:>>, être converti en liste, ... par exemple """
        for line in  [ [ self[x,y] for y in range(self.n) ] for x in range(self.n) ] :
            yield line
    def __setitem__(self, elements, value):
        """ Modifie à "value" la dissimilarité entre les éléments x et y"""
        x,y = elements          # Il s'agit des éléments ! pas des positions
        if x >= self.n or x < 0 or y >= self.n or y < 0:
            raise IndexError("L'indice sort des limites du tableau")
        if x > y: self.matrix[x-1 + (self.n-1)*(self.n-2)//2 - (self.n-y-1)*(self.n-y-2)//2] = value
        elif x < y: self.matrix[y-1 + (self.n-1)*(self.n-2)//2 - (self.n-x-1)*(self.n-x-2)//2] = value
    def index(self,x):
        """ Retourne l'indice/position de l'élément x """
        return self.order.index(x)
    def __reversed__(self):
        """ Renverse l'ordre des éléments """
        self.order.reverse()
    def __str__(self):
        """ Renvoie une chaine de caractère de la matrice de la dissimilarité """
        block_size = max([ len(str(v)) for v in self.matrix ])
        spacing = 2
        output = ""
        for i in range(0,self.n):
            output += "[" + " "*(block_size + spacing)*(i+1)
            if i != self.n-1:
                for j in range(i+1,self.n):
                    c = str(self[self.order[i],self.order[j]])
                    output += " "*spacing + c + " "*(block_size - len(c))
                output += " "*spacing + "]\n"
            else: output += " "*spacing + "]"
        return output
    def __add__(self, other):
        """ Additionne deux dissimilarités, sur leurs valeurs avec une paire d'éléments, indépendamment de l'ordre des éléments """
        if self.n != other.n: raise IndexError("Le nombre d'éléments des deux matrices à additionner est différent")
        else: return Dissimilarity( [ self.matrix[i] + other.matrix[i] for i in range(len(self.matrix)) ] )
    def __pow__(self, p):
        """ Élève toutes les valeurs de la dissimilarité à la puissance p"""
        result_matrix = [0]*len(self.matrix)
        for i in range(len(self.matrix)):
            result_matrix[i] = self.matrix[i] ** p
        return Dissimilarity(result_matrix, self.order)
    def permute(self, i, j,by="index"):
        """ Permute les éléments de la dissimilarité
        modes : 
        element : premute les éléments selon leur nom (entier)
        index : permute les éléments selon leur indice (position dans "order")
        """
        if by == "element":
            i,j = self.index(i),self.index(j)
        self.order[j],self.order[i] = self.order[i],self.order[j]
    def isCompatibleOrder(self):
        """ Test si l'ordre est un ordre compatible pour la dissimilarité """
        isCompatible = True
        precision=0.0001
        for i in range(self.n):
            for j in range(self.n-1):
                if i < j: isCompatible = (self[self.order[i],self.order[j]] <= self[self.order[i],self.order[j+1]] + precision)
                elif i > j: isCompatible = (precision + self[self.order[i],self.order[j]] >= self[self.order[i],self.order[j+1]])
                if not(isCompatible): return False
        return True
    def update(self, new_order=[]):
        """ Met à jour les valeurs de "matrix", et réinitialise l'ordre """
        if not new_order: new_order = list(range(self.n))
        if len(new_order) != len(self.order): raise IndexError("Le nombre d'éléments dans la liste fournie ne correspond pas à celui de cette dissimilarité")
        # Copie de la matrice
        old_matrix = [[ self[x,y] for x in range(self.n) ] for y in range(self.n) ]
        for x in range(self.n):
            for y in range(self.n):
                self[new_order[self.index(x)],new_order[self.index(y)]] = old_matrix[x][y]
        self.order = new_order

## Tests ##
if __name__=="__main__":
    # D = Dissimilarity([ 1,2,3,4,5,1,2,3,4,1,2,3,1,2,1 ])
    # D = Dissimilarity([ 1,2,3,4,1,2,3,1,2,1 ])
    D = Dissimilarity(list(range(10*9//2)))
    # D = Dissimilarity(random.sample(range(100),55))
    print(D)
    print(D.isCompatibleOrder())
