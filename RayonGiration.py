#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Auteur : HENRIQUES Adrien et ARBES Hugo 
date : 28/04/2017
Programme Rayon de Giration : calcul du rayon de giration pour chaque conformation
et la distance séparant chaque centre de masse d'un résidu du centre de masse de la protéine
"""

from math import sqrt
from matplotlib import pyplot as plt

from CentresDeMasse import CentreMasseProt
from CentresDeMasse import CentreMasseResidu


# Fonction prenant en argument le dictionnaire des conformations de la proteine et renvoyant son rayon de giration
def RayonProt(dico_proteine):
	
	CdM_prot = CentreMasseProt(dico_proteine)							# On recupere le centre de masse pour les calculs de distances
	dico_rayon = {}														# On initialise le dictionnaire contenant les rayon de giration de chaque conformation
	
	for Proteine in dico_proteine.keys():								# Pour chaque conformation de la proteine
		dico_rayon[Proteine] = -1										# On initialise la valeur a -1 car une distance est toujours positive (necessaire pour un test plus bas)
		for Chaine in dico_proteine[Proteine].keys():					# Pour chaque chaine de la conformation
			for Residu in dico_proteine[Proteine][Chaine].keys():		# Pour chaque residu de la chaine
				
				# On stocke son centre de masse de ce residu
				parcourt = CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])
				
				# On calcule la distance entre le centre de masse et chaque residu
				distant = sqrt(pow(CdM_prot[Proteine]['x']-parcourt['x'],2) + pow(CdM_prot[Proteine]['y']-parcourt['y'],2) + pow(CdM_prot[Proteine]['z']-parcourt['z'],2))
				
				# Si cette distance est superieure a la plus grande trouvee pour les conformations observees precedemment
				if distant > dico_rayon[Proteine] :						# Pour la premiere valeur on rentre toujours dans le test car distant est positif
					dico_rayon[Proteine] = distant						# On stocke cette distance (correspondant au rayon de giration)
	
	return dico_rayon													# On retourne la plus grande distance (rayon de giration)


# Fonction prenant en argument le dictionnaire des rayons de giration et affichant le graphe de ces rayons pour chaque conformation de la proteine
def GraphDesRayons (dico_rayon_giration):
	
	x = list()															# Liste qui contiendra les indices des conformations de la proteine
	y = list()															# Liste qui contiendra les valeurs des rayons de giration de chaque conformation
	
	for prot in sorted(dico_rayon_giration.keys()): 					# Parcourt les conformations dans le dictionnaire des rayons de giration
		x.append(prot)													# Concatenation des indices des conformations rencontrees dans la liste x
		y.append(dico_rayon_giration[prot]) 							# Concatenation des rayons de giration dans la liste y
	
	plt.plot(x, y)  													# Trace le graphique des angles de giration en fonction des conformations
	plt.title("Graphique des Rayons de Giration")
	plt.xlabel("Conformations")
	plt.ylabel("Rayons de Giration")
	plt.show()
	
	return ()
