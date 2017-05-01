#!/usr/bin/env python
#-*- coding : utf8 -*-


from math import sqrt
from matplotlib import pyplot as plt

from CentresDeMasse import CentreMasseProt
from CentresDeMasse import CentreMasseResidu


# Fonction prenant en argument le dictionnaire des conformations de la proteine
def DistanceCentre(dico_proteine):
	
	compteur_prot = 0
	distance_centre = {}
	centre_masse_prot = CentreMasseProt(dico_proteine)
	
	for Proteine in dico_proteine.keys():								# Pour chaque conformation de la proteine
		compteur_prot += 1
		for Chaine in dico_proteine[Proteine].keys():					# Pour chaque chaine de cette conformation
			if Chaine not in distance_centre.keys():					# Si elle n'existe pas
				distance_centre[Chaine] = {}							# On cree son dictionnaire correspondant
			
			for Residu in dico_proteine[Proteine][Chaine].keys():		# Pour chaque residu de la chaine
				if Residu not in distance_centre[Chaine].keys():		# S'il n'existe pas
					distance_centre[Chaine][Residu] = 0					# On cree son dictionnaire correspondant
				parcourt = CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])		# On stocke la valeur du centre de masse de ce residu
				
				
				# On calcule la distance au centre des differentes conformation de la proteine pour chaque residu et on fait leur moyenne
				# On pondere la moyenne obtenue au tour de boucle precedent pour eviter tout biais dans le calcul
				# On utilise cette formule pour ne pas avoir a refaire une boucle pour tout parcourir et diviser par le nombre de chaines de la proteine
				distance_centre[Chaine][Residu] = ((distance_centre[Chaine][Residu]*(compteur_prot-1)) +
													sqrt(
													pow(centre_masse_prot['x']-parcourt['x'],2) +
													pow(centre_masse_prot['y']-parcourt['y'],2) +
													pow(centre_masse_prot['z']-parcourt['z'],2))
													) / compteur_prot
	
	return distance_centre												# On retourne le dictionnaire des distances au centre moyennes de chaque residu

#Fonction prenant en argument le dictionnaire des distances au centre et tracant le graphe de ces distances moyennes pour chaque position de residu
def GraphDesDistances(dico_distances_centre):
	
	x = list()
	y = list()
	
	for chaine in sorted(dico_distances_centre.keys()): 				# Pour chaque chaine de la proteine
		for position in sorted(dico_distances_centre[chaine].keys()):	# Pour chaque position de residu
			x.append(position)											# Concatenation de la position du residu dans la liste x
			y.append(dico_distances_centre[chaine][position])			# Concatenation de chaque distance au centre moyenne dans la liste y
	
	plt.plot(x, y)														# Graphique des distances au centre moyennes en fonction des positions pour chaque residu 
	plt.title("Enfouissement des Residus")
	plt.xlabel("Positions")
	plt.ylabel("Eloignement du Centre")
	plt.show()
	
	return ()	
