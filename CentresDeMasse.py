#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Auteur : HENRIQUES Adrien et ARBES Hugo 
date : 28/04/2017
Programme Centres de Masse : calcule des centres de masse pour chaque conformation et pour chaque residu"""

# Fonction prenant en argument le dictionnaire des conformations de la proteine et retournant son centre de masse
def CentreMasseProt(dico_proteine):
	
	nbr_residus = 0
	x = 0
	y = 0
	z = 0
	centre_masse_prot = {}
	
	# Ajout de toutes les coordonnees des centres de masse de tous les residus
	for Proteine in dico_proteine.keys():								# Pour chaque conformation de la proteine
		centre_masse_prot[Proteine] = {}								# On cree un dictionnaire contenant son centre de masse
		for Chaine in dico_proteine[Proteine].keys():					# Pour chaque chaine de cette conformation
			for Residu in dico_proteine[Proteine][Chaine].keys():		# Pourchaque residu de cette chaine
				
				# On ajoute les coordonnees du centre de masse
				x += CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])['x']
				y += CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])['y']
				z += CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])['z']
				nbr_residus += 1
				
		# On fait la moyenne des coordonnees
		centre_masse_prot[Proteine]['x'] = x / nbr_residus
		centre_masse_prot[Proteine]['y'] = y / nbr_residus
		centre_masse_prot[Proteine]['z'] = z / nbr_residus
	
	return centre_masse_prot											# On retourne le centre de masse

# Fonction prenant en argument le dictionnaire d'un residu et qui calcule son centre de masse
def CentreMasseResidu(dico_residu):
	
	masse_molaire = dict()												# Creation d'un dictionnaire pour les references de poids atomiques
	masse_molaire = {	"H": float(1.0079),		"C": float(12.0107),	# Initialisation
						"O": float(15.9994),	"N": float(14.0067),
						"S": float(32.065)}
	x = 0
	y = 0
	z = 0
	masse_totale_residu = 0
	centre_masse_residu = {}
	
	# Ajout de toutes les coordonnees des atomes du residu
	for atom_compare in dico_residu.keys():								# Pour chaque atome du residu observe
		masse_atome = 0
		for atom_reference in masse_molaire.keys():						# Pour chaque atome dans le dictionnaire de reference
			if atom_reference in atom_compare:							# Si l'atome de reference fait partie du nom de l'atome qu'on observe
																		# On utilise "in" car dans le cas des H on peut avoir H, OH ou NH
				masse_atome += masse_molaire[atom_reference]			# On ajoute le poids de cet atome au poids du residu
				
			
		# On multiplie les coordonnees par le poids pour avoir le centre de masse du residu
		x += masse_atome * float(dico_residu[atom_compare]['x'])
		y += masse_atome * float(dico_residu[atom_compare]['y'])
		z += masse_atome * float(dico_residu[atom_compare]['z'])
		masse_totale_residu += masse_atome
	
	# On divise par la masse cumulee des atomes pour avoir les coordonnees du centre de masse du residu
	centre_masse_residu['x'] = x / masse_totale_residu
	centre_masse_residu['y'] = y / masse_totale_residu
	centre_masse_residu['z'] = z / masse_totale_residu
	
	return centre_masse_residu											# On retourne le centre de masse
