#!/usr/bin/env python
#-*- coding : utf8 -*-

from math import sqrt



def CentreMasseProt(dico_proteine):
	
	nbr_residus = 0
	x = 0
	y = 0
	z = 0
	centre_masse_prot = {}
	
	# Ajout de toutes les coordonnees des centres de masse de tous les residus
	for Proteine in dico_proteine.keys():
		for Chaine in dico_proteine[Proteine].keys():
			for Residu in dico_proteine[Proteine][Chaine].keys():
				x += CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])['x']
				y += CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])['y']
				z += CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])['z']
				nbr_residus += 1
				
	centre_masse_prot['x'] = x / nbr_residus
	centre_masse_prot['y'] = y / nbr_residus
	centre_masse_prot['z'] = z / nbr_residus
	
	return centre_masse_prot


def CentreMasseResidu(dico_residu):
	
	masse_molaire = dict()
	masse_molaire = {	"H": float(1.0079),		"C": float(12.0107),
						"O": float(15.9994),	"N": float(14.0067),
						"S": float(32.065)}
	x = 0
	y = 0
	z = 0
	masse_totale_residu = 0
	centre_masse_residu = {}
	
	# Ajout de toutes les coordonnees des atomes du residu
	for atom_compare in dico_residu.keys():
		masse_atome = 0
		for atom_reference in masse_molaire.keys():
			if atom_reference in atom_compare:
				masse_atome += masse_molaire[atom_reference]
				
		x += masse_atome * float(dico_residu[atom_compare]['x'])
		y += masse_atome * float(dico_residu[atom_compare]['y'])
		z += masse_atome * float(dico_residu[atom_compare]['z'])
		masse_totale_residu += masse_atome
	
	# On divise par la masse cumulee des atomes pour avoir les coordonnees du centre de masse du residu
	centre_masse_residu['x'] = x / masse_totale_residu
	centre_masse_residu['y'] = y / masse_totale_residu
	centre_masse_residu['z'] = z / masse_totale_residu
	
	return centre_masse_residu



def RayonProt(dico_proteine):
	
	CdM_prot = CentreMasseProt(dico_proteine)
	le_plus_loin = 0
	
	for Proteine in dico_proteine.keys():
		for Chaine in dico_proteine[Proteine].keys():
			for Residu in dico_proteine[Proteine][Chaine].keys():
				parcourt = CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])
				distant = sqrt(pow(CdM_prot['x']-parcourt['x'],2) + pow(CdM_prot['y']-parcourt['y'],2) + pow(CdM_prot['z']-parcourt['z'],2))
				if distant > le_plus_loin :
					le_plus_loin = distant
	
	return le_plus_loin



def DistanceCentre(dico_complexe_prot):
	
	compteur_prot = 0
	distance_centre = {}
	centre_masse_prot = CentreMasseProt(dico_complexe_prot)
	
	for Proteine in dico_complexe_prot.keys():							# Pour chaque proteine
		compteur_prot += 1
		for Chaine in dico_complexe_prot[Proteine].keys():				# Pour chaque chaine de la proteine
			if Chaine not in distance_centre.keys():					# Si elle n'existe pas
				distance_centre[Chaine] = {}							# On cree son dico correspondant
			
			for Residu in dico_complexe_prot[Proteine][Chaine].keys():	# Pour chaque residu de la chaine
				if Residu not in distance_centre[Chaine].keys():		# S'il n'existe pas
					distance_centre[Chaine][Residu] = 0					# On le cree
				parcourt = CentreMasseResidu(dico_complexe_prot[Proteine][Chaine][Residu])		# On stocke la valeur du centre de masse de ce residu
				
				# On calcule la distance au centre des differentes conformation de la proteine pour chaque residu et on fait leur moyenne
				# On utilise cette formule pour ne pas avoir a refaire une boucle pour tout parcourir et diviser par le nombre de chaines de la proteine
				distance_centre[Chaine][Residu] = ((distance_centre[Chaine][Residu]*(compteur_prot-1)) +
													sqrt(
													pow(centre_masse_prot['x']-parcourt['x'],2) +
													pow(centre_masse_prot['y']-parcourt['y'],2) +
													pow(centre_masse_prot['z']-parcourt['z'],2))
													) / compteur_prot
	return distance_centre
