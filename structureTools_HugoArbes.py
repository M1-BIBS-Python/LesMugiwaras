#!/usr/bin/env python
#-*- coding : utf8 -*-


from math import sqrt

from matplotlib import pyplot as plt




# Fonction prenant en argument le dictionnaire des conformations de la proteine et retournant son centre de masse
def CentreMasseProt(dico_proteine):
	
	nbr_residus = 0
	x = 0
	y = 0
	z = 0
	centre_masse_prot = {}
	
	# Ajout de toutes les coordonnees des centres de masse de tous les residus
	for Proteine in dico_proteine.keys():								# Pour chaque conformation de la proteine
		for Chaine in dico_proteine[Proteine].keys():					# Pour chaque chaine de cette conformation
			for Residu in dico_proteine[Proteine][Chaine].keys():		# Pourchaque residu de cette chaine
				
				# On ajoute les coordonnees du centre de masse
				x += CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])['x']
				y += CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])['y']
				z += CentreMasseResidu(dico_proteine[Proteine][Chaine][Residu])['z']
				nbr_residus += 1
				
	# On fait la moyenne des coordonnees
	centre_masse_prot['x'] = x / nbr_residus
	centre_masse_prot['y'] = y / nbr_residus
	centre_masse_prot['z'] = z / nbr_residus
	
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
				distant = sqrt(pow(CdM_prot['x']-parcourt['x'],2) + pow(CdM_prot['y']-parcourt['y'],2) + pow(CdM_prot['z']-parcourt['z'],2))
				
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


