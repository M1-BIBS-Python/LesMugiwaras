#!/usr/bin/env python
#-*- coding : utf8 -*-

from math import sqrt
#### PARCER PDB #######

def ParserPDB(a):
	""" Ce programme prace un fichier au format PDB en format lisible par Python """
	contenu=list()
	mon_fichier=open(a,"r")
	for line in mon_fichier.readlines():
		contenu.append(line.strip())    #met le contenu du fichier pdb dans la liste "contenu"

	acidea=dict()
	


	for chain in range(len(contenu)): #On parcourt cette liste contenant tout le fichier pdb
		if contenu[chain][0:5]=="TITLE":
			newProt = contenu[chain][56:76]
			
			if newProt not in acidea.keys():
				acidea[newProt]={}
		
		
		if contenu[chain][0:4]=="ATOM":   #Si la ligne commence par "ATOM" 
			chaine = contenu[chain][21]
			
			if chaine not in acidea[newProt].keys(): 	#Si la chaine ( A, B, ... ) n'existe pas deja 
				acidea[newProt][chaine] = {}     		#creation du dictionnaire qui a pour nom les caracteres a la ligne 21 ( chaine)
			
			residu = contenu[chain][24:26]
			if residu not in acidea[newProt][chaine].keys(): #Si la residution pour une chaine n'existe pas deja (ex : -3 dans la chaine A)
				acidea[newProt][chaine][residu]={} # creation du dictionnaire poisition dans le dictionnaire chaine 
			
			atome = contenu[chain][13:16]
			if atome not in acidea[newProt][chaine][residu].keys(): #si le atome n'existe pas deja pour une chaine et une residution donnee (ex : un CO de la chaine A a la residution -3)
				acidea[newProt][chaine][residu][atome]= {}  #Creation du dictionnaire nom de l'atome, contenu dans le dictionnaire residution lui meme contenu dans le dictionnaire chaine	
			
			#repartition de l'information dans le dictionnaire.
			acidea[newProt][chaine][residu][atome]["x"] = contenu[chain][32:38] #Mise des information de X dans le dictionnaire atome
			acidea[newProt][chaine][residu][atome]["y"] = contenu[chain][40:46] #Mise des information de Y dans le dictionnaire atome
			acidea[newProt][chaine][residu][atome]["z"] = contenu[chain][48:54] #Meme chose pour Z
			acidea[newProt][chaine][residu][atome]["Id"] = contenu[chain][9:11] #Meme chose pour Identifiant

	return(acidea)


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
				if atom_reference is "H":	# On effectue ce test car il peut y avoir des H seuls ou NH ou OH
					masse_atome += 1.0079
				else:
					masse_atome = masse_molaire[atom_reference]
				
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
	
	for Proteine in dico_complexe_prot.keys():
		compteur_prot += 1
		#~ print ("prot :"+Proteine)
		#~ print dico_complexe_prot[Proteine].keys()
		for Chaine in dico_complexe_prot[Proteine].keys():
			#~ print Chaine
			if Chaine not in distance_centre.keys():
				#~ print "initialise"
				distance_centre[Chaine] = {}
			
			for Residu in dico_complexe_prot[Proteine][Chaine].keys():
				if Residu not in distance_centre[Chaine].keys():
					#~ print "initialise"
					distance_centre[Chaine][Residu] = 0
				
				parcourt = CentreMasseResidu(dico_complexe_prot[Proteine][Chaine][Residu])
				#~ print parcourt
				distance_centre[Chaine][Residu] += sqrt(pow(centre_masse_prot['x']-parcourt['x'],2) + pow(centre_masse_prot['y']-parcourt['y'],2) + pow(centre_masse_prot['z']-parcourt['z'],2)) / compteur_prot
		#~ print compteur_prot
	return distance_centre


if __name__ == '__main__':
	
	import json
	print(json.dumps(ParserPDB("ClientTest.pdb"), indent = 4))
	
	#~ fichier = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
	
	dico_complexe_prot = ParserPDB("ClientTest.pdb")
	
	CentreDeMasseProt = CentreMasseProt(dico_complexe_prot)
	print CentreDeMasseProt
	rayon_proteine = RayonProt(dico_complexe_prot)
	print(rayon_proteine)
	dico_distances = DistanceCentre(dico_complexe_prot)
	print dico_distances
