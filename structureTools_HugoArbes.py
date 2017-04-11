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
		
		
		else:							# POUR LES TEEEEEEEEESTS MAIS A ENLEEEEEVEEEEEEEEEEEEER
			newProt="1"
			if newProt not in acidea.keys():
				acidea[newProt]={}
				
		
		
		if contenu[chain][0:4]=="ATOM":   #Si la ligne commence par "ATOM" 
			chaine = contenu[chain][21]
			
			if chaine not in acidea[newProt].keys(): #Si la chaine ( A, B ... ) existe pas deja 
				acidea[newProt][chaine] = {}     #creation du dictionnaire qui a pour nom les caracteres a la ligne 21 ( chaine)
			
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


def CentreMasseProt(dictionnaire_proteine):
	
	nbr_residus = 0
	x = 0
	y = 0
	z = 0
	centre_masse_prot = {}
	
	# Ajout de toutes les coordonnees des centres de masse de tous les residus
	for Proteine in dictionnaire_proteine.keys():
		for Chaine in dictionnaire_proteine[Proteine].keys():
			for Residu in dictionnaire_proteine[Proteine][Chaine].keys():
				x += CentreMasseResidu(dictionnaire_proteine[Proteine][Chaine][Residu])['x']
				y += CentreMasseResidu(dictionnaire_proteine[Proteine][Chaine][Residu])['y']
				z += CentreMasseResidu(dictionnaire_proteine[Proteine][Chaine][Residu])['z']
				nbr_residus += 1
				
	centre_masse_prot['x'] = x / nbr_residus
	centre_masse_prot['y'] = y / nbr_residus
	centre_masse_prot['z'] = z / nbr_residus
	
	return centre_masse_prot

def CentreMasseResidu(dictionnaire_residu):
	
	atomicMass = dict()
	atomicMass = {"C": float(12.0107), "S": float(32.065), 
				  "O ": float(15.9994), "N ": float(14.0067),
				  "OH": float(17.0073), "NH": float(15.0146)}
	
	x = 0
	y = 0
	z = 0
	masse_totale_residu = 0
	centre_masse_residu = {}
	
	# Ajout de toutes les coordonnees des atomes du residu
	for atom_compare in dictionnaire_residu.keys():
			for atom_reference in atomicMass.keys():
				if atom_reference in atom_compare:
					masse_atome = atomicMass[atom_reference]
				
			x += masse_atome * float(dictionnaire_residu[atom_compare]['x'])
			y += masse_atome * float(dictionnaire_residu[atom_compare]['y'])
			z += masse_atome * float(dictionnaire_residu[atom_compare]['z'])
			masse_totale_residu += masse_atome
			
	# On divise par la masse totale des atomes pour avoir les coordonnees du centre de masse du residu
	centre_masse_residu['x'] = x / masse_totale_residu
	centre_masse_residu['y'] = y / masse_totale_residu
	centre_masse_residu['z'] = z / masse_totale_residu
	
	return centre_masse_residu



def RayonProt(dictionnaire_proteine):
	
	CdM_prot = CentreMasseProt(dictionnaire_proteine)
	le_plus_loin = 0
	
	for Proteine in dictionnaire_proteine.keys():
		for Chaine in dictionnaire_proteine[Proteine].keys():
			for Residu in dictionnaire_proteine[Proteine][Chaine].keys():
				parcourt = CentreMasseResidu(dictionnaire_proteine[Proteine][Chaine][Residu])
				distant = sqrt(pow(CdM_prot['x']-parcourt['x'],2) + pow(CdM_prot['y']-parcourt['y'],2) + pow(CdM_prot['z']-parcourt['z'],2))
				if distant > le_plus_loin :
					le_plus_loin = distant
	
	return le_plus_loin



if __name__ == '__main__':
	
	import json
	print(json.dumps(ParserPDB("1EJH.pdb"), indent = 4))
	
	#~ fichier = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
	
	dico = ParserPDB("1EJH.pdb")
	CentreDeMasseProt = CentreMasseProt(dico)
	print CentreDeMasseProt
	rayon_proteine = RayonProt(dico)
	print(rayon_proteine)
