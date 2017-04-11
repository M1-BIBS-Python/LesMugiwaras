#!/usr/bin/env python
#-*- coding : utf8 -*-


"""
Author: Hugo ARBES
Contact: hugo.arbes@u-psud.fr
Date: 13/03/2017
Description: Script containing many useful functions for parsing and processing PDF files (i.e. 3D structure of proteins)
Licence (que je garde car elle me semble parfaite): DSSL (http://dssl.flyounet.net/)
"""



def Parcer_PDB(entree) :
	"""
	Fonction permettant de creer un dictionnaire definissant une proteine en prenant un fichier PDB en parametre

	Hugo ARBES le 06/03/2017 : Parser de fichier PDB
	"""
	
	try :
		pdb = open(entree,"r")
		lines = pdb.readlines()
		pdb.close()
		
	except :
		print "Le fichier ne s'est pas ouvert correctement"


	nbr_ligne = len(lines)
	dico = {}
	dico["Chaine"] = []
	
	for i in range(nbr_ligne):
		
		if(lines[i][0:4] == "ATOM"):
			if((lines[i][16] == "A") or (lines[i][16] == " ")):
				
				chain = lines[i][21:22]
				
				if not chain in dico["Chaine"]:
					dico["Chaine"].append(chain)
					dico[chain]= {}
					dico[chain]["Residus"] = []
						
				residu = lines[i][24:26]
				
				if not residu in dico[chain]["Residus"]:
					dico[chain]["Residus"].append(residu)
					dico[chain][residu]= {}
					dico[chain][residu]["Atome"] = []
				
				atom = lines[i][13:15]
				
				if not atom in dico[chain][residu]["Atome"]:
					dico[chain][residu]["Atome"].append(atom)
					dico[chain][residu][atom] = {}
				
				dico[chain][residu][atom]["Name"] = lines[i][13:15]
				dico[chain][residu][atom]["x"] = lines[i][31:38]
				dico[chain][residu][atom]["y"] = lines[i][39:46]
				dico[chain][residu][atom]["z"] = lines[i][47:54]
	
	return dico



def CentreMasseResidu(dictionnaire_residu):	# Mettre dico de 1 residu
	
	atomicMass = dict()
	atomicMass = {"C": float(12.0107), "N": float(14.0067), 
				  "O": float(15.9994), "S": float(32.065),
				  "OH": float(17.0073), "NH": float(15.0146)}
	
	x = 0
	y = 0
	z = 0
	mass_totale = 0
	cm_res = {}
	
	# Ajout de toutes les coordonnees des atomes du residu
	for atom_compare in dictionnaire_residu.keys():
			for atom_reference in atomicMass.keys():
				print(atom_reference)
				print(atom_compare + "\n")
				if (atom_reference == atom_compare) :
					print("FUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUCK")
					mass = atomicMass[atom_reference]
				
			x += mass*dictionnaire_residu[atom_compare]['x']
			y += mass*dictionnaire_residu[atom_compare]['y']
			z += mass*dictionnaire_residu[atom_compare]['z']
			mass_totale += mass
			
	# On divise par la masse totale des atomes pour avoir les coordonnees du centre de masse du residu:
	cm_res['x'] = x/mass_totale
	cm_res['y'] = y/mass_totale
	cm_res['z'] = z/mass_totale
	
	return cm_res






if __name__ == '__main__':
	
	#~ import json
	#~ print("Donnees pour l'arginine : \n"+json.dumps(ParcerPDB("1EJH.pdb"), indent = 4))
	
	#~ fichier = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
	
	dico = ParcerPDB("1EJH.pdb")
	print(dico)
	print("\n")
	print(dico.keys())
	print(dico["A"].keys())
	CentreDeMasse = CentreMasseResidu(dico["A"]["03"])



