#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

""" Programme RMSD """
from math import sqrt
#### PARCER PDB #######
from TestParcer import ParserPDB
from TestParcer import Affichage

def RMSD(fichier):
	resultat=list()
	distance=0
	RmsdFinal=list()
	RmsdFichier=list()
	RmsdTot=0
	for prot1 in fichier.keys():
		for prot2 in fichier.keys():
			RmsdTot=0
			if prot2 not in prot1:
				
				for residu in fichier[prot1].keys():
					RmsdFinal=list()
					RmsdPosition=0
					for posi in fichier[prot1][residu].keys():
						compteur=0
						resultat=list()
						for atome in fichier[prot1][residu][posi].keys():
							distance=0
							for atome2 in fichier[prot2][residu][posi].keys():
								if atome2==atome:
									mesure=0
									mesure=sqrt((float(fichier[prot1][residu][posi][atome]["x"])-float(fichier[prot2][residu][posi][atome2]["x"]))**2+(float(fichier[prot1][residu][posi][atome]["y"])-float(fichier[prot2][residu][posi][atome2]["y"]))**2+(float(fichier[prot1][residu][posi][atome]["z"])-float(fichier[prot2][residu][posi][atome2]["z"]))**2)
									compteur+=1
									distance+= mesure**2
							resultat.append(distance)
						RmsdPosition=sqrt(sum(resultat)/compteur)
						
					RmsdFinal.append(RmsdPosition)
					RmsdTot= sum(RmsdFinal)
					print (RmsdTot)
				RmsdFichier.append(str(prot1)+" "+str(prot2)+" "+str(RmsdTot)+" , ")
	return (RmsdFichier)
	
	
	
	
	



fichier1 = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
result_fichier1 = ParserPDB(fichier1)
#Affichage(result_fichier1)

calcul= RMSD(result_fichier1)
print calcul
