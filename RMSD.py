#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

""" Programme RMSD """
from math import sqrt
#### PARCER PDB #######
from ParcerAdrienHenriques import ParcerPDB
from ParcerAdrienHenriques import Affichage

def RMSD(fichier1, fichier2):
	resultat=list()
	distance=0
	for residu in fichier1.keys():
		for posi in fichier1[residu].keys():
					compteur=0
					mesure=0
					for atome in fichier1[residu][posi].keys():
						for atome2 in fichier2[residu][posi].keys():
							if atome2==atome:
								mesure=sqrt((float(fichier1[residu][posi][atome]["x"])-float(fichier2[residu][posi][atome2]["x"]))**2+(float(fichier1[residu][posi][atome]["y"])-float(fichier2[residu][posi][atome2]["y"]))**2+(float(fichier1[residu][posi][atome]["z"])-float(fichier2[residu][posi][atome2]["z"]))**2)
								compteur+=1
								distance+= mesure**2
		distance=sqrt(distance/compteur)
		resultat.append(distance)
	return (resultat)
	
	
	
	
	



fichier1 = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
result_fichier1 = ParcerPDB(fichier1)
#Affichage(result_fichier1)

fichier2 = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
result_fichier2 = ParcerPDB(fichier2)
#Affichage(result_fichier2)

calcul= RMSD(result_fichier1, result_fichier2)
print calcul
