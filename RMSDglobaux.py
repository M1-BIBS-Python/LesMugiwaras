#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

""" Programme RMSD """
from math import sqrt
#### PARCER PDB #######
from TestParcer import ParserPDB
from TestParcer import Affichage

def RMSD(fichierref, fichier2):
	#Creation et initialisation des differentes listes dont on a besoin
	resultat=list()
	distance=0
	RmsdFinal=list()
	RmsdFichier=list()
	RmsdTot=0
	for prot1 in fichier2.keys(): #On parcours les differentes proteines contenues dans le fichier 2
		for prot2 in fichierref.keys(): #On parcours la proteine contenu dans le fichier de reference pour pouvoir la comparée au autres proteines.
			RmsdTot=0
			if prot2 not in prot1: #Si la proteine de reference est comprise dans le fichier 2.
				
				for chaine in fichier2[prot1].keys(): #Parcours des chaine
					RmsdFinal=list()
					RmsdPosition=0
					for posi in fichier2[prot1][chaine].keys(): #Parcours des positions 
						compteur=0
						resultat=list()
						for atome in fichier2[prot1][chaine][posi].keys(): #Parcours des Atomes dans le fichier 2
							distance=0
							for atome2 in fichierref[prot2][chaine][posi].keys(): #Parcours des Atomes dans le fichier 1
								if atome2==atome: #Si les atomes sont les memes dans les 2 fichiers 
									#Calcules :
									mesure=0
									mesure=sqrt((float(fichier2[prot1][chaine][posi][atome]["x"])-float(fichierref[prot2][chaine][posi][atome2]["x"]))**2+(float(fichier2[prot1][chaine][posi][atome]["y"])-float(fichierref[prot2][chaine][posi][atome2]["y"]))**2+(float(fichier2[prot1][chaine][posi][atome]["z"])-float(fichierref[prot2][chaine][posi][atome2]["z"]))**2)
									compteur+=1
									distance+= mesure**2 #Distance au carre
							resultat.append(distance) #Recuperation des distances au carre
						RmsdPosition=sqrt(sum(resultat)/compteur) #Somme des distances au carré divisee par le nombre de paires d'atomes rencontre (Formule poly)
						
					RmsdFinal.append(RmsdPosition) #Pour chaque positions on concatene les RMSD
					RmsdTot= sum(RmsdFinal) #On somme les RMSD de chaque positions
				RmsdFichier.append(str(prot1)+" "+str(prot2)+" "+str(RmsdTot)+" , ") #Creation de la liste de sortie 
	return (RmsdFichier)
	
	
	
	
	



fichier1 = raw_input ("Saississez le nom de votre fichier de reference avec le format (ex: arginine.pdb):")
result_fichier1 = ParserPDB(fichier1)
#Affichage(result_fichier1)

fichier2 = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
result_fichier2 = ParserPDB(fichier2)
#Affichage(result_fichier1)
calcul= RMSD(result_fichier1,result_fichier2 )
print calcul
