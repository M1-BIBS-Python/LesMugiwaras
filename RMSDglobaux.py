#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Auteur : HENRIQUES Adrien et ARBES Hugo 
date : 28/04/2017
Programme RMSDglobaux : calcule la RMSD de chaques conformation .
Comprend aussi une fonction graphique pour observer la RMSD graphiquement. """

####################### IMPORT #####################
from math import sqrt
#### PARCER PDB #######
from ParcerPDB import ParserPDB
from ParcerPDB import Affichage

#### librairie graphique ####
from matplotlib import pyplot as plt
import numpy as np

################# Fonction de calcul de RMSD #######################
def RMSD(fichierref, fichier2):
	""" Cette fonction calcul la RMSD globale des differentes conformations """
	
	#Creation et initialisation des differentes listes dont on a besoin
	resultat=list()
	distance=0
	RmsdFinal=list()
	RMSDglob=dict()
	RmsdTot=0
	for prot1 in fichier2.keys(): #On parcours les differentes proteines contenues dans le fichier 2
		if prot1 not in RMSDglob.keys():
			RMSDglob[prot1]={}
		
		for prot2 in fichierref.keys(): #On parcours la proteine contenu dans le fichier de reference pour pouvoir la comparée au autres proteines.

			RmsdTot=0
			if prot2 not in RMSDglob[prot1].keys():
				RMSDglob[prot1][prot2]={}
				
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
			RMSDglob[prot1][prot2]=RmsdTot
	return (RMSDglob)


##################### Fonction graphique ####################
def graphglob (dico):
	""" Cette fonction creer un graphique a partir du dictionnaire de sortie de la fonction RMSD """
	
	x= list()
	y=list()
	for prot1 in sorted(dico.keys()): #parcourt les conformations dans le dictionnaire de sortie de la fonction RMSD
		for prot2 in sorted(dico[prot1].keys()):#Parcourt les conformations (2ieme clé du dictionnaire)
			x.append(prot1)			#concatenation des identifiants des conformations rencontrées
			y.append(dico[prot1][prot2]) # Concaténation des valeurs de RMSD globaux dans Y 


	plt.plot(x, y)  #trace le graphique des valeurs de RMSD en fonction des conformations
	plt.title("Graphique des RMSD globaux")
	plt.xlabel('Conformation')
	plt.ylabel('RMSD globales')
	plt.show()
	
	return ()	
	
	
################### MAIN ###########################
if __name__ == "main":
	try:
		fichier1 = raw_input ("Saississez le nom de votre fichier de reference avec le format (ex: arginine.pdb):")
		result_fichier1 = ParserPDB(fichier1)
		#Affichage(result_fichier1)
	except:
		print("Ce fichier n'existe pas")

	
	try:
		fichier2 = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
		result_fichier2 = ParserPDB(fichier2)
		#Affichage(result_fichier1)
	except:
		print("Ce fichier n'existe pas")

	calcul= RMSD(result_fichier1,result_fichier2 )
	print calcul
	graphglob(calcul)
	

