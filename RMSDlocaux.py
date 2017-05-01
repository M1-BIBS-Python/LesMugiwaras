#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
Auteur : HENRIQUES Adrien et ARBES Hugo 
date : 28/04/2017
Programme RMSDlocaux : calcule la RMSD locale de chaques conformation et 
fait la moyenne de ces RMSD a chaques positions pour avoir une RMSD moyenne local 
Comprend aussi une fonction graphique pour observer la RMSD moyenne locale graphiquement
et pouvoir determiner les regions flexibles. """


############### IMPORT #######################
from math import sqrt
#### PARCER PDB #######
from ParcerPDB import ParserPDB
from ParcerPDB import Affichage

#### librairie graphique ####
from matplotlib import pyplot as plt
import numpy as np

################## fonction RMSDlocaux ##########################

def RMSDlocaux(fichierref, fichier2):
	""" Cette fonction calcul la RMSD locale moyenne au differentes positions """
	#Creation et initialisation des differentes listes dont on a besoin
	resultat=list()
	distance=0
	RmsdFinal=list()
	RmsdTot=0
	RMSDloc=dict()
	for prot1 in fichier2.keys(): #On parcours les differentes proteines contenues dans le fichier 2
		if prot1 not in RMSDloc.keys():
			RMSDloc[prot1]={}
		for prot2 in fichierref.keys(): #On parcours la proteine contenu dans le fichier de reference pour pouvoir la comparée au autres proteines.
			RmsdTot=0
			for chaine in fichier2[prot1].keys(): #Parcours des chaine
				if chaine not in RMSDloc[prot1].keys():
					RMSDloc[prot1][chaine]={}
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
					RMSDloc[prot1][chaine][posi]= RmsdPosition

#BOUCLE POUR ADDITIONNER LES VALEURS DE CHAQUES POSITIONS POUR POUVOIR FAIRE UNE MOYENNE
	RMSDmoy=dict()
	cmpt=0
	for prot_dict in RMSDloc.keys(): #Parcourt des conformations non referentes contenue dans le dictionnaire obtenu
		cmpt+=1						 #Compteur servant a compter le nombre de conformations (pour pouvoir faire la moyenne)
		for chain in RMSDloc[prot_dict].keys():  #Parcourt des differentes chaines dans le dictionnaire
			if chain not in RMSDmoy.keys():
				RMSDmoy[chain]={}
			for residu in RMSDloc[prot_dict][chain].keys(): #Parcourt des residus dans le dictionnaire
				if residu not in RMSDmoy[chain].keys():
					RMSDmoy[chain][residu]=0
				RMSDmoy[chain][residu]=RMSDmoy[chain][residu]+RMSDloc[prot_dict][chain][residu] #Somme des RMSD locales pour chaques conformations que l'on met dans un dictionnaire
	
#BOUCLE POUR DIVISER A CHAQUES POSITIONS POUR AVOIR LA MOYENNE
	for chain in RMSDmoy.keys():  #On parcourt le dictionnaire de nouveau obtenu
		for residu in RMSDmoy[chain].keys(): #on parcourt les residu du dictionnaire
			RMSDmoy[chain][residu]=RMSDmoy[chain][residu]/cmpt #on divise les sommes des RMSDlocales par le nombre de conformation pour avoir la moyenne
		

	return (RMSDmoy)

############## Fonction Graphique ######################

def graph (dico):
	""" cette fonction trace le graphique des RMSD moyenne locale en fonction de la position """
	
	x= list()
	y=list()
	for chaine in dico.keys(): #On parcourt le dictionnaire mit en paramètre
		for position in sorted(dico[chaine].keys()): #on parcourt son 2ieme niveau de dictionnaire correspondant aux positions differentes
			x.append(position) #On concatene la position dans la liste X
			y.append(dico[chaine][position]) #On concatene les valeurs des positions dans la liste Y
	
	plt.plot(x, y) #Graphique des RMSDlocales moyenne en fonction des positions
	plt.title("Regions flexibles")
	plt.xlabel('Position')
	plt.ylabel('RMSD')
	plt.show()
	
	return ()	


################## MAIN ###########################

if __name__ == 'main':
	try:
		fichier1 = raw_input ("Saississez le nom de votre fichier de reference avec le format (ex: arginine.pdb):")
		fichier1 = open("start_prot_only.pdb")
		result_fichier1 = ParserPDB(fichier1)
		#Affichage(result_fichier1)
	except:
		print("Ce fichier n'existe pas")


	try:
		fichier2 = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
		fichier2 = open("_md_prot_only_skip100.pdb")
		result_fichier2 = ParserPDB(fichier2)
		#Affichage(result_fichier1)
	except:
		print("Ce fichier n'existe pas")

	calcul= RMSDlocaux(result_fichier1,result_fichier2 )
	print calcul

	graph(calcul)
