#!/usr/bin/env python
#-*- coding:utf8 -*-

"""
Auteur : HENRIQUES Adrien et ARBES Hugo 
date : 28/04/2017
Programme BarStar : Appel l'ensemble des fonctions crées pour realiser le projet BarStar.
"""
########### IMPORT ################"
from Scribe import ScribeGlobal
from Scribe import ScribeLocal
from RMSDlocaux import RMSDlocaux
from RMSDlocaux import graph
from RMSDglobaux import RMSD
from RMSDglobaux import graphglob
from ParcerPDB import ParserPDB
from ParcerPDB import Affichage


########### MAIN ###############

try:
	fichier1 = raw_input ("Saississez le nom de votre fichier de reference avec le format (ex: arginine.pdb):")
	result_fichier1 = ParserPDB(fichier1)
	#Affichage(result_fichier1)
	fichier2 = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
	result_fichier2 = ParserPDB(fichier2)
	#Affichage(result_fichier1)
		
	#Appel des fonction de calcul de RMSD global + affichage
	RmsdGlobal= RMSD(result_fichier1,result_fichier2 )
	graphglob(RmsdGlobal)
	ScribeGlobal(RmsdGlobal) #########################################AJOUTER LE DICTIONNAIRE VOIR SCRIBE
	
	#Appel des fonction de calcul de RMSD locale + affichage
	RmsdLocal= RMSDlocaux(result_fichier1,result_fichier2 )
	graph(RmsdLocal)
	ScribeLocal(RmsdLocal) #########################################AJOUTER LE DICTIONNAIRE VOIR SCRIBE

	
except:
	print("L'un des deux fichier specifié n'existe pas ! ")

## Revoir la gestion d'erreur car elle m'affiche toujours ce message d'erreur meme si l'erreur survient apres l'entree des fichier 
