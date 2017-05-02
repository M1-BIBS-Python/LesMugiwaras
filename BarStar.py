#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Auteur : HENRIQUES Adrien et ARBES Hugo 
date : 28/04/2017
Programme BarStar : Appel l'ensemble des fonctions crees pour realiser le projet BarStar.
"""
########### IMPORT ################

from RayonGiration import RayonProt
from RayonGiration import GraphDesRayons
from DistancesCentre import DistanceCentre
from DistancesCentre import GraphDesDistances
from DistancesCentre import GraphEnfouissementRMSD

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
	try:
		fichier1 = raw_input ("Saississez le nom de votre fichier de reference avec le format (ex: start_prot_only.pdb):")
		result_fichier1 = ParserPDB(fichier1)
	
	except:
		print("Ce fichier n'existe pas!")
	
	try:
		fichier2 = raw_input ("Saississez le nom de votre fichier avec le format (ex: md_prot_only_skip100.pdb):")
		result_fichier2 = ParserPDB(fichier2)
	except:
		print("Ce fichier n'existe pas !")
	
	
except:
	print("Le parsing des fichiers n'a pas pu s'effectuer correctement")


try:
	# Appel des fonctions de calcul des centres de masse de la proteine de reference et des conformation + affichage pour chaque conformation
	rayon_proteine = RayonProt(result_fichier2)
	GraphDesRayons(rayon_proteine)

	#Appel des fonction de calcul de RMSD global + affichage
	RmsdGlobal= RMSD(result_fichier1,result_fichier2)
	graphglob(RmsdGlobal)

	ScribeGlobal(RmsdGlobal,rayon_proteine)

except:
	print("L'analyse globale n'a pas pu s'effectuer correctement")


try:
	# Appel des fonctions de calcul des distances de chaque residu par rapport au centre et affichage des ces distances
	dico_distances = DistanceCentre(result_fichier2)
	GraphDesDistances(dico_distances)

	#Appel des fonction de calcul de RMSD locale + affichage
	RmsdLocal= RMSDlocaux(result_fichier1,result_fichier2)
	graph(RmsdLocal)
	
	# Appel de la fonction d'affichage des RMSD en fonction de l'enfouissement
	GraphEnfouissementRMSD(RmsdLocal,dico_distances)
	
	# Appel de la fonction d'ecriture dans le fichier de sortie
	ScribeLocal(RmsdLocal,dico_distances)

except:
	print("L'analyse locale n'a pas pu s'effectuer correctement")
	

