#!/usr/bin/env python
#-*- coding:utf8 -*-

"""
Auteur : HENRIQUES Adrien et ARBES Hugo 
date : 28/04/2017
Programme BarStar : Appel l'ensemble des fonctions crées pour realiser le projet BarStar.
"""
########### IMPORT ################"
from math import sqrt
from structureTools_HugoArbes import CentreMasseProt
from structureTools_HugoArbes import CentreMasseResidu
from structureTools_HugoArbes import RayonProt
from structureTools_HugoArbes import GraphDesRayons
from structureTools_HugoArbes import DistanceCentre
from structureTools_HugoArbes import GraphDesDistances

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
	# Ouverture des fichiers d'entree
	fichier1 = raw_input ("Saississez le nom de votre fichier de reference avec le format (ex: reference.pdb):")
	fichier2 = raw_input ("Saississez le nom de votre fichier avec le format (ex: prot_test1.pdb):")

except:
	print("Au moins l'un des deux fichiers specifiés n'existe pas ! ")


try:
	# Parsing des fichiers .pdb
	result_fichier1 = ParserPDB(fichier1)
	result_fichier2 = ParserPDB(fichier2)

except:
	print("Le parsing des fichiers n'a pas fonctionne correctement")


try:
	# Appel des fonctions de calcul des centres de masse de la proteine de reference et des conformation + affichage pour chaque conformation
	rayon_proteine = RayonProt(result_fichier2)
	GraphDesRayons(rayon_proteine)

	#Appel des fonction de calcul de RMSD global + affichage
	RmsdGlobal= RMSD(result_fichier1,result_fichier2)
	graphglob(RmsdGlobal)
	
	ScribeGlobal(RmsdGlobal,rayon_proteine)

except:
	print("L'analyse globale n'a pas ete correctement effectuee")


try:
	# Appel des fonctions de calcul des distances de chaque residu par rapport au centre et affichage des ces distances
	dico_distances = DistanceCentre(result_fichier2)
	GraphDesDistances(dico_distances)
	
	#Appel des fonction de calcul de RMSD locale + affichage
	RmsdLocal= RMSDlocaux(result_fichier1,result_fichier2)
	graph(RmsdLocal)
	
	ScribeLocal(RmsdLocal,dico_distances)

except:
	print("L'analyse locale n'a pas ete correctement effectuee")
	

