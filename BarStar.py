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
	fichier1 = raw_input ("Saississez le nom de votre fichier de reference avec le format (ex: reference.pdb):")
	fichier2 = raw_input ("Saississez le nom de votre fichier avec le format (ex: prot_test1.pdb):")
	#~ monfichier1 = open("start_prot_only.pdb",'r')
	#~ print "AAAAAAAAAAAAAAAA"
	#~ fichier1 = monfichier.read()
	#~ print "UUUUUUUUUUUUUU"
	#Affichage(result_fichier1)
	#~ fichier2 = open("md_prot_only_skip100.pdb",'r')
	#Affichage(result_fichier1)
	#~ print "OOOOOOOOO"
	
except:
	print("Au moins l'un des deux fichiers specifies n'existe pas ! ")


try:
	result_fichier1 = ParserPDB(fichier1)
	result_fichier2 = ParserPDB(fichier2)
except:
	print("Le parsing des fichiers n'a pas fonctionne correctement")


try:
	print "AAAAAAAAAAAAAA1"
	# Appel des fonctions de calcul des centres de masse de la proteine de reference et des conformation + affichage pour chaque conformation
	#~ CentreDeMasseProt = CentreMasseProt(result_fichier2)
	rayon_proteine = RayonProt(result_fichier2)
	print "PPPPPPPPPPPPP1"
	GraphDesRayons(rayon_proteine)

	print "AAAAAAAAAAAAAA2"
	#Appel des fonction de calcul de RMSD global + affichage
	RmsdGlobal= RMSD(result_fichier1,result_fichier2)
	print "PPPPPPPPPPPPP2"
	graphglob(RmsdGlobal)

	print "QQQQQQQQQQQQ"
	ScribeGlobal(RmsdGlobal,rayon_proteine)

except:
	print("L'analyse globale n'a pas ete correctement effectuee")


try:
	print "BBBBBBBBBBBB1"
	# Appel des fonctions de calcul des distances de chaque residu par rapport au centre et affichage des ces distances
	dico_distances = DistanceCentre(result_fichier2)
	print "FFFFFFFFFFFF1"
	GraphDesDistances(dico_distances)

	print "BBBBBBBBBBBB2"
	#Appel des fonction de calcul de RMSD locale + affichage
	RmsdLocal= RMSDlocaux(result_fichier1,result_fichier2)
	print "FFFFFFFFFFFF2"
	graph(RmsdLocal)

	print "HHHHHHHHHHHH"
	ScribeLocal(RmsdLocal,dico_distances)

	print "CCCCCCCCCCCCCC"

except:
	print("L'analyse locale n'a pas ete correctement effectuee")
	

