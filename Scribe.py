#!/usr/bin/env python
#-*- coding:utf8 -*-

"""
Auteur : HENRIQUES Adrien et ARBES Hugo 
date : 28/04/2017
Programme Scribe : Permet d'ecrire les differents dictionnaires dans des fichiers de sortie.
"""

import os

def ScribeGlobal(dicoRMSDglob):
	""" ouvre un fichier de sortie et ecrit dedant les informations selectionnées a propos de la RMSD globale et angle de giration """
	##### A AJOUTER LE RAYON DE GIRATION DANS LES PARAMETRES ET PARCOURIR LE FICHIER POUR ECRIRE LES DONNEES 
	sortieGlobale= open("Out_Global.txt",'w') #Ouverture d'un fichier de sortie en mode ecriture
	
	donnees ="ConformationRef\tConformationTmps\tRMSD\tAngleGiration\n" #1er ligne qui sera dans le fichier de sortie
	
	for prot1 in sorted(dicoRMSDglob.keys()): #Parcourt des conformations du dictionnaire de RMSDglobale
		for prot2 in sorted(dicoRMSDglob[prot1].keys()): #Parcourt de proteine referente du dictionnaire de RMSD globale
			donnees += str(prot2)+"\t"+str(prot1)+"\t"+str(round(dicoRMSDglob[prot1][prot2],3))+"\n" #Ce qu'on va ecrire a chaque ligne
			
	
	sortieGlobale.write(donnees) #Ecriture dans le fichier de sortie 
	sortieGlobale.close() #Fermeture du fichier de sortie


def ScribeLocal(dicoRMSDloc):
	"""Ouvre un fichier de sortie et ecrit les informations retenue dans la RMSD locale et angle de giration """
	############# A AJOUTER LE RAYON DE GIRATION !!
	sortieLocale= open("Out_local.txt",'w')
	
	donneesloc ="Position\tRMSD\tAngleGiration\n"
	
	for prot in dicoRMSDloc.keys(): #On parcourt les differentes proteines de reference du dictionnaire de RMSDlocale
		for position in sorted(dicoRMSDloc[prot].keys()): #on parcourt son 2ieme niveau de dictionnaire correspondant aux positions differentes
	
			donneesloc += str(position)+"\t"+str(round(dicoRMSDloc[prot][position],3))+"\n"
			
	
	sortieLocale.write(donneesloc)
	sortieLocale.close()
