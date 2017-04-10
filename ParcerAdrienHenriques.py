#!/usr/bin/env python
#-*- coding:utf8 -*-

""" Programme  de parcing pour lire un fichier PDB et en faire un fichier python
Fait le 20 fevrier """

def ParcerPDB(a):
	""" Ce programme prace un fichier au format PDB en format lisible par Python """
	contenu=list()
	mon_fichier=open(a,"r")
	for line in mon_fichier.readlines():
		contenu.append(line.strip())    #met le contenu du fichier pdb dans la liste "contenu"

	acidea=dict()

	for chain in range(len(contenu)): #On parcourt cette liste contenant tout le fichier pdb
		if contenu[chain][0:4]=="ATOM":   #Si la ligne commence par "ATOM" 
			Chaine = contenu[chain][21]
			
			if Chaine not in acidea.keys(): #Si la chaine ( A, B ... ) existe pas deja 
				acidea[Chaine] = {}     #creation du dictionnaire qui a pour nom les caractères a la ligne 21 ( Chaine)
			
			Posi = contenu[chain][24:26]
			if Posi not in acidea[Chaine].keys(): #Si la position pour une chaine n'existe pas deja (ex : -3 dans la chaine A)
				acidea[Chaine][Posi]={} # creation du dictionnaire poisition dans le dictionnaire chaine 
			
			residu = contenu[chain][12:16]
			if residu not in acidea[Chaine][Posi].keys(): #si le residu n'existe pas deja pour une chaine et une position donnée (ex : un CO de la chaine A a la position -3)
				acidea[Chaine][Posi][residu]= {}  #Creation du dictionnaire nom de l'atome, contenu dans le dictionnaire position lui meme contenu dans le dictionnaire chaine	
			
			#repartition de l'information dans le dictionnaire.
			acidea[Chaine][Posi][residu]["x"] = contenu[chain][32:38] #Mise des information de X dans le dictionnaire atome
			acidea[Chaine][Posi][residu]["y"] = contenu[chain][40:46] #Mise des information de Y dans le dictionnaire atome
			acidea[Chaine][Posi][residu]["z"] = contenu[chain][48:54] #Meme chose pour Z
			acidea[Chaine][Posi][residu]["Id"] = contenu[chain][9:11] #Meme chose pour Identifiant

	return( acidea)


def Affichage(resultats):
	cmpt=0
	for chain in resultats.keys():
		print(chain+"\n")
		for residu in resultats[chain].keys():
			print "  "+residu+","
			for atome in resultats[chain][residu].keys():
				if cmpt==0:
					print "\n"
					cmpt=1
				print " "+"	"+atome+",",
			
				for info in resultats[chain][residu][atome].keys():
					print " "+" "+" "+info+" : "+resultats[chain][residu][atome][info],
				print "\n"
	

# Utilisation du parcer (main)

if __name__ == "main":

	try:
		fichier = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
		ParcerPDB(fichier)

		dddd_result = ParcerPDB(fichier)

		Affichage(dddd_result)
	except:
		print("Ce fichier n'existe pas")



