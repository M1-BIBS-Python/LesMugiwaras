#!/usr/bin/env python
#-*- coding:utf8 -*-

""" Programme  de parcing pour lire un fichier PDB et en faire un fichier python
Fait le 20 fevrier """

def ParserPDB(a):
	""" Ce programme prace un fichier au format PDB en format lisible par Python """
	contenu=list()
	mon_fichier=open(a,"r")
	for line in mon_fichier.readlines():
		contenu.append(line.strip())    #met le contenu du fichier pdb dans la liste "contenu"

	acidea=dict()
	


	for chain in range(len(contenu)): #On parcourt cette liste contenant tout le fichier pdb
		if contenu[chain][0:5]=="TITLE":
			newProt = contenu[chain][56:76]
			
			if newProt not in acidea.keys():
				acidea[newProt]={}
				
		if contenu[chain][0:4]=="ATOM":   #Si la ligne commence par "ATOM" 
			Chaine = contenu[chain][21]
			
			if Chaine not in acidea[newProt].keys(): #Si la chaine ( A, B ... ) existe pas deja 
				acidea[newProt][Chaine] = {}     #creation du dictionnaire qui a pour nom les caractères a la ligne 21 ( Chaine)
			
			Posi = contenu[chain][24:26]
			if Posi not in acidea[newProt][Chaine].keys(): #Si la position pour une chaine n'existe pas deja (ex : -3 dans la chaine A)
				acidea[newProt][Chaine][Posi]={} # creation du dictionnaire poisition dans le dictionnaire chaine 
			
			residu = contenu[chain][12:16]
			if residu not in acidea[newProt][Chaine][Posi].keys(): #si le residu n'existe pas deja pour une chaine et une position donnée (ex : un CO de la chaine A a la position -3)
				acidea[newProt][Chaine][Posi][residu]= {}  #Creation du dictionnaire nom de l'atome, contenu dans le dictionnaire position lui meme contenu dans le dictionnaire chaine	
			
			#repartition de l'information dans le dictionnaire.
			acidea[newProt][Chaine][Posi][residu]["x"] = float(contenu[chain][32:38]) #Mise des information de X dans le dictionnaire atome
			acidea[newProt][Chaine][Posi][residu]["y"] = float(contenu[chain][40:46]) #Mise des information de Y dans le dictionnaire atome
			acidea[newProt][Chaine][Posi][residu]["z"] = float(contenu[chain][48:54]) #Meme chose pour Z
			acidea[newProt][Chaine][Posi][residu]["Id"] = contenu[chain][9:11] #Meme chose pour Identifiant

	return( acidea)


def Affichage(resultats):
	cmpt=0
	for prot in resultats.keys():
		print(prot+"\n")
		for chain in resultats[prot].keys():
			print(chain+"\n")
			for residu in resultats[prot][chain].keys():
				print "  "+residu+","
				for atome in resultats[prot][chain][residu].keys():
					if cmpt==0:
						print "\n"
						cmpt=1
					print " "+"	"+atome+",",
				
					for info in resultats[prot][chain][residu][atome].keys():
						print " "+" "+" "+info+" : "+resultats[prot][chain][residu][atome][info],
				print "\n"
	

# Utilisation du parcer (main)


if __name__ == "main":

	try:
		fichier = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
		ParserPDB(fichier)

		dddd_result = ParserPDB(fichier)

		Affichage(dddd_result)
	except:
		print("Ce fichier n'existe pas")



