# LesMugiwaras
Projet BARstar de ARBES HUGO et HENRIQUES ADRIEN

# INTRODUCTION 

Le groupe LesMugiwaras fut créé afin de resoudre un problème de docking posé durant notre cursus scolaire a l'université Paris-Sud (M1 BIBS). En effet la plus part des algorithmes de dockings considèrent que les proteines sont rigides et ne suibissent pas de changements conformationnels. Or en realité, certaines proteines telle que la proteine Barstar suibit des changements conformationnels lors de l'interaction avec la Barnase, formant ainsi le complexe Barstar-Barnase.
Le but de ce projet est donc de créer un algorithme qui va nous permettre d’étudier la stabilité de la structure de la protéine Barstar au cours du temps ainsi que ses changements conformationnels. 

# Exemple De Code 

Pour executer le code, il suffit d'executer le programme nommer BarStar.py.

```
  fichier1 = raw_input ("Saississez le nom de votre fichier de reference avec le format (ex: arginine.pdb):")
	result_fichier1 = ParserPDB(fichier1)
	#Affichage(result_fichier1)
	fichier2 = raw_input ("Saississez le nom de votre fichier avec le format (ex: arginine.pdb):")
	result_fichier2 = ParserPDB(fichier2)
	#Affichage(result_fichier1)
		
	#Appel des fonction de calcul de RMSD global + affichage
	RmsdGlobal= RMSD(result_fichier1,result_fichier2 )
	graphglob(RmsdGlobal)
	ScribeGlobal(RmsdGlobal)
	
	#Appel des fonction de calcul de RMSD locale + affichage
	RmsdLocal= RMSDlocaux(result_fichier1,result_fichier2 )
	graph(RmsdLocal)
	ScribeLocal(RmsdLocal)
```

# API 
Le projet est composé de plusieurs fichiers qui sont des fonctions du programme. Celui-ci comprend :
- ParcerPDB : contient la fonction qui sert a parser le fichier
- RMSDglobaux : contient la fonction de calcule et d'affichage graphique de la RMSD globale
- RMSDlocaux : contient la fonction de calcule et d'affchage graphique de la RMSD locale 
- Scribe : contient les fonctions d'ecriture des données générées dans des fichiers de sorties



# License
