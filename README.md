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
		
	#Calculs globaux 
	rayon_proteine = RayonProt(result_fichier2)
	GraphDesRayons(rayon_proteine)
	RmsdGlobal= RMSD(result_fichier1,result_fichier2 )
	graphglob(RmsdGlobal)
	ScribeGlobal(RmsdGlobal)
	
	#Calculs locaux
	dico_distances = DistanceCentre(result_fichier2)
	GraphDesDistances(dico_distances)
	RmsdLocal= RMSDlocaux(result_fichier1,result_fichier2 )
	graph(RmsdLocal)
	GraphEnfouissementRMSD(RmsdLocal,dico_distances)
	ScribeLocal(RmsdLocal)
	
```

# API 
Le projet est composé de plusieurs fichiers qui sont des fonctions du programme. Celui-ci comprend :
- ParcerPDB : contient la fonction qui sert a parser le fichier
- RMSDglobaux : contient la fonction de calcule et d'affichage graphique de la RMSD globale
- RMSDlocaux : contient la fonction de calcule et d'affchage graphique de la RMSD locale 
- Scribe : contient les fonctions d'ecriture des données générées dans des fichiers de sorties
- CentresDeMasse : contient les fonctions de calcul de centre de masse de chaque conformation et pour un residu 
- RayonGiration : contient les fonctions de calcul de rayon de giration de chaque conformations et contient le fonction pour tracer le graph
- DistancesCentre : contient les fonction de calcul de distance entre le centre de masse d'un residu et le centre de masse de la proteine pour chaque conformation et une fonction permettant de tracer le graph de la distance en fonction du residu ainsi que la fonction permettant de tracer le graph de la RMSD moyenne en fonction de l'enfouissement des residus.

# Contributeurs 
LesMugiwaras : ARBES HUGO et HENRIQUES ADRIEN

# License
Le projet Barstar de la team LesMugiwaras est libre de droits 
