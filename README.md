# Friction-code-CD
The firction code writted in python.

Hugo Martin
version du 16 mai 2017



PREREQUIS:
- avoir Python version 3.4 avec les bonnes librairies
- avoir Mosek
- avoir créé dans le dossier des fichiers du code, les répertoires "temporary_files" et "outputs".



Rappel perso
Pour l activer sur mon pc perso il faut lancer l environnement py34 via la commande (mieux vaut faire un export après):
> source activate py34
> export PS1=""




Ce code fonctionne sur python 3.4
Il faut installer anaconda avec les librairies numpy, scipy et matplotlib.pyplot
Il faut installer la librairie Mosek pour anaconda version python 3.4 (utiliser conda) voir le site de Mosek

A noter la présence du fichier matrices.c indispensable à travers le fichier matrices.so
Pour modifier ce fichier .c utiliser le makefile avec l'option "make clean" pour écraser le fichier .so et 
"make all" pour recompiler le fichier .c

On notera immédiatement que dans cette version on a le principe que l'on boucle la
pluspart du temps sur des problèmes mal posés.

Cette version fontionne sur les principaux cas tests. Dans la version suivante,
j'aimerais que l'on puisse utiliser le principe de l'ajout d'une variable ce qui
permet d'avoir une fonctionnelle convexe immédiatement.
