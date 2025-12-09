## Projet Avancé Groupe 3 : Noa Juguera, Roman Begorre, Matteo Lafitte

* Les structures utilisés pour le projets sont dans types.h
* Les fonctions de lecture des fichiers TSP sont dans Tsp_io .c et .h.
* Pour les distance euclidienne et GEO, il y a toujours une petite différence entre ce qui est attendu par les tests python et le code C, à cause de l'arrondissement.
* Les fichiers tsp sont séparés du code source, ils se trouvent dans TestTsp, il faut précisé le chemin d'accés lors de l'exécution, exemple : ./tsp -f ../TestTsp/att10.tsp -m bf
* Le fichier sortie.txt correspond au fichier de sortie de l'option -o, en mode append, exemple : ./tsp -f ../TestTsp/att10.tsp -m bf -o ../sortie.txt
* Pour éxécuter le brute force sur des instances de dimension >= 12 il faut ajouter l'option -F.
* A cause de la gestion des interruptions du Ctrl_C, pour arrêter l'éxécution du programme il faut fermer le terminal, notamment pour des instances très longue à calculer.
* Pour utiiser le force brute avec une matrice de distances pré-calculé, il faut ajouter -M, cela réduit le temps processeur.
* Sur des grosses instances, les tests Python ne donnet pas la force brute car c'est beaucoup trop long.
* Nous avons un aussi un main_dbg pour "debug" qui contient nos affichages personnalisés avec l'option -X, (ne prend pas les mêmes options).
* Nous avons modifié les requirements pour régler nos problèmes d'installation dans l'environnement virtuel.
