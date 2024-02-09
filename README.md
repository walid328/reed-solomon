# reed-solomon

Une implémentation des codes de Reed Solomon sur F_p.

Pour compiler le projet on utilise make. 
Le fichier main définit une interface en ligne de commande pour encoder et décoder des codes de Reed Solomon avec les méthodes rapides.
Des tests sont fournis, on peut les lancer avec make test. On peut aussi faire des tests de performances avec make perf. Cela permet d'écrire dans un fichier .csv les temps de calculs.