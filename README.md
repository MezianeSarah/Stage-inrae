# 🌲 Simulation numérique de la fracture du bois cellulaire (Travail en cours)

Ce dépôt contient les scripts Python développés dans le cadre de mon stage de Master 2 à l’Université de Montpellier.  
Je travaille sur la modélisation numérique de la rupture du bois à l’échelle cellulaire, en combinant :

- la **méthode pérydinamique** (modélisation de la fissuration),
- la **génération de microstructures** (diagrammes de Voronoï),
- et le **traitement d’images réelles de bois** (microscopie).

##  Objectif

Simuler la propagation de fissures dans des structures inspirées du bois naturel,  
et **comparer les résultats avec des images réelles** pour valider les comportements observés.

##  Ce que je fais actuellement

- Génération de cellules :
  -  Voronoï hexagonaux (réguliers, type nid d’abeille)
  -  Voronoï aléatoires (non hexagonaux, plus proches du bois réel)
- Simulation de la rupture :
  - Méthode **bond-based péridynamique**
  - Rupture automatique des liens selon un seuil d’élongation
- Traitement d’images réelles :
  - Extraction des contours cellulaires à partir de photos microscopiques
  - Comparaison avec les structures générées (forme, taille, distribution)

##  Outils et bibliothèques

- Python (NumPy, Matplotlib)
- Voronoi via Scipy ou CGAL
- Traitement d’image : OpenCV, Scikit-image
- Méthode péridynamique avec intégration explicite (Velocity-Verlet)

##  En cours

- Post-traitement et visualisation des champs de contraintes
- Validation par comparaison images vs modèles
- Amélioration des structures aléatoires à partir des images réelles

##  Analyse mécanique

À la fin des simulations, j’extrais le **module de Young effectif** du matériau modélisé en analysant la courbe contrainte-déformation.  
Ce module est ensuite **comparé au modèle théorique de Gibson-Ashby**, qui décrit les propriétés mécaniques des matériaux cellulaires.

###  Objectif de la comparaison :

- Vérifier si les **modèles Voronoï générés** (hexagonaux ou aléatoires) se comportent comme prévu.
- Évaluer l’influence de la **densité**, **forme des cellules**, et **type de maillage** sur les propriétés élastiques.
- Valider l’approche numérique par rapport à un **cadre théorique reconnu**.


