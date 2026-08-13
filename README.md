# riemann_ndim_bench

Banc de recherche numérique en Rust autour de l'hypothèse de Riemann, conçu pour séparer strictement :

1. les identités mathématiques exactes ;
2. les constructions numériques vérifiables ;
3. les hypothèses géométriques expérimentales.

## Statut scientifique

Ce dépôt **ne contient pas de preuve de l'hypothèse de Riemann**.

L'hypothèse de Riemann affirme que tous les zéros non triviaux de `zeta(s)` ont une partie réelle égale à `1/2`. Le Clay Mathematics Institute la classe toujours parmi les Millennium Prize Problems non résolus.

Référence : https://www.claymath.org/millennium/riemann-hypothesis/

## Fondation mathématique utilisée

La fonction xi de Riemann peut être écrite

```text
xi(s) = (1/2) s (s - 1) Gamma(s/2) pi^(-s/2) zeta(s)
```

et satisfait

```text
xi(s) = xi(1 - s).
```

Référence : NIST Digital Library of Mathematical Functions, §25.4, équations 25.4.3 et 25.4.4 :
https://dlmf.nist.gov/25.4

Pour `s = sigma + i t`, le module du seul facteur `pi^(-s/2)` vaut exactement

```text
|pi^(-s/2)| = pi^(-sigma/2).
```

Le dépôt introduit donc comme **coordonnée expérimentale**

```text
R_pi(sigma) = pi^(-sigma/2)
```

et sa version normalisée sur la ligne critique

```text
rho(sigma) = R_pi(sigma) / R_pi(1/2)
           = pi^((1 - 2 sigma)/4).
```

Cette définition donne exactement

```text
rho(1/2) = 1
rho(1 - sigma) = 1 / rho(sigma)
log rho(1 - sigma) = -log rho(sigma).
```

**Important :** interpréter `R_pi` ou `rho` comme un rayon géométrique ou physique est une hypothèse de travail du banc, pas un résultat connu sur les zéros de zeta.

## Pourquoi ne pas utiliser une simple matrice de Gram ?

Une matrice de la forme

```text
G_ij = integral phi_i(t) conjugate(phi_j(t)) dt
```

est positive semi-définie par construction, puisque pour tout vecteur `c`,

```text
c* G c = integral |sum_i c_i phi_i(t)|^2 dt >= 0.
```

La positivité d'une telle matrice ne peut donc pas, à elle seule, distinguer `sigma = 1/2` de `sigma != 1/2`. Le banc évitera de confondre un échec numérique de Cholesky avec une propriété de la fonction zeta.

## Direction N-dimensionnelle

La cible du projet est une suite de modèles finis `E_N` dont la dimension augmente, avec une forme quadratique ou un opérateur dont la positivité n'est **pas automatique**.

Une direction documentée est la positivité de Weil et ses formulations spectrales. Connes et Consani étudient précisément une interprétation hilbertienne de cette positivité :

- Alain Connes, Caterina Consani, *Weil positivity and Trace formula, the archimedean place* (2020): https://arxiv.org/abs/2006.13771
- Alain Connes, *The Riemann Hypothesis: Past, Present and a Letter Through Time* (2026): https://arxiv.org/abs/2602.04022

Ces références orientent la phase N-dimensionnelle ; elles ne valident pas notre interprétation radiale.

## Feuille de route

### Phase 0 — invariants élémentaires

- [x] représentation de `s = sigma + i t`
- [x] réflexion fonctionnelle `s -> 1 - s`
- [x] réflexion géométrique autour de `Re(s)=1/2`
- [x] coordonnée radiale `R_pi`
- [x] coordonnée normalisée `rho`
- [x] tests d'involution et de réciprocité

### Phase 1 — évaluation numérique contrôlée

- [ ] implémentation de `zeta(s)` avec plusieurs régimes numériques
- [ ] implémentation de `xi(s)`
- [ ] tests de l'équation fonctionnelle
- [ ] contrôle d'erreur et précision arbitraire pour les validations sensibles
- [ ] comparaison à des valeurs de référence indépendantes

### Phase 2 — espace N-dimensionnel

- [ ] définir explicitement la base de fonctions tests `E_N`
- [ ] construire une forme quadratique non trivialement positive
- [ ] calculer valeurs propres / LDL^T sans inverser explicitement les matrices
- [ ] suivre la plus petite valeur propre quand `N` augmente
- [ ] tester la stabilité par changement de base et de précision

### Phase 3 — géométrie radiale expérimentale

- [ ] comparer la coordonnée `rho` aux observables de la forme quadratique
- [ ] tester si elle apporte une information indépendante ou seulement un changement de variable
- [ ] rejeter le modèle radial s'il ne produit aucun invariant nouveau

### Phase 4 — calcul haute performance

Cible principale : machine ARM64 `aarch64`, 14 CPU, 122 GiB RAM, SIMD Advanced SIMD + SVE/SVE2.

L'optimisation ne sera entreprise qu'après validation mathématique du noyau numérique.

## Principe du projet

Un phénomène numérique n'est jamais considéré comme une preuve. Toute signature observée doit survivre :

- à l'augmentation de `N` ;
- à l'augmentation de la précision ;
- au changement de discrétisation ;
- au changement de base ;
- à une implémentation ou référence indépendante.
