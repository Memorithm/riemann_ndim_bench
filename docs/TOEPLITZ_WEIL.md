# Phase 2 — discrétisation Toeplitz de la positivité de Weil

## Statut

Cette phase implémente **l'infrastructure finie** d'une discrétisation utilisée par Alain Connes et Caterina Consani dans *Weil positivity and Trace formula, the archimedean place*.

Elle ne calcule pas encore le noyau arithmético-spectral complet et ne constitue pas un test de l'hypothèse de Riemann.

Référence primaire :

- Alain Connes, Caterina Consani, *Weil positivity and Trace formula, the archimedean place*, arXiv:2006.13771, notamment §6.1 et les équations (105)–(106): https://arxiv.org/abs/2006.13771

## 1. Réseau logarithmique

Le papier remplace le groupe multiplicatif continu par un sous-groupe discret `q^Z` avec `q -> 1+` et pose

```text
omega = log(q).
```

Sur un intervalle logarithmique `[0,a]`, les points discrets sont

```text
x_j = j omega,
```

avec

```text
N = floor(a / omega),
j = 0,...,N.
```

Le code `LogLattice` reproduit exactement cette convention. La dimension de la matrice finie est donc `N + 1`.

Pour l'échelle utilisée dans le papier avec

```text
a = log(2),
omega = 10^-3,
```

on obtient

```text
N = 693,
dimension = 694.
```

## 2. Structure Toeplitz

Après discrétisation, le coefficient d'une entrée dépend uniquement de `|i-j|`. La matrice est donc Toeplitz et réelle symétrique dans la version considérée ici.

Nous définissons la fonction normalisée

```text
chi(x) = (Q epsilon)(exp(x)) / (2 epsilon'(1+)).
```

Le banc construit alors les coefficients finis

```text
T[i,j] = omega * chi(|i-j| omega).
```

C'est cette seule étape de discrétisation que `SymmetricToeplitz::sample_normalized_kernel` automatise.

## 3. Ce qui n'est PAS encore implémenté

La fonction `chi` n'est pas inventée par le banc.

Dans la référence primaire, `epsilon` est construite à partir de fonctions sphéroïdales prolates et l'analyse de `Q epsilon` fait intervenir cette construction. Par conséquent, la phase 2 exige que le noyau `chi` soit fourni explicitement par l'appelant.

Les noyaux simples utilisés dans les tests (`constant`, matrice 2x2, etc.) sont uniquement des **tests d'infrastructure informatique**. Ils ne représentent pas la fonction `Q epsilon` de Connes–Consani.

## 4. Pourquoi cette matrice est différente de l'ancien Gram

Une matrice de Gram `G = Phi Phi*` est positive semi-définie par construction. Elle ne peut donc pas servir à détecter une positivité non triviale.

Ici, `SymmetricToeplitz` n'impose aucune positivité. Le test unitaire avec première ligne

```text
[0, 1]
```

produit une forme quadratique positive dans la direction `(1,1)` et négative dans la direction `(1,-1)`. Le code vérifie ainsi que l'infrastructure ne fabrique pas artificiellement la propriété que nous voudrons éventuellement étudier.

## 5. Deux chemins de calcul

Le banc fournit deux représentations complémentaires :

- `apply`: application matricielle sans matérialiser les `N^2` coefficients ; seul le premier rang Toeplitz est stocké ;
- `dense` + diagonalisation auto-adjointe `faer`: chemin de validation pour dimensions finies raisonnables.

L'optimisation ARM/SVE viendra après validation scientifique du noyau.

## 6. Prochaine étape scientifique

Avant d'utiliser les valeurs propres à des fins de recherche, il faut reproduire indépendamment et fidèlement la construction de `epsilon`, puis de `Q epsilon`, à partir des fonctions sphéroïdales prolates décrites dans la référence primaire.

Jusqu'à cette étape, aucune valeur propre issue d'un noyau synthétique ne doit être interprétée comme information sur les zéros de `zeta(s)`.
