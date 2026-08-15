# riemann_ndim_bench

Banc de recherche numérique en Rust autour de l'hypothèse de Riemann, conçu pour séparer strictement :

1. les identités mathématiques exactes ;
2. les constructions numériques vérifiables ;
3. les hypothèses et asymptotiques expérimentales.

## Statut scientifique

Ce dépôt **ne contient pas de preuve de l'hypothèse de Riemann**.

Les calculs actuels portent sur des compressions finies d'opérateurs de type prolate/semilocal inspirés des travaux de Connes, Consani et Moscovici. Les paramètres de croisement de ces matrices finies ne sont pas identifiés aux zéros de `zeta`, et aucune extrapolation numérique n'est considérée comme une preuve.

## État du banc

### Phase 0–2 — fondation et noyau numérique

- représentation de `s = sigma + i t` ;
- géométrie expérimentale `pi^(-sigma/2)` avec tests d'involution ;
- évaluation contrôlée de `zeta`, `xi` et de l'équation fonctionnelle ;
- matrices de Toeplitz symétriques et solveur spectral en Rust.

### Phase 3 — opérateur archimédien source-locké

La construction `Q_epsilon` / prolate de Connes–Consani a été reproduite numériquement, notamment le premier benchmark publié et les valeurs propres de la compression de Toeplitz de référence.

Documentation :

- [`docs/Q_EPSILON_SPEC.md`](docs/Q_EPSILON_SPEC.md)
- [`docs/TOEPLITZ_WEIL.md`](docs/TOEPLITZ_WEIL.md)
- [`docs/PHASE3_CHECKLIST.md`](docs/PHASE3_CHECKLIST.md)

### Phase 4 — perturbation semilocale `q = 1/p`

La phase active utilise le premier coefficient exact de la récurrence de Jacobi semilocale pour construire la dérivée à `q=0` de la matrice prolate généralisée.

Résultats actuellement établis à dimension finie dans le banc :

- formule fermée de `K'(0)` dérivée du coefficient publié ;
- lemme de signe fini : toutes les dérivées de croisement `W+` sont négatives et toutes les `W-` positives ;
- validation dense ↔ EVD tridiagonale ;
- reproduction Rust directe jusqu'à `m=4096` ;
- calcul eigenvalues-only à haute dimension avec contrôle explicite de l'annulation numérique ;
- checkpoints homogènes validés jusqu'à `m=16384` pour la variation totale de premier ordre ;
- dérivation formelle du coefficient principal `1/(2*pi^2)` de la croissance `(log m)^2` ;
- première correction sans fit `3/[32*pi^2*(log 2)^2]` dans le coefficient centré `A_m` ;
- factorisation exacte par état fondamental et noyau d'inverse fini explicite ;
- validation spectrale du modèle de Laplacien libre au soft edge, avec phase asymptotique `j+1` et espacement constant.

Documentation Phase 4 :

- [`docs/PHASE4_SEMILOCAL_FINDINGS.md`](docs/PHASE4_SEMILOCAL_FINDINGS.md)
- [`docs/PHASE4_FIRST_ORDER_DERIVATION.md`](docs/PHASE4_FIRST_ORDER_DERIVATION.md)
- [`docs/PHASE4_FIRST_ORDER_SIGN_LEMMA.md`](docs/PHASE4_FIRST_ORDER_SIGN_LEMMA.md)
- [`docs/PHASE4_RUST_VALIDATION_2026-08-14.md`](docs/PHASE4_RUST_VALIDATION_2026-08-14.md)
- [`docs/PHASE4_NUMERICAL_CHECKPOINTS_2026-08-15.md`](docs/PHASE4_NUMERICAL_CHECKPOINTS_2026-08-15.md)
- [`docs/PHASE4_LOCAL_IMPLEMENTATION_MANIFEST.md`](docs/PHASE4_LOCAL_IMPLEMENTATION_MANIFEST.md) — inventaire du code Phase 4 validé sur Thor mais pas encore upstreamé proprement.
- [`docs/PHASE4_LOG2_ASYMPTOTIC_HEURISTIC.md`](docs/PHASE4_LOG2_ASYMPTOTIC_HEURISTIC.md)
- [`docs/PHASE4_FIRST_FINITE_SIZE_CORRECTION.md`](docs/PHASE4_FIRST_FINITE_SIZE_CORRECTION.md)
- [`docs/PHASE4_EXACT_GROUND_STATE_FACTORISATION.md`](docs/PHASE4_EXACT_GROUND_STATE_FACTORISATION.md)
- [`docs/PHASE4_FREE_LAPLACIAN_SPECTRAL_VALIDATION_2026-08-15.md`](docs/PHASE4_FREE_LAPLACIAN_SPECTRAL_VALIDATION_2026-08-15.md)

## Checkpoints numériques actuels

Pour

```text
S(m) = sum_j |lambda'_j(0)|,
```

les checkpoints Rust retenus sont :

| m | S(m) |
|---:|---:|
| 128 | 3.970845543531 |
| 256 | 4.640481894221 |
| 512 | 5.359223651882 |
| 1024 | 6.126883687871 |
| 2048 | 6.943355708182 |
| 4096 | 7.808580171428 |
| 8192 | 8.72252476599664 |
| 16384 | 9.68517075774107 |

Les deux derniers points utilisent une différence finie mode-par-mode des valeurs propres, avec validation indépendante à l'intérieur de la fenêtre de pas stable. Les extrapolations contaminées par l'annulation numérique sont explicitement rejetées dans la documentation.

Le coefficient centré

```text
A_m = [S(2m) - 2S(m) + S(m/2)] / [2(log 2)^2]
```

donne à présent

```text
A_4096 = 0.05070228504273...
A_8192 = 0.05068278870666...
```

La constante

```text
1/(2*pi^2) = 0.05066059182116...
```

est désormais motivée par les coefficients locaux du soft edge et pas seulement par un fit numérique. La dérivation globale reste formelle tant que la trace singulière `Tr(K^(-1/2) H)` n'est pas contrôlée uniformément.

## Prochaine cible

La priorité est maintenant d'exploiter la factorisation exacte

```text
T_m^(-1) = R_m R_m^T
```

avec l'opérateur triangulaire de Hardy/Copson explicite `R_m`, ainsi que la validation du Laplacien libre au bas du spectre, pour obtenir un contrôle uniforme de

```text
Tr(T_m^(-1/2) H_m).
```

## Références principales

- Alain Connes, Caterina Consani, Henri Moscovici, *On q-series and the moment problem associated to local factors*, arXiv:2403.01247.
- Alain Connes, Caterina Consani, Henri Moscovici, *Zeta zeros and prolate wave operators*, arXiv:2310.18423.
- Alain Connes, Caterina Consani, *Weil positivity and Trace formula, the archimedean place*, arXiv:2006.13771.
- Grzegorz Świderski, *Periodic perturbations of unbounded Jacobi matrices III: The soft edge regime*, arXiv:1707.06486.

## Principe du projet

Un phénomène numérique n'est jamais considéré comme une preuve. Toute signature observée doit, autant que possible, survivre :

- à l'augmentation de la dimension ;
- à l'augmentation de la précision ;
- au changement de discrétisation ;
- au changement de méthode de sommation ;
- à une implémentation ou référence indépendante.
