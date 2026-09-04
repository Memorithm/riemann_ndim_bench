# riemann_ndim_bench

Banc de recherche numérique en Rust autour de l'hypothèse de Riemann, conçu pour séparer strictement :

1. les identités mathématiques exactes ;
2. les constructions numériques vérifiables ;
3. les hypothèses et asymptotiques expérimentales.

## Statut scientifique

Ce dépôt **ne contient pas de preuve de l'hypothèse de Riemann**.

Les calculs actuels portent sur des compressions finies d'opérateurs de type prolate/semilocal inspirés des travaux de Connes, Consani et Moscovici. Les paramètres de croisement de ces matrices finies ne sont pas identifiés aux zéros de `zeta`, et aucune extrapolation numérique n'est considérée comme une preuve.

Le graphe de dépendances vers RH et les frontières de portée du banc sont audités dans [`docs/RH_PROOF_DEPENDENCY_AUDIT.md`](docs/RH_PROOF_DEPENDENCY_AUDIT.md).

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
- validation spectrale du modèle de Laplacien libre au soft edge ;
- loi fixe-rang numérique `lambda_j * 4 sqrt(pi m)/pi -> j+1` ;
- première correction de bord compatible avec `m -> m + (4 epsilon + 1)/8`, laissant un résidu fixe-rang compatible avec `O(m^-2)` ;
- formules exactes de cavité de Schur à shift nul et de sa première réponse au shift.

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
- [`docs/PHASE4_SOFT_EDGE_PROOF_ROUTE.md`](docs/PHASE4_SOFT_EDGE_PROOF_ROUTE.md)
- [`docs/PHASE4_SOFT_EDGE_RUST_VALIDATION_2026-08-15.md`](docs/PHASE4_SOFT_EDGE_RUST_VALIDATION_2026-08-15.md)
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

est motivée par les coefficients locaux du soft edge et pas seulement par un fit numérique. La dérivation globale reste formelle. Le contrôle générique des résolvants/Jacobi lentement variables n'est toutefois plus le verrou directeur de ce dépôt : cette théorie est extraite vers la ligne opératorielle TDI-10.xx lorsqu'elle peut être formulée indépendamment de Riemann.

## Prochaine cible Riemann

La priorité est maintenant de **source-locker le pont semilocal vers la positivité de Weil**.

Le banc doit déterminer précisément :

1. quel espace de fonctions test et quelle fenêtre de support sont concernés ;
2. quelle inégalité ou positivité opératorielle semilocale serait suffisante ;
3. quel rôle exact jouent le Sonin space, les projections de cutoff, `Q epsilon`, l'opérateur prolate semilocal et l'augmentation de l'ensemble fini de places `S` ;
4. quelles flèches sont des théorèmes de la littérature et lesquelles restent conjecturales ;
5. quel premier énoncé Riemann-spécifique peut être attaqué ou falsifié par les objets déjà implémentés.

Aucune nouvelle asymptotique de matrice ne doit devenir prioritaire tant que son implication vers ce pont n'est pas documentée.

## Politique d'extraction

Lorsqu'un algorithme ou un lemme découvert ici peut être formulé et testé sans mentionner Riemann, zeta, Weil, Sonin ou la construction prolate semilocale spécifique, il devient candidat à transfert vers le banc approprié.

En particulier :

- théorie générique Jacobi/résolvant/soft-edge -> TDI-10.xx ;
- noyaux numériques scientifiques génériques mûrs -> SciRust ;
- recherche automatique d'algorithmes -> Forge ;
- kernels CUDA/NVIDIA -> NNIS.

RiemannBench conserve l'adaptateur spécifique, les hypothèses à vérifier et la conséquence vers la chaîne RH.

## Références principales

- Alain Connes, Caterina Consani, *The Scaling Hamiltonian*, arXiv:1910.14368.
- Alain Connes, Caterina Consani, *Weil positivity and Trace formula, the archimedean place*, arXiv:2006.13771.
- Alain Connes, Caterina Consani, Henri Moscovici, *Zeta zeros and prolate wave operators*, arXiv:2310.18423.
- Alain Connes, Caterina Consani, Henri Moscovici, *On q-series and the moment problem associated to local factors*, arXiv:2403.01247.
- Grzegorz Świderski, *Periodic perturbations of unbounded Jacobi matrices III: The soft edge regime*, arXiv:1707.06486.

## Principe du projet

Un phénomène numérique n'est jamais considéré comme une preuve. Toute signature observée doit, autant que possible, survivre :

- à l'augmentation de la dimension ;
- à l'augmentation de la précision ;
- au changement de discrétisation ;
- au changement de méthode de sommation ;
- à une implémentation ou référence indépendante.
