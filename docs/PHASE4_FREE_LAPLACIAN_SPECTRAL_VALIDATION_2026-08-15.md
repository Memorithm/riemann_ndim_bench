# Phase 4 free-Laplacian spectral validation — 2026-08-15

This note records the Rust spectral check of the soft-edge continuum model inferred from the exact ground-state / Hardy factorisation.

## Setup

For each parity block and block size `m`, the test computes the positive square-root spectrum `lambda_j` of the archimedean generalized-prolate block at `q=0`, together with the exact discrete Liouville length

```text
L_m = sum_{i<m} sqrt(B_d / [(d+1/2)(d+3/2)]),
d = 2i + epsilon.
```

The leading asymptotic length is

```text
L_m ~ 4 sqrt(pi m).
```

The free one-dimensional Laplacian model predicts asymptotically constant spacing in the square-root spectrum.

## Main numerical result

The measured scaled spacings

```text
(lambda_{j+1}-lambda_j) L_m / pi
```

are nearly independent of rank and converge to 1 as `m` increases.

### W+

```text
m=1024   mean_scaled_spacing = 0.9754724148198896
m=4096   mean_scaled_spacing = 0.9878090889897554
m=16384  mean_scaled_spacing = 0.9939026771501716
```

### W-

```text
m=1024   mean_scaled_spacing = 0.9696409903039694
m=4096   mean_scaled_spacing = 0.9848443657784803
m=16384  mean_scaled_spacing = 0.9924185037937292
```

The RMS spacing error over the first eight gaps decreases accordingly:

```text
W+ : 2.4529e-2 -> 1.2191e-2 -> 6.0973e-3
W- : 3.0359e-2 -> 1.5156e-2 -> 7.5815e-3
```

This is strong finite-size evidence for the free-Laplacian spacing law.

## Boundary phase: correction of the earlier expectation

The phase is not the initially guessed Neumann–Dirichlet value `j+1/2`.

At `m=16384`:

```text
W+ : lambda_0 L_m/pi = 0.9939034484480231
W- : lambda_0 L_m/pi = 0.9924191935620239
```

and the higher ranks are close to successive integers. Thus the observed finite spectrum is compatible with

```text
lambda_j L_m/pi -> j+1
```

rather than `j+1/2`.

This corrects the previous informal boundary-phase statement. The singular left endpoint / Liouville transformation must be treated before assigning a regular Neumann condition to the transformed free operator.

## A sharper finite-size observation

The spacing error is almost entirely explained by the difference between the exact discrete accumulated length `L_m` and its leading asymptotic `4 sqrt(pi m)`.

The measured length ratios are

```text
W+ :
  m=1024   L_m/[4 sqrt(pi m)] = 0.9757663385137011
  m=4096   L_m/[4 sqrt(pi m)] = 0.9878373904467580
  m=16384  L_m/[4 sqrt(pi m)] = 0.9939072509788144

W- :
  m=1024   0.9700728450707514
  m=4096   0.9849296411558522
  m=16384  0.9924381195814174
```

Dividing the measured mean spacing by these ratios gives

```text
W+ :
  0.9996987765592946
  0.9999713500852708
  0.9999953981333386

W- :
  0.9995548223321822
  0.9999134198283729
  0.9999802347498538
```

Likewise the first scaled eigenvalue divided by the same length ratio is

```text
W+ : 0.9999387928489577, 0.9999847316210562, 0.9999961741593217
W- : 0.9996947998049540, 0.9999237036070180, 0.9999809297738368
```

Therefore the leading spectral scale is numerically much closer to

```text
lambda_j * 4 sqrt(pi m) / pi -> j+1
```

than to a finite-size law using the uncorrected discrete sum `L_m` directly.

Interpretation: `L_m` has an additive boundary correction of order one, hence a relative `O(m^-1/2)` deficit. The spectral phase/endpoint correction compensates that term to high accuracy, leaving a much smaller residual.

This observation is numerical and should not yet be promoted to a theorem.

## Scientific boundary

The result concerns the low spectrum of finite archimedean generalized-prolate blocks in the Phase 4 model. It does not identify those finite spectral points with Riemann-zeta zeros and does not imply the Riemann hypothesis.
