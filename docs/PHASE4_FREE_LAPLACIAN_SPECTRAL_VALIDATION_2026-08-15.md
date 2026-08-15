# Phase 4 free-Laplacian spectral validation — 2026-08-15

This note records the Rust spectral checks of the soft-edge continuum model inferred from the exact ground-state / Hardy factorisation.

## Setup

For each parity block and block size `m`, the tests compute the positive square-root spectrum `lambda_j` of the archimedean generalized-prolate block at `q=0`, together with the exact discrete Liouville length

```text
L_m = sum_{i<m} sqrt(B_d / [(d+1/2)(d+3/2)]),
d = 2i + epsilon.
```

The leading asymptotic length is

```text
L_m ~ 4 sqrt(pi m).
```

The free one-dimensional Laplacian model predicts asymptotically constant spacing in the square-root spectrum.

## 1. Constant-spacing test

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

## 2. Boundary phase: correction of the earlier expectation

The phase is not the initially guessed Neumann-Dirichlet value `j+1/2`.

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

This corrects the earlier informal boundary-phase statement. The singular left endpoint / Liouville transformation must be treated before assigning a regular boundary condition to the transformed free operator.

## 3. Cleaner leading length scale

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

Therefore the leading fixed-rank law is numerically much cleaner in the form

```text
lambda_j * 4 sqrt(pi m) / pi -> j+1.
```

## 4. First parity-dependent boundary correction

The dedicated test `validates_fixed_rank_free_laplacian_first_correction` examines

```text
m * [ lambda_j * 4 sqrt(pi m) / (pi (j+1)) - 1 ].
```

For the first mode the data are consistent with

```text
W+ : -> -1/16,
W- : -> -5/16.
```

This is equivalent to the effective size shift

```text
m_eff = m + (4 epsilon + 1)/8,
```

namely

```text
W+ : m_eff = m + 1/8,
W- : m_eff = m + 5/8.
```

Using

```text
lambda_j
 * 4 sqrt(pi m_eff)
 / [pi (j+1)]
```

removes the full `1/m` boundary correction for fixed rank.

At `m=16384`, rank zero gives

```text
W+ shifted error = -1.116528292044450e-8
W- shifted error = +2.714537039594234e-9.
```

For ranks `1,...,7`, the shifted residuals remain small and their decrease from `m=4096` to `m=16384` is approximately the factor expected for an `O(m^-2)` term. Representative rank-7 values are

```text
W+ : -8.209717912044390e-6 -> -4.863910976204977e-7
W- : -6.492141007119479e-6 -> -4.327108952262293e-7.
```

Thus the first correction behaves like a parity-dependent boundary displacement, while the next term is rank dependent and is compatible with the first discrete/non-uniform dispersion correction.

The rank-zero quantity multiplied by `m` is sensitive to eigensolver roundoff at the largest dimension because it amplifies relative errors by `m`; the unamplified shifted ratio remains stable at `1e-8` to `1e-9`.

## 5. Current asymptotic conjecture for fixed low rank

The finite data support

```text
lambda_j(m)
 = pi (j+1)
   / [4 sqrt(pi (m + (4 epsilon + 1)/8))]
   * [1 + O_j(m^-2)]
```

for fixed low rank `j`, separately in the two parity blocks.

This formula is a numerical asymptotic conjecture, not yet a theorem. The next target is to derive the rank-dependent `m^-2` coefficient from the exact weighted-difference / Hardy representation.

## Scientific boundary

The result concerns the low spectrum of finite archimedean generalized-prolate blocks in the Phase 4 model. It does not identify those finite spectral points with Riemann-zeta zeros and does not imply the Riemann hypothesis.
