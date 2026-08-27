//! Exact first-order `q=1/p` derivative of the finite semilocal prolate blocks.
//!
//! This module implements the source-derived coefficient recorded in
//! `docs/PHASE4_FIRST_ORDER_DERIVATION.md`. It is a finite-compression
//! perturbation calculation only; it does not identify crossings with zeta
//! zeros and carries no direct RH implication.

use std::cmp::Ordering;
use std::f64::consts::{PI, SQRT_2};
use std::fmt;

use faer::linalg::solvers::{EvdError, SelfAdjointEigen};
use faer::{Mat, Side};

/// Parity block of the generalized prolate compression.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProlateParity {
    /// `W+`, with degrees `d = 2 i`.
    WPlus,
    /// `W-`, with degrees `d = 2 i + 1`.
    WMinus,
}

impl ProlateParity {
    /// Degree associated with row `i` of this parity block.
    #[inline]
    pub const fn degree(self, i: usize) -> usize {
        match self {
            Self::WPlus => 2 * i,
            Self::WMinus => 2 * i + 1,
        }
    }

    /// Sign that turns `K'(0)` into the positive-definite matrix of the sign lemma.
    #[inline]
    pub const fn sign_correction(self) -> f64 {
        match self {
            Self::WPlus => -1.0,
            Self::WMinus => 1.0,
        }
    }
}

/// Symmetric tridiagonal matrix stored by diagonal and upper subdiagonal.
#[derive(Clone, Debug)]
pub struct SymmetricTridiagonal {
    diag: Vec<f64>,
    off_diag: Vec<f64>,
}

impl SymmetricTridiagonal {
    fn new(diag: Vec<f64>, off_diag: Vec<f64>) -> Self {
        debug_assert!(diag.is_empty() || off_diag.len() + 1 == diag.len());
        Self { diag, off_diag }
    }

    /// Matrix dimension.
    #[inline]
    pub fn len(&self) -> usize {
        self.diag.len()
    }

    /// Whether the matrix is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.diag.is_empty()
    }

    /// Diagonal entries.
    #[inline]
    pub fn diagonal(&self) -> &[f64] {
        &self.diag
    }

    /// Upper/lower off-diagonal entries.
    #[inline]
    pub fn off_diagonal(&self) -> &[f64] {
        &self.off_diag
    }

    /// Convert to a dense symmetric `faer::Mat`.
    pub fn to_dense(&self) -> Mat<f64> {
        let n = self.len();
        Mat::from_fn(n, n, |i, j| {
            if i == j {
                self.diag[i]
            } else if i + 1 == j {
                self.off_diag[i]
            } else if j + 1 == i {
                self.off_diag[j]
            } else {
                0.0
            }
        })
    }

    /// Evaluate `u^T M u` in O(n) using the tridiagonal structure.
    fn quadratic_form_column(&self, u: faer::MatRef<'_, f64>, column: usize) -> f64 {
        let mut value = 0.0;

        for i in 0..self.len() {
            let ui = u.read(i, column);
            value += self.diag[i] * ui * ui;
        }

        for i in 0..self.off_diag.len() {
            value += 2.0
                * self.off_diag[i]
                * u.read(i, column)
                * u.read(i + 1, column);
        }

        value
    }
}

/// Error returned by the finite crossing derivative calculation.
#[derive(Debug)]
pub enum SemilocalError {
    /// `faer` self-adjoint eigensolver failure.
    Eigensolver(EvdError),
    /// `K(0)` unexpectedly produced a non-positive eigenvalue.
    NonPositiveEigenvalue { index: usize, value: f64 },
}

impl fmt::Display for SemilocalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Eigensolver(error) => {
                write!(f, "self-adjoint eigensolver failed: {error:?}")
            }
            Self::NonPositiveEigenvalue { index, value } => write!(
                f,
                "K(0) eigenvalue at index {index} is not strictly positive: {value:e}"
            ),
        }
    }
}

impl std::error::Error for SemilocalError {}

impl From<EvdError> for SemilocalError {
    fn from(value: EvdError) -> Self {
        Self::Eigensolver(value)
    }
}

/// One finite crossing and its first `q=0` derivative.
#[derive(Clone, Copy, Debug)]
pub struct CrossingDerivative {
    pub parity: ProlateParity,
    pub parity_index: usize,
    pub mu: f64,
    pub lambda: f64,
    pub mu_prime: f64,
    pub lambda_prime: f64,
}

/// Aggregates of the merged normalized response `v_j = lambda'_j(0)/sqrt(m)`.
#[derive(Clone, Copy, Debug)]
pub struct ResponseStats {
    pub mean_abs: f64,
    pub trimmed_mean_abs: f64,
    pub rms: f64,
    pub linf: f64,
}

/// Generate `alpha_0, ..., alpha_{count-1}` stably by recurrence.
///
/// `alpha_0 = 1` and
/// `alpha_{n+1} = -((n+1/2)/(n+1)) alpha_n`.
pub fn alpha_sequence(count: usize) -> Vec<f64> {
    if count == 0 {
        return Vec::new();
    }

    let mut alpha = vec![0.0; count];
    alpha[0] = 1.0;

    for n in 0..count - 1 {
        let nf = n as f64;
        alpha[n + 1] = -((nf + 0.5) / (nf + 1.0)) * alpha[n];
    }

    alpha
}

/// Stable single-value form of [`alpha_sequence`].
pub fn alpha(n: usize) -> f64 {
    let mut value = 1.0;

    for k in 0..n {
        let kf = k as f64;
        value *= -((kf + 0.5) / (kf + 1.0));
    }

    value
}

/// Exact archimedean coefficient squared at `q=0`:
/// `a_n(0)^2 = (n+1/2)(n+1)`.
#[inline]
pub fn archimedean_a2(n: usize) -> f64 {
    let n = n as f64;
    (n + 0.5) * (n + 1.0)
}

/// Source-derived first derivative
/// `(a_n^2)'(0) = -(1/sqrt(2))(2n+1)(4n+3) alpha_n`.
#[inline]
pub fn archimedean_a2_prime(n: usize, alpha_n: f64) -> f64 {
    let n = n as f64;
    -((2.0 * n + 1.0) * (4.0 * n + 3.0) * alpha_n) / SQRT_2
}

#[inline]
fn b_degree(d: usize) -> f64 {
    2.0 * PI * (4 * d + 1) as f64
}

#[inline]
fn r_n(alpha_n: f64, alpha_n_plus_1: f64) -> f64 {
    SQRT_2 * (alpha_n_plus_1 - alpha_n)
}

/// Build the exact archimedean finite matrix `K(0)` for one parity block.
pub fn build_k0(block_size: usize, parity: ProlateParity) -> SymmetricTridiagonal {
    if block_size == 0 {
        return SymmetricTridiagonal::new(Vec::new(), Vec::new());
    }

    let mut diag = Vec::with_capacity(block_size);
    let mut off_diag = Vec::with_capacity(block_size.saturating_sub(1));

    for i in 0..block_size {
        let d = parity.degree(i);
        let left = if d == 0 {
            0.0
        } else {
            archimedean_a2(d - 1)
        };
        let right = archimedean_a2(d);
        diag.push((left + right + 0.25) / b_degree(d));

        if i + 1 < block_size {
            let numerator = (archimedean_a2(d) * archimedean_a2(d + 1)).sqrt();
            let denominator = (b_degree(d) * b_degree(d + 2)).sqrt();
            off_diag.push(numerator / denominator);
        }
    }

    SymmetricTridiagonal::new(diag, off_diag)
}

/// Build `K'(0)` from the closed forms of the Phase-4 derivation note.
pub fn build_kprime_closed(
    block_size: usize,
    parity: ProlateParity,
) -> SymmetricTridiagonal {
    if block_size == 0 {
        return SymmetricTridiagonal::new(Vec::new(), Vec::new());
    }

    let max_degree = parity.degree(block_size - 1);
    let alphas = alpha_sequence(max_degree + 1);
    let mut diag = Vec::with_capacity(block_size);
    let mut off_diag = Vec::with_capacity(block_size.saturating_sub(1));

    for i in 0..block_size {
        let d = parity.degree(i);
        let alpha_d = alphas[d];
        diag.push(-3.0 * alpha_d / (2.0 * SQRT_2 * PI));

        if i + 1 < block_size {
            let d_f = d as f64;
            let radical = ((2.0 * d_f + 1.0) * (2.0 * d_f + 3.0)
                / ((d_f + 1.0)
                    * (d_f + 2.0)
                    * (4.0 * d_f + 1.0)
                    * (4.0 * d_f + 9.0)))
                .sqrt();
            let value =
                -SQRT_2 * (4.0 * d_f + 5.0) * alpha_d * radical / (16.0 * PI);
            off_diag.push(value);
        }
    }

    SymmetricTridiagonal::new(diag, off_diag)
}

/// Build `K'(0)` directly from the unsimplified source-derived formulas.
///
/// This is intentionally separate from [`build_kprime_closed`] so tests can
/// detect algebraic regressions in the closed-form simplification.
pub fn build_kprime_unsimplified(
    block_size: usize,
    parity: ProlateParity,
) -> SymmetricTridiagonal {
    if block_size == 0 {
        return SymmetricTridiagonal::new(Vec::new(), Vec::new());
    }

    let max_degree = parity.degree(block_size - 1);
    let alphas = alpha_sequence(max_degree + 3);
    let mut diag = Vec::with_capacity(block_size);
    let mut off_diag = Vec::with_capacity(block_size.saturating_sub(1));

    for i in 0..block_size {
        let d = parity.degree(i);
        let left = if d == 0 {
            0.0
        } else {
            archimedean_a2_prime(d - 1, alphas[d - 1])
        };
        let right = archimedean_a2_prime(d, alphas[d]);
        diag.push((left + right) / b_degree(d));

        if i + 1 < block_size {
            let a_product = (archimedean_a2(d) * archimedean_a2(d + 1)).sqrt();
            let logarithmic_derivative = r_n(alphas[d], alphas[d + 1])
                + r_n(alphas[d + 1], alphas[d + 2]);
            off_diag.push(
                a_product * logarithmic_derivative
                    / (b_degree(d) * b_degree(d + 2)).sqrt(),
            );
        }
    }

    SymmetricTridiagonal::new(diag, off_diag)
}

/// Minimum strict diagonal-dominance margin of the sign-corrected `K'(0)`.
///
/// Returns `+infinity` for an empty block.
pub fn sign_corrected_min_diagonal_dominance_margin(
    block_size: usize,
    parity: ProlateParity,
) -> f64 {
    let kp = build_kprime_closed(block_size, parity);
    if kp.is_empty() {
        return f64::INFINITY;
    }

    let sign = parity.sign_correction();
    let mut min_margin = f64::INFINITY;

    for i in 0..kp.len() {
        let diagonal = sign * kp.diagonal()[i];
        let mut row_sum = 0.0;

        if i > 0 {
            row_sum += kp.off_diagonal()[i - 1].abs();
        }
        if i + 1 < kp.len() {
            row_sum += kp.off_diagonal()[i].abs();
        }

        min_margin = min_margin.min(diagonal - row_sum);
    }

    min_margin
}

/// Compute every finite crossing derivative for one parity block using
/// `faer::linalg::solvers::SelfAdjointEigen` and its `U()` / `S()` factors.
pub fn crossing_derivatives(
    block_size: usize,
    parity: ProlateParity,
) -> Result<Vec<CrossingDerivative>, SemilocalError> {
    if block_size == 0 {
        return Ok(Vec::new());
    }

    let k0 = build_k0(block_size, parity);
    let kprime = build_kprime_closed(block_size, parity);
    let dense = k0.to_dense();
    let evd = SelfAdjointEigen::new(dense.as_ref(), Side::Lower)?;
    let u = evd.U();
    let eigenvalues = evd.S().column_vector();

    let mut out = Vec::with_capacity(block_size);

    for j in 0..block_size {
        let mu = eigenvalues.read(j);
        if mu <= 0.0 || !mu.is_finite() {
            return Err(SemilocalError::NonPositiveEigenvalue {
                index: j,
                value: mu,
            });
        }

        let lambda = mu.sqrt();
        let mu_prime = kprime.quadratic_form_column(u, j);
        let lambda_prime = mu_prime / (2.0 * lambda);

        out.push(CrossingDerivative {
            parity,
            parity_index: j,
            mu,
            lambda,
            mu_prime,
            lambda_prime,
        });
    }

    Ok(out)
}

/// Merge the two parity spectra and sort them by increasing crossing `lambda`.
pub fn merged_crossing_derivatives(
    block_size: usize,
) -> Result<Vec<CrossingDerivative>, SemilocalError> {
    let mut merged = crossing_derivatives(block_size, ProlateParity::WPlus)?;
    merged.extend(crossing_derivatives(block_size, ProlateParity::WMinus)?);
    merged.sort_by(|left, right| {
        left.lambda
            .partial_cmp(&right.lambda)
            .unwrap_or(Ordering::Equal)
    });
    Ok(merged)
}

/// Statistics of the merged normalized response.
///
/// One eighth of the merged spectrum is removed from **each** edge for the
/// trimmed mean, matching the Phase-4 regression tables.
pub fn merged_response_stats(block_size: usize) -> Result<ResponseStats, SemilocalError> {
    if block_size == 0 {
        return Ok(ResponseStats {
            mean_abs: 0.0,
            trimmed_mean_abs: 0.0,
            rms: 0.0,
            linf: 0.0,
        });
    }

    let merged = merged_crossing_derivatives(block_size)?;
    let normalization = (block_size as f64).sqrt();
    let values: Vec<f64> = merged
        .iter()
        .map(|crossing| crossing.lambda_prime / normalization)
        .collect();

    let count = values.len() as f64;
    let mean_abs = values.iter().map(|value| value.abs()).sum::<f64>() / count;
    let rms = (values.iter().map(|value| value * value).sum::<f64>() / count).sqrt();
    let linf = values
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);

    let trim = values.len() / 8;
    let trimmed = &values[trim..values.len() - trim];
    let trimmed_mean_abs = if trimmed.is_empty() {
        0.0
    } else {
        trimmed.iter().map(|value| value.abs()).sum::<f64>() / trimmed.len() as f64
    };

    Ok(ResponseStats {
        mean_abs,
        trimmed_mean_abs,
        rms,
        linf,
    })
}

/// Total unnormalised first-order response `sum_j |lambda'_j(0)|`.
pub fn merged_total_abs_derivative(block_size: usize) -> Result<f64, SemilocalError> {
    Ok(merged_crossing_derivatives(block_size)?
        .iter()
        .map(|crossing| crossing.lambda_prime.abs())
        .sum())
}
