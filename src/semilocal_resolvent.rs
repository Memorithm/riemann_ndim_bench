//! O(m) shifted resolvent tools for the finite semilocal prolate blocks.
//!
//! The Phase-4 soft-edge proof route rewrites the singular square-root trace
//! through the exact finite-dimensional resolvent identity. This module
//! supplies finite tridiagonal kernels without forming a dense inverse.

use std::fmt;

use crate::semilocal::{ProlateParity, build_k0, build_kprime_closed};

/// Error returned by the shifted tridiagonal resolvent calculations.
#[derive(Debug)]
pub enum ResolventTraceError {
    /// The shift in `K + t I` must be finite and non-negative.
    InvalidShift { value: f64 },
    /// A row-local frozen model requires both neighboring edges.
    InvalidInteriorRow { row: usize, len: usize },
    /// A non-positive or non-finite Schur/LDL denominator would violate
    /// positive definiteness of the shifted finite matrix.
    NonPositivePivot { index: usize, value: f64 },
    /// The frozen Toeplitz denominator must stay strictly positive.
    NonPositiveFrozenDiscriminant { row: usize, value: f64 },
}

impl fmt::Display for ResolventTraceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidShift { value } => {
                write!(
                    f,
                    "resolvent shift must be finite and non-negative: {value:e}"
                )
            }
            Self::InvalidInteriorRow { row, len } => write!(
                f,
                "frozen resolvent model requires an interior row: row={row}, len={len}"
            ),
            Self::NonPositivePivot { index, value } => write!(
                f,
                "shifted tridiagonal Schur/LDL denominator at index {index} is not positive: {value:e}"
            ),
            Self::NonPositiveFrozenDiscriminant { row, value } => write!(
                f,
                "frozen resolvent symbol at row {row} is not strictly positive: {value:e}"
            ),
        }
    }
}

impl std::error::Error for ResolventTraceError {}

fn checked_shift(shift: f64) -> Result<(), ResolventTraceError> {
    if shift.is_finite() && shift >= 0.0 {
        Ok(())
    } else {
        Err(ResolventTraceError::InvalidShift { value: shift })
    }
}

fn checked_pivot(index: usize, value: f64) -> Result<f64, ResolventTraceError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(ResolventTraceError::NonPositivePivot { index, value })
    }
}

/// Exact finite Green-function data reconstructed from left/right Schur
/// complements of a shifted positive tridiagonal matrix.
#[derive(Clone, Debug)]
pub struct CavityGreenBands {
    left_denominators: Vec<f64>,
    right_denominators: Vec<f64>,
    diagonal: Vec<f64>,
    off_diagonal: Vec<f64>,
}

impl CavityGreenBands {
    /// Left-to-right cavity denominators.
    #[inline]
    pub fn left_denominators(&self) -> &[f64] {
        &self.left_denominators
    }

    /// Right-to-left cavity denominators.
    #[inline]
    pub fn right_denominators(&self) -> &[f64] {
        &self.right_denominators
    }

    /// Diagonal of `(K+tI)^(-1)`.
    #[inline]
    pub fn diagonal(&self) -> &[f64] {
        &self.diagonal
    }

    /// First upper/lower off-diagonal of `(K+tI)^(-1)`.
    #[inline]
    pub fn off_diagonal(&self) -> &[f64] {
        &self.off_diagonal
    }
}

/// Reusable coefficients for
/// `Tr[(K(0) + shift I)^(-1) H]` on one finite parity block.
///
/// Here `H=-K'_+(0)` on `W+` and `H=K'_-(0)` on `W-`, so `H` is the
/// positive sign-corrected perturbation from the finite-block sign lemma.
/// Constructing the kernel is O(m); each subsequent shifted trace is O(m) and
/// uses O(m) auxiliary storage without rebuilding `K(0)` or `H`.
#[derive(Clone, Debug)]
pub struct SignCorrectedResolventTraceKernel {
    diagonal: Vec<f64>,
    off_diagonal: Vec<f64>,
    h_diagonal: Vec<f64>,
    h_off_diagonal: Vec<f64>,
}

impl SignCorrectedResolventTraceKernel {
    /// Build the reusable tridiagonal coefficients for one parity block.
    pub fn new(block_size: usize, parity: ProlateParity) -> Self {
        let k0 = build_k0(block_size, parity);
        let kprime = build_kprime_closed(block_size, parity);
        let sign = parity.sign_correction();

        Self {
            diagonal: k0.diagonal().to_vec(),
            off_diagonal: k0.off_diagonal().to_vec(),
            h_diagonal: kprime.diagonal().iter().map(|value| sign * value).collect(),
            h_off_diagonal: kprime
                .off_diagonal()
                .iter()
                .map(|value| sign * value)
                .collect(),
        }
    }

    /// Matrix dimension represented by this kernel.
    #[inline]
    pub fn len(&self) -> usize {
        self.diagonal.len()
    }

    /// Whether this kernel represents an empty block.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.diagonal.is_empty()
    }

    fn selected_inverse_bands(
        &self,
        shift: f64,
    ) -> Result<(Vec<f64>, Vec<f64>), ResolventTraceError> {
        checked_shift(shift)?;
        if self.is_empty() {
            return Ok((Vec::new(), Vec::new()));
        }

        let block_size = self.len();
        let mut pivots = Vec::with_capacity(block_size);
        let mut multipliers = Vec::with_capacity(block_size.saturating_sub(1));

        pivots.push(checked_pivot(0, self.diagonal[0] + shift)?);
        for i in 1..block_size {
            let off = self.off_diagonal[i - 1];
            let multiplier = off / pivots[i - 1];
            let pivot = self.diagonal[i] + shift - multiplier * off;
            multipliers.push(multiplier);
            pivots.push(checked_pivot(i, pivot)?);
        }

        let mut inverse_diagonal = vec![0.0; block_size];
        let mut inverse_off_diagonal = vec![0.0; block_size.saturating_sub(1)];
        let last = block_size - 1;
        inverse_diagonal[last] = 1.0 / pivots[last];

        for i in (0..last).rev() {
            let multiplier = multipliers[i];
            let next_inverse_diagonal = inverse_diagonal[i + 1];
            inverse_off_diagonal[i] = -multiplier * next_inverse_diagonal;
            inverse_diagonal[i] = 1.0 / pivots[i] + multiplier * multiplier * next_inverse_diagonal;
        }

        Ok((inverse_diagonal, inverse_off_diagonal))
    }

    /// Reconstruct the exact finite Green diagonal and first off-diagonal from
    /// independent left/right Schur-complement continued fractions.
    ///
    /// For shifted diagonal `a_i` and edge `b_i`, define
    ///
    /// `L_0=a_0`, `L_i=a_i-b_{i-1}^2/L_{i-1}`
    ///
    /// and
    ///
    /// `R_{m-1}=a_{m-1}`, `R_i=a_i-b_i^2/R_{i+1}`.
    ///
    /// Then the exact finite Green diagonal is
    ///
    /// `G_ii = 1/(a_i-b_{i-1}^2/L_{i-1}-b_i^2/R_{i+1})`,
    ///
    /// with absent boundary terms omitted, and
    ///
    /// `G_i,i+1 = -b_i G_ii/R_{i+1}`.
    ///
    /// The two cavity recurrences explicitly retain both finite-section
    /// boundaries. No frozen/local approximation is used.
    pub fn cavity_green_bands(&self, shift: f64) -> Result<CavityGreenBands, ResolventTraceError> {
        checked_shift(shift)?;
        if self.is_empty() {
            return Ok(CavityGreenBands {
                left_denominators: Vec::new(),
                right_denominators: Vec::new(),
                diagonal: Vec::new(),
                off_diagonal: Vec::new(),
            });
        }

        let n = self.len();
        let mut left = vec![0.0; n];
        let mut right = vec![0.0; n];

        left[0] = checked_pivot(0, self.diagonal[0] + shift)?;
        for i in 1..n {
            let edge = self.off_diagonal[i - 1];
            let denominator = self.diagonal[i] + shift - edge * edge / left[i - 1];
            left[i] = checked_pivot(i, denominator)?;
        }

        let last = n - 1;
        right[last] = checked_pivot(last, self.diagonal[last] + shift)?;
        for i in (0..last).rev() {
            let edge = self.off_diagonal[i];
            let denominator = self.diagonal[i] + shift - edge * edge / right[i + 1];
            right[i] = checked_pivot(i, denominator)?;
        }

        let mut diagonal = vec![0.0; n];
        for i in 0..n {
            let mut denominator = self.diagonal[i] + shift;
            if i > 0 {
                let edge = self.off_diagonal[i - 1];
                denominator -= edge * edge / left[i - 1];
            }
            if i + 1 < n {
                let edge = self.off_diagonal[i];
                denominator -= edge * edge / right[i + 1];
            }
            diagonal[i] = 1.0 / checked_pivot(i, denominator)?;
        }

        let mut off_diagonal = vec![0.0; n.saturating_sub(1)];
        for i in 0..off_diagonal.len() {
            off_diagonal[i] = -self.off_diagonal[i] * diagonal[i] / right[i + 1];
        }

        Ok(CavityGreenBands {
            left_denominators: left,
            right_denominators: right,
            diagonal,
            off_diagonal,
        })
    }

    /// Exact row contributions to `Tr[(K(0) + shift I)^(-1) H]`.
    ///
    /// Row `i` is `(R H)_{ii}` for `R=(K(0)+shift I)^(-1)`:
    /// `R_ii H_ii + R_i,i-1 H_i-1,i + R_i,i+1 H_i+1,i`.
    /// Summing the returned vector therefore recovers the full weighted trace
    /// exactly up to floating-point summation order. This row resolution is a
    /// finite identity; no local-symbol or large-m approximation is used.
    pub fn row_contributions(&self, shift: f64) -> Result<Vec<f64>, ResolventTraceError> {
        let (inverse_diagonal, inverse_off_diagonal) = self.selected_inverse_bands(shift)?;
        let mut rows = Vec::with_capacity(self.len());

        for i in 0..self.len() {
            let mut contribution = self.h_diagonal[i] * inverse_diagonal[i];
            if i > 0 {
                contribution += self.h_off_diagonal[i - 1] * inverse_off_diagonal[i - 1];
            }
            if i + 1 < self.len() {
                contribution += self.h_off_diagonal[i] * inverse_off_diagonal[i];
            }
            rows.push(contribution);
        }

        Ok(rows)
    }

    /// Closed frozen-row Toeplitz model for the weighted resolvent density.
    ///
    /// At an interior row, freeze the two neighboring `K` edges to their
    /// arithmetic mean `b` and the two neighboring `H` edges to `o`. After
    /// alternating conjugation the scalar symbols are
    ///
    /// `k(theta) = a - 2 b cos(theta)` and
    /// `h(theta) = d - 2 o cos(theta)`,
    ///
    /// where `a=K_ii+shift` and `d=H_ii`. The returned value is the exact
    /// infinite-Toeplitz row integral
    ///
    /// `(1/(2 pi)) integral h(theta)/k(theta) dtheta`.
    ///
    /// This is a local comparison model only. The method makes no claim that
    /// the finite-section row contribution converges to it at a particular
    /// rate.
    pub fn frozen_row_resolvent_density(
        &self,
        row: usize,
        shift: f64,
    ) -> Result<f64, ResolventTraceError> {
        checked_shift(shift)?;
        if row == 0 || row >= self.len().saturating_sub(1) {
            return Err(ResolventTraceError::InvalidInteriorRow {
                row,
                len: self.len(),
            });
        }

        let b = 0.5 * (self.off_diagonal[row - 1] + self.off_diagonal[row]);
        let o = 0.5 * (self.h_off_diagonal[row - 1] + self.h_off_diagonal[row]);
        let a = self.diagonal[row] + shift;
        let c = 2.0 * b;
        let discriminant = (a - c) * (a + c);
        if !discriminant.is_finite() || discriminant <= 0.0 {
            return Err(ResolventTraceError::NonPositiveFrozenDiscriminant {
                row,
                value: discriminant,
            });
        }

        let i0 = 1.0 / discriminant.sqrt();
        let i1 = (a * i0 - 1.0) / c;
        Ok(self.h_diagonal[row] * i0 - 2.0 * o * i1)
    }

    /// Compute `Tr[(K(0) + shift I)^(-1) H]` in O(m) time and storage.
    pub fn trace(&self, shift: f64) -> Result<f64, ResolventTraceError> {
        Ok(self.row_contributions(shift)?.into_iter().sum())
    }
}

/// Convenience wrapper for a single shifted trace evaluation.
pub fn sign_corrected_resolvent_trace_tridiagonal(
    block_size: usize,
    parity: ProlateParity,
    shift: f64,
) -> Result<f64, ResolventTraceError> {
    SignCorrectedResolventTraceKernel::new(block_size, parity).trace(shift)
}
