//! O(m) shifted resolvent trace for the finite semilocal prolate blocks.
//!
//! The Phase-4 soft-edge proof route rewrites the singular square-root trace
//! through the exact finite-dimensional resolvent identity. This module
//! supplies the corresponding finite tridiagonal kernel without forming an
//! inverse matrix or computing eigenvectors.

use std::fmt;

use crate::semilocal::{ProlateParity, build_k0, build_kprime_closed};

/// Error returned by the shifted tridiagonal resolvent trace.
#[derive(Debug)]
pub enum ResolventTraceError {
    /// The shift in `K + t I` must be finite and non-negative.
    InvalidShift { value: f64 },
    /// A non-positive or non-finite LDL pivot would violate positive
    /// definiteness of the shifted finite matrix.
    NonPositivePivot { index: usize, value: f64 },
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
            Self::NonPositivePivot { index, value } => write!(
                f,
                "shifted tridiagonal LDL pivot at index {index} is not positive: {value:e}"
            ),
        }
    }
}

impl std::error::Error for ResolventTraceError {}

fn checked_pivot(index: usize, value: f64) -> Result<f64, ResolventTraceError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(ResolventTraceError::NonPositivePivot { index, value })
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
        if !shift.is_finite() || shift < 0.0 {
            return Err(ResolventTraceError::InvalidShift { value: shift });
        }
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
