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

    /// Compute `Tr[(K(0) + shift I)^(-1) H]`.
    ///
    /// The implementation uses an LDL factorization followed by selected
    /// inversion recurrences for only the diagonal and first off-diagonal
    /// entries of the inverse. Since `H` is tridiagonal, those entries are
    /// sufficient for the trace.
    pub fn trace(&self, shift: f64) -> Result<f64, ResolventTraceError> {
        if !shift.is_finite() || shift < 0.0 {
            return Err(ResolventTraceError::InvalidShift { value: shift });
        }
        if self.is_empty() {
            return Ok(0.0);
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
        let last = block_size - 1;
        inverse_diagonal[last] = 1.0 / pivots[last];

        let mut trace = self.h_diagonal[last] * inverse_diagonal[last];

        for i in (0..last).rev() {
            let multiplier = multipliers[i];
            let next_inverse_diagonal = inverse_diagonal[i + 1];
            let inverse_off_diagonal = -multiplier * next_inverse_diagonal;

            inverse_diagonal[i] = 1.0 / pivots[i] + multiplier * multiplier * next_inverse_diagonal;

            trace += self.h_diagonal[i] * inverse_diagonal[i]
                + 2.0 * self.h_off_diagonal[i] * inverse_off_diagonal;
        }

        Ok(trace)
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
