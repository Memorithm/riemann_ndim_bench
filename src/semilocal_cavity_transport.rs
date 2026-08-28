//! Exact transport identity for finite-cavity errors relative to the local
//! frozen Toeplitz fixed points.
//!
//! This module does not prove that the transport factors are uniformly
//! contractive or that the coefficient-drift forcing is summable.  It exposes
//! the exact finite recurrence in the form needed for those estimates.

use crate::semilocal::{ProlateParity, build_k0};
use crate::semilocal_frozen_cavity::frozen_row_cavity_fixed_point;
use crate::semilocal_resolvent::{ResolventTraceError, SignCorrectedResolventTraceKernel};

/// Exact one-row decomposition of the left/right finite-cavity errors.
#[derive(Clone, Copy, Debug)]
pub struct CavityErrorTransport {
    left_error: f64,
    left_previous_error: f64,
    left_transport_factor: f64,
    left_drift: f64,
    right_error: f64,
    right_next_error: f64,
    right_transport_factor: f64,
    right_drift: f64,
}

impl CavityErrorTransport {
    /// Current left-cavity error `L_i-q_i`.
    #[inline]
    pub fn left_error(self) -> f64 {
        self.left_error
    }

    /// Previous left-cavity error `L_{i-1}-q_{i-1}`.
    #[inline]
    pub fn left_previous_error(self) -> f64 {
        self.left_previous_error
    }

    /// Exact multiplier of the previous left error.
    #[inline]
    pub fn left_transport_factor(self) -> f64 {
        self.left_transport_factor
    }

    /// Left coefficient-drift forcing.
    #[inline]
    pub fn left_drift(self) -> f64 {
        self.left_drift
    }

    /// Current right-cavity error `R_i-q_i`.
    #[inline]
    pub fn right_error(self) -> f64 {
        self.right_error
    }

    /// Next right-cavity error `R_{i+1}-q_{i+1}`.
    #[inline]
    pub fn right_next_error(self) -> f64 {
        self.right_next_error
    }

    /// Exact multiplier of the next right error.
    #[inline]
    pub fn right_transport_factor(self) -> f64 {
        self.right_transport_factor
    }

    /// Right coefficient-drift forcing.
    #[inline]
    pub fn right_drift(self) -> f64 {
        self.right_drift
    }

    /// Reconstructed left error from transport plus drift.
    #[inline]
    pub fn reconstructed_left_error(self) -> f64 {
        self.left_transport_factor * self.left_previous_error + self.left_drift
    }

    /// Reconstructed right error from transport plus drift.
    #[inline]
    pub fn reconstructed_right_error(self) -> f64 {
        self.right_transport_factor * self.right_next_error + self.right_drift
    }
}

/// Decompose the exact finite cavity errors at one row into propagated error
/// and local coefficient drift.
///
/// Let
///
/// `L_i = a_i - e_{i-1}^2/L_{i-1}`
///
/// and
///
/// `R_i = a_i - e_i^2/R_{i+1}`
///
/// be the exact finite left/right Schur denominators, and let `q_i` denote the
/// positive fixed point of the row-frozen symmetric Toeplitz model.  Then
/// exactly
///
/// `L_i-q_i = alpha_i^- (L_{i-1}-q_{i-1}) + delta_i^-`
///
/// and
///
/// `R_i-q_i = alpha_i^+ (R_{i+1}-q_{i+1}) + delta_i^+`,
///
/// where
///
/// `alpha_i^- = e_{i-1}^2/(L_{i-1} q_{i-1})`,
///
/// `delta_i^- = a_i - e_{i-1}^2/q_{i-1} - q_i`,
///
/// `alpha_i^+ = e_i^2/(R_{i+1} q_{i+1})`,
///
/// and
///
/// `delta_i^+ = a_i - e_i^2/q_{i+1} - q_i`.
///
/// The row must have two neighbors on each side.  In addition, the frozen
/// Toeplitz symbols at `i-1`, `i`, and `i+1` must be strictly positive so
/// that all three positive fixed points exist.  Failure of that local
/// positivity condition is reported as `NonPositiveFrozenDiscriminant`
/// rather than extending the frozen model into the first exceptional rows.
pub fn cavity_error_transport(
    block_size: usize,
    row: usize,
    parity: ProlateParity,
    shift: f64,
) -> Result<CavityErrorTransport, ResolventTraceError> {
    if row < 2 || row.saturating_add(2) >= block_size {
        return Err(ResolventTraceError::InvalidInteriorRow {
            row,
            len: block_size,
        });
    }

    let k0 = build_k0(block_size, parity);
    let kernel = SignCorrectedResolventTraceKernel::new(block_size, parity);
    let cavity = kernel.cavity_green_bands(shift)?;

    let q_previous =
        frozen_row_cavity_fixed_point(block_size, row - 1, parity, shift)?.cavity_denominator();
    let q_current =
        frozen_row_cavity_fixed_point(block_size, row, parity, shift)?.cavity_denominator();
    let q_next =
        frozen_row_cavity_fixed_point(block_size, row + 1, parity, shift)?.cavity_denominator();

    let left_previous = cavity.left_denominators()[row - 1];
    let left_current = cavity.left_denominators()[row];
    let right_current = cavity.right_denominators()[row];
    let right_next = cavity.right_denominators()[row + 1];
    let left_edge = k0.off_diagonal()[row - 1];
    let right_edge = k0.off_diagonal()[row];
    let shifted_diagonal = k0.diagonal()[row] + shift;

    let left_previous_error = left_previous - q_previous;
    let left_transport_factor = left_edge * left_edge / (left_previous * q_previous);
    let left_drift = shifted_diagonal - left_edge * left_edge / q_previous - q_current;
    let left_error = left_current - q_current;

    let right_next_error = right_next - q_next;
    let right_transport_factor = right_edge * right_edge / (right_next * q_next);
    let right_drift = shifted_diagonal - right_edge * right_edge / q_next - q_current;
    let right_error = right_current - q_current;

    Ok(CavityErrorTransport {
        left_error,
        left_previous_error,
        left_transport_factor,
        left_drift,
        right_error,
        right_next_error,
        right_transport_factor,
        right_drift,
    })
}
