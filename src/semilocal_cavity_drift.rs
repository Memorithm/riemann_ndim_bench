//! Exact factorization of the finite-to-frozen cavity transport terms.
//!
//! The transport identity separates each cavity error into a propagated
//! boundary contribution and a coefficient-drift forcing.  This module
//! decomposes those two pieces further, without introducing an asymptotic
//! approximation.

use crate::semilocal::{ProlateParity, build_k0};
use crate::semilocal_cavity_transport::{CavityErrorTransport, cavity_error_transport};
use crate::semilocal_frozen_cavity::frozen_row_cavity_fixed_point;
use crate::semilocal_resolvent::{ResolventTraceError, SignCorrectedResolventTraceKernel};

/// Exact algebraic components behind one row of the cavity transport identity.
#[derive(Clone, Copy, Debug)]
pub struct CavityDriftFactorization {
    left_edge_drift: f64,
    left_fixed_point_drift: f64,
    left_local_contraction: f64,
    left_edge_ratio: f64,
    left_cavity_ratio: f64,
    right_edge_drift: f64,
    right_fixed_point_drift: f64,
    right_local_contraction: f64,
    right_edge_ratio: f64,
    right_cavity_ratio: f64,
}

impl CavityDriftFactorization {
    /// Part of `delta_i^-` caused by replacing the incoming edge by the
    /// current row's frozen averaged edge.
    #[inline]
    pub fn left_edge_drift(self) -> f64 {
        self.left_edge_drift
    }

    /// Part of `delta_i^-` caused by the change from `q_{i-1}` to `q_i`.
    #[inline]
    pub fn left_fixed_point_drift(self) -> f64 {
        self.left_fixed_point_drift
    }

    /// Frozen contraction factor `kappa_{i-1}`.
    #[inline]
    pub fn left_local_contraction(self) -> f64 {
        self.left_local_contraction
    }

    /// Incoming-edge correction `e_{i-1}^2 / b_{i-1}^2` relative to the
    /// frozen averaged edge at row `i-1`.
    #[inline]
    pub fn left_edge_ratio(self) -> f64 {
        self.left_edge_ratio
    }

    /// Finite-cavity correction `q_{i-1}/L_{i-1}`.
    #[inline]
    pub fn left_cavity_ratio(self) -> f64 {
        self.left_cavity_ratio
    }

    /// Part of `delta_i^+` caused by replacing the outgoing edge by the
    /// current row's frozen averaged edge.
    #[inline]
    pub fn right_edge_drift(self) -> f64 {
        self.right_edge_drift
    }

    /// Part of `delta_i^+` caused by the change from `q_i` to `q_{i+1}`.
    #[inline]
    pub fn right_fixed_point_drift(self) -> f64 {
        self.right_fixed_point_drift
    }

    /// Frozen contraction factor `kappa_{i+1}`.
    #[inline]
    pub fn right_local_contraction(self) -> f64 {
        self.right_local_contraction
    }

    /// Outgoing-edge correction `e_i^2 / b_{i+1}^2` relative to the frozen
    /// averaged edge at row `i+1`.
    #[inline]
    pub fn right_edge_ratio(self) -> f64 {
        self.right_edge_ratio
    }

    /// Finite-cavity correction `q_{i+1}/R_{i+1}`.
    #[inline]
    pub fn right_cavity_ratio(self) -> f64 {
        self.right_cavity_ratio
    }

    /// Reconstruct `delta_i^-` from edge asymmetry and fixed-point variation.
    #[inline]
    pub fn reconstructed_left_drift(self) -> f64 {
        self.left_edge_drift + self.left_fixed_point_drift
    }

    /// Reconstruct `delta_i^+` from edge asymmetry and fixed-point variation.
    #[inline]
    pub fn reconstructed_right_drift(self) -> f64 {
        self.right_edge_drift + self.right_fixed_point_drift
    }

    /// Reconstruct `alpha_i^-` from frozen contraction, edge correction, and
    /// finite-cavity correction.
    #[inline]
    pub fn reconstructed_left_transport_factor(self) -> f64 {
        self.left_local_contraction * self.left_edge_ratio * self.left_cavity_ratio
    }

    /// Reconstruct `alpha_i^+` from frozen contraction, edge correction, and
    /// finite-cavity correction.
    #[inline]
    pub fn reconstructed_right_transport_factor(self) -> f64 {
        self.right_local_contraction * self.right_edge_ratio * self.right_cavity_ratio
    }
}

/// Factor the exact drift and transport multipliers at one admissible row.
///
/// Write `b_i` for the arithmetic mean of the two exact edges adjacent to row
/// `i`, and let `q_i` be the positive fixed point of the corresponding frozen
/// Toeplitz cavity map.  Since
///
/// `q_i = a_i - b_i^2/q_i`,
///
/// the left drift from the transport identity satisfies exactly
///
/// `delta_i^- = b_i^2/q_i - e_{i-1}^2/q_{i-1}`
///
/// and hence
///
/// `delta_i^-`
/// `= (b_i^2-e_{i-1}^2)/q_i`
/// `  + e_{i-1}^2 (1/q_i - 1/q_{i-1})`.
///
/// The right drift has the analogous decomposition
///
/// `delta_i^+`
/// `= (b_i^2-e_i^2)/q_i`
/// `  + e_i^2 (1/q_i - 1/q_{i+1})`.
///
/// The exact transport factors also split as
///
/// `alpha_i^-`
/// `= kappa_{i-1}`
/// `  * (e_{i-1}^2/b_{i-1}^2)`
/// `  * (q_{i-1}/L_{i-1})`,
///
/// and
///
/// `alpha_i^+`
/// `= kappa_{i+1}`
/// `  * (e_i^2/b_{i+1}^2)`
/// `  * (q_{i+1}/R_{i+1})`,
///
/// where `kappa_j=b_j^2/q_j^2` is the frozen local contraction factor.
/// These identities isolate the three quantities that a subsequent estimate
/// must control; no rate or uniform contraction bound is asserted here.
pub fn cavity_drift_factorization(
    block_size: usize,
    row: usize,
    parity: ProlateParity,
    shift: f64,
) -> Result<(CavityErrorTransport, CavityDriftFactorization), ResolventTraceError> {
    let transport = cavity_error_transport(block_size, row, parity, shift)?;
    let k0 = build_k0(block_size, parity);
    let kernel = SignCorrectedResolventTraceKernel::new(block_size, parity);
    let cavity = kernel.cavity_green_bands(shift)?;

    let previous = frozen_row_cavity_fixed_point(block_size, row - 1, parity, shift)?;
    let current = frozen_row_cavity_fixed_point(block_size, row, parity, shift)?;
    let next = frozen_row_cavity_fixed_point(block_size, row + 1, parity, shift)?;

    let q_previous = previous.cavity_denominator();
    let q_current = current.cavity_denominator();
    let q_next = next.cavity_denominator();
    let frozen_edge_current = current.edge_coefficient();
    let frozen_edge_previous = previous.edge_coefficient();
    let frozen_edge_next = next.edge_coefficient();
    let incoming_edge = k0.off_diagonal()[row - 1];
    let outgoing_edge = k0.off_diagonal()[row];
    let left_cavity = cavity.left_denominators()[row - 1];
    let right_cavity = cavity.right_denominators()[row + 1];

    let left_edge_drift =
        (frozen_edge_current * frozen_edge_current - incoming_edge * incoming_edge) / q_current;
    let left_fixed_point_drift =
        incoming_edge * incoming_edge * (1.0 / q_current - 1.0 / q_previous);
    let left_edge_ratio =
        incoming_edge * incoming_edge / (frozen_edge_previous * frozen_edge_previous);
    let left_cavity_ratio = q_previous / left_cavity;

    let right_edge_drift =
        (frozen_edge_current * frozen_edge_current - outgoing_edge * outgoing_edge) / q_current;
    let right_fixed_point_drift =
        outgoing_edge * outgoing_edge * (1.0 / q_current - 1.0 / q_next);
    let right_edge_ratio =
        outgoing_edge * outgoing_edge / (frozen_edge_next * frozen_edge_next);
    let right_cavity_ratio = q_next / right_cavity;

    let factorization = CavityDriftFactorization {
        left_edge_drift,
        left_fixed_point_drift,
        left_local_contraction: previous.contraction_factor(),
        left_edge_ratio,
        left_cavity_ratio,
        right_edge_drift,
        right_fixed_point_drift,
        right_local_contraction: next.contraction_factor(),
        right_edge_ratio,
        right_cavity_ratio,
    };

    Ok((transport, factorization))
}
