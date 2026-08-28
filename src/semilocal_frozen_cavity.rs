//! Frozen Toeplitz cavity fixed point for the semilocal resolvent model.
//!
//! This module bridges the exact finite left/right cavity recurrences with the
//! row-frozen infinite Toeplitz model.  It only describes the frozen local
//! problem; it does not assert an error bound for the variable-coefficient
//! finite section.

use crate::semilocal::{ProlateParity, build_k0, build_kprime_closed};
use crate::semilocal_resolvent::ResolventTraceError;

/// Positive cavity fixed point and the corresponding infinite-Toeplitz Green
/// data for one frozen interior row.
#[derive(Clone, Copy, Debug)]
pub struct FrozenCavityFixedPoint {
    diagonal_coefficient: f64,
    edge_coefficient: f64,
    weight_diagonal: f64,
    weight_edge: f64,
    cavity_denominator: f64,
    green_diagonal: f64,
    green_off_diagonal: f64,
    contraction_factor: f64,
    weighted_density: f64,
}

impl FrozenCavityFixedPoint {
    /// Frozen shifted diagonal coefficient `a`.
    #[inline]
    pub fn diagonal_coefficient(self) -> f64 {
        self.diagonal_coefficient
    }

    /// Arithmetic mean `b` of the two neighboring `K(0)` edges.
    #[inline]
    pub fn edge_coefficient(self) -> f64 {
        self.edge_coefficient
    }

    /// Frozen sign-corrected perturbation diagonal `d`.
    #[inline]
    pub fn weight_diagonal(self) -> f64 {
        self.weight_diagonal
    }

    /// Arithmetic mean `o` of the neighboring sign-corrected perturbation edges.
    #[inline]
    pub fn weight_edge(self) -> f64 {
        self.weight_edge
    }

    /// Positive solution `q` of `q = a - b^2 / q`.
    #[inline]
    pub fn cavity_denominator(self) -> f64 {
        self.cavity_denominator
    }

    /// Infinite-Toeplitz diagonal Green entry.
    #[inline]
    pub fn green_diagonal(self) -> f64 {
        self.green_diagonal
    }

    /// Infinite-Toeplitz first off-diagonal Green entry in the original basis.
    #[inline]
    pub fn green_off_diagonal(self) -> f64 {
        self.green_off_diagonal
    }

    /// Derivative of the frozen cavity map at the positive fixed point.
    ///
    /// For `F(x)=a-b^2/x`, this is `F'(q)=b^2/q^2`.
    #[inline]
    pub fn contraction_factor(self) -> f64 {
        self.contraction_factor
    }

    /// Frozen weighted resolvent density `d G_00 + 2 o G_01`.
    #[inline]
    pub fn weighted_density(self) -> f64 {
        self.weighted_density
    }
}

/// Construct the positive frozen cavity fixed point for one interior row.
///
/// With frozen shifted diagonal `a` and edge `b`, the half-line cavity map is
///
/// `F(x) = a - b^2 / x`.
///
/// Strict positivity of the frozen Toeplitz symbol gives
/// `a^2 - 4 b^2 > 0`, and the positive stable fixed point is
///
/// `q = (a + sqrt(a^2 - 4 b^2)) / 2`.
///
/// The corresponding two-sided Green bands are
///
/// `G_00 = 1 / sqrt(a^2 - 4 b^2)`
///
/// and
///
/// `G_01 = -b G_00 / q`.
///
/// No finite-section convergence estimate is asserted here.
pub fn frozen_row_cavity_fixed_point(
    block_size: usize,
    row: usize,
    parity: ProlateParity,
    shift: f64,
) -> Result<FrozenCavityFixedPoint, ResolventTraceError> {
    if !shift.is_finite() || shift < 0.0 {
        return Err(ResolventTraceError::InvalidShift { value: shift });
    }
    if row == 0 || row >= block_size.saturating_sub(1) {
        return Err(ResolventTraceError::InvalidInteriorRow {
            row,
            len: block_size,
        });
    }

    let k0 = build_k0(block_size, parity);
    let kprime = build_kprime_closed(block_size, parity);
    let sign = parity.sign_correction();

    let diagonal_coefficient = k0.diagonal()[row] + shift;
    let edge_coefficient = 0.5 * (k0.off_diagonal()[row - 1] + k0.off_diagonal()[row]);
    let weight_diagonal = sign * kprime.diagonal()[row];
    let weight_edge = 0.5 * sign * (kprime.off_diagonal()[row - 1] + kprime.off_diagonal()[row]);

    let twice_edge = 2.0 * edge_coefficient;
    let discriminant = (diagonal_coefficient - twice_edge) * (diagonal_coefficient + twice_edge);
    if !discriminant.is_finite() || discriminant <= 0.0 {
        return Err(ResolventTraceError::NonPositiveFrozenDiscriminant {
            row,
            value: discriminant,
        });
    }

    let square_root = discriminant.sqrt();
    let cavity_denominator = 0.5 * (diagonal_coefficient + square_root);
    let green_diagonal = 1.0 / square_root;
    let green_off_diagonal = -edge_coefficient * green_diagonal / cavity_denominator;
    let edge_ratio = edge_coefficient / cavity_denominator;
    let contraction_factor = edge_ratio * edge_ratio;
    let weighted_density =
        weight_diagonal * green_diagonal + 2.0 * weight_edge * green_off_diagonal;

    Ok(FrozenCavityFixedPoint {
        diagonal_coefficient,
        edge_coefficient,
        weight_diagonal,
        weight_edge,
        cavity_denominator,
        green_diagonal,
        green_off_diagonal,
        contraction_factor,
        weighted_density,
    })
}
