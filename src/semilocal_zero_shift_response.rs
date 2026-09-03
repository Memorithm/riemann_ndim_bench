//! Exact first shift derivative of the semilocal left Schur cavity at zero.
//!
//! The zero-shift left cavity has an exact closed form. Differentiating the
//! finite Schur recurrence shows that its relative first shift derivative is
//! the same at every row and in both parity sectors. This is a first-order
//! identity only; no uniform higher-order remainder is asserted here.

use std::f64::consts::PI;

use crate::semilocal::ProlateParity;
use crate::semilocal_zero_shift_cavity::zero_shift_left_cavity_closed_form;

/// Exact zero-shift left-cavity value and its first derivative with respect to
/// the resolvent shift.
#[derive(Clone, Copy, Debug)]
pub struct ZeroShiftLeftCavityResponse {
    degree: usize,
    denominator: f64,
    shift_derivative: f64,
}

impl ZeroShiftLeftCavityResponse {
    /// Parity degree `d = 2 i + epsilon`.
    #[inline]
    pub fn degree(self) -> usize {
        self.degree
    }

    /// Exact `L_i(0)`.
    #[inline]
    pub fn denominator(self) -> f64 {
        self.denominator
    }

    /// Exact derivative `L_i'(0)` with respect to the additive shift in
    /// `K(0) + t I`.
    #[inline]
    pub fn shift_derivative(self) -> f64 {
        self.shift_derivative
    }

    /// Exact logarithmic derivative `L_i'(0)/L_i(0) = 8 pi / 3`.
    #[inline]
    pub fn relative_shift_derivative(self) -> f64 {
        self.shift_derivative / self.denominator
    }
}

/// Exact first shift response of one left cavity row.
///
/// If
///
/// `L_0(t) = K_00 + t`
///
/// and
///
/// `L_i(t) = K_ii + t - K_i,i-1^2/L_{i-1}(t)`,
///
/// then
///
/// `L_i'(0) = (8 pi / 3) L_i(0)`
///
/// for every row and both parity sectors. Equivalently,
///
/// `d/dt log L_i(t)|_{t=0} = 8 pi / 3`.
///
/// The closed derivative can also be written directly in the parity degree:
///
/// `L_i'(0) = (2d+1)(2d+3) / [3(4d+1)]`.
pub fn zero_shift_left_cavity_response(
    row: usize,
    parity: ProlateParity,
) -> ZeroShiftLeftCavityResponse {
    let closed = zero_shift_left_cavity_closed_form(row, parity);
    let shift_derivative = (8.0 * PI / 3.0) * closed.denominator();

    ZeroShiftLeftCavityResponse {
        degree: closed.degree(),
        denominator: closed.denominator(),
        shift_derivative,
    }
}

/// Generate exact first shift derivatives for a finite left-cavity prefix.
pub fn zero_shift_left_cavity_shift_derivatives(
    block_size: usize,
    parity: ProlateParity,
) -> Vec<f64> {
    (0..block_size)
        .map(|row| zero_shift_left_cavity_response(row, parity).shift_derivative())
        .collect()
}

/// Universal exact relative first shift derivative of the zero-shift left
/// cavity.
#[inline]
pub fn zero_shift_left_cavity_relative_derivative() -> f64 {
    8.0 * PI / 3.0
}
