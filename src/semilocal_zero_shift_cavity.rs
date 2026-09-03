//! Exact zero-shift left Schur cavity for the semilocal `K(0)` blocks.
//!
//! The left-to-right cavity recurrence admits a closed form inherited from the
//! exact ground-state factorisation. This is a finite algebraic identity at
//! shift zero, not a large-degree approximation.

use std::f64::consts::PI;

use crate::semilocal::ProlateParity;

/// Exact closed-form data for one zero-shift left cavity row.
#[derive(Clone, Copy, Debug)]
pub struct ZeroShiftLeftCavityClosedForm {
    degree: usize,
    denominator: f64,
    incoming_schur_correction: f64,
}

impl ZeroShiftLeftCavityClosedForm {
    /// Parity degree `d = 2 i + epsilon`.
    #[inline]
    pub fn degree(self) -> usize {
        self.degree
    }

    /// Exact zero-shift left Schur denominator
    /// `(d+1/2)(d+3/2) / [2 pi (4d+1)]`.
    #[inline]
    pub fn denominator(self) -> f64 {
        self.denominator
    }

    /// Exact incoming Schur correction `d(d-1) / [2 pi (4d+1)]`.
    ///
    /// This is zero at the left endpoint in both parity sectors.
    #[inline]
    pub fn incoming_schur_correction(self) -> f64 {
        self.incoming_schur_correction
    }

    /// Reconstruct the exact `K(0)` diagonal as left denominator plus incoming
    /// Schur correction.
    #[inline]
    pub fn reconstructed_diagonal(self) -> f64 {
        self.denominator + self.incoming_schur_correction
    }
}

/// Exact zero-shift left cavity at one parity row.
///
/// For parity degree `d = 2 i + epsilon`, let `L_i` satisfy
///
/// `L_0 = K_00`,
///
/// `L_i = K_ii - K_i,i-1^2 / L_{i-1}`.
///
/// Then exactly
///
/// `L_i = (d+1/2)(d+3/2) / [2 pi (4d+1)]`.
///
/// For `i>0`, substituting the previous closed form into the incoming edge
/// gives
///
/// `K_i,i-1^2 / L_{i-1} = d(d-1) / [2 pi (4d+1)]`,
///
/// which closes the induction. The formula also holds at the left endpoint,
/// where the incoming correction is absent and evaluates to zero.
pub fn zero_shift_left_cavity_closed_form(
    row: usize,
    parity: ProlateParity,
) -> ZeroShiftLeftCavityClosedForm {
    let degree = parity.degree(row);
    let degree_f = degree as f64;
    let normalization = 2.0 * PI * (4.0 * degree_f + 1.0);
    let denominator = (degree_f + 0.5) * (degree_f + 1.5) / normalization;
    let incoming_schur_correction = degree_f * (degree_f - 1.0) / normalization;

    ZeroShiftLeftCavityClosedForm {
        degree,
        denominator,
        incoming_schur_correction,
    }
}

/// Generate the exact zero-shift left cavity prefix for one finite parity
/// block. The values are independent of the block's right endpoint; the block
/// size only chooses how many rows are returned.
pub fn zero_shift_left_cavity_denominators(block_size: usize, parity: ProlateParity) -> Vec<f64> {
    (0..block_size)
        .map(|row| zero_shift_left_cavity_closed_form(row, parity).denominator())
        .collect()
}
