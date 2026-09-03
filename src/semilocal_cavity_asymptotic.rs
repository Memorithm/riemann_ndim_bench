//! Formal large-degree asymptotics for the row-frozen semilocal cavity model.
//!
//! The formulas in this module are derived from the exact `K(0)` coefficients
//! and the positive frozen cavity fixed point. They are local asymptotic
//! statements only: no uniform finite-section error bound or global trace
//! asymptotic is asserted here.

use std::f64::consts::PI;

use crate::semilocal::ProlateParity;
use crate::semilocal_frozen_cavity::frozen_row_cavity_fixed_point;
use crate::semilocal_resolvent::ResolventTraceError;

/// Diagnostic comparing the exact frozen contraction with its second-order
/// large-degree asymptotic expansion.
#[derive(Clone, Copy, Debug)]
pub struct FrozenContractionAsymptoticDiagnostic {
    degree: usize,
    exact_contraction: f64,
    second_order_contraction: f64,
}

impl FrozenContractionAsymptoticDiagnostic {
    /// Parity degree `d = 2 i + epsilon` associated with the frozen row.
    #[inline]
    pub fn degree(self) -> usize {
        self.degree
    }

    /// Exact frozen cavity contraction `kappa_d = b_d^2 / q_d^2`.
    #[inline]
    pub fn exact_contraction(self) -> f64 {
        self.exact_contraction
    }

    /// Formal second-order approximation `1 - 1/d + 3/(4 d^2)`.
    #[inline]
    pub fn second_order_contraction(self) -> f64 {
        self.second_order_contraction
    }

    /// Exact minus second-order approximation.
    #[inline]
    pub fn residual(self) -> f64 {
        self.exact_contraction - self.second_order_contraction
    }

    /// Scaled first-order gap `d (1-kappa_d)`, formally tending to `1`.
    #[inline]
    pub fn scaled_first_order_gap(self) -> f64 {
        self.degree as f64 * (1.0 - self.exact_contraction)
    }

    /// Scaled second-order remainder
    /// `d^2 (kappa_d - 1 + 1/d)`, formally tending to `3/4`.
    #[inline]
    pub fn scaled_second_order_remainder(self) -> f64 {
        let degree = self.degree as f64;
        degree * degree * (self.exact_contraction - 1.0 + 1.0 / degree)
    }
}

/// Formal expansion of the frozen diagonal coefficient through `d^-2`.
///
/// `a_d = d/(4 pi) + 1/(16 pi) + 5/(64 pi d) - 5/(256 pi d^2)`.
pub fn frozen_diagonal_second_order(degree: usize) -> Option<f64> {
    let degree = nonzero_degree(degree)?;
    Some(
        degree / (4.0 * PI) + 1.0 / (16.0 * PI) + 5.0 / (64.0 * PI * degree)
            - 5.0 / (256.0 * PI * degree * degree),
    )
}

/// Formal expansion of the arithmetic-mean frozen edge through `d^-2`.
///
/// `b_d = d/(8 pi) + 1/(32 pi) + 3/(128 pi d) - 3/(512 pi d^2)`.
pub fn frozen_edge_second_order(degree: usize) -> Option<f64> {
    let degree = nonzero_degree(degree)?;
    Some(
        degree / (8.0 * PI) + 1.0 / (32.0 * PI) + 3.0 / (128.0 * PI * degree)
            - 3.0 / (512.0 * PI * degree * degree),
    )
}

/// Formal soft-edge symbol gap `a_d - 2 b_d` through `d^-2`.
///
/// `a_d - 2 b_d = 1/(32 pi d) - 1/(128 pi d^2)`.
pub fn frozen_soft_gap_second_order(degree: usize) -> Option<f64> {
    let degree = nonzero_degree(degree)?;
    Some(1.0 / (32.0 * PI * degree) - 1.0 / (128.0 * PI * degree * degree))
}

/// Formal expansion of the positive frozen cavity denominator through `d^-2`.
///
/// `q_d = d/(8 pi) + 3/(32 pi) + 5/(128 pi d) - 49/(512 pi d^2)`.
pub fn frozen_cavity_denominator_second_order(degree: usize) -> Option<f64> {
    let degree = nonzero_degree(degree)?;
    Some(
        degree / (8.0 * PI) + 3.0 / (32.0 * PI) + 5.0 / (128.0 * PI * degree)
            - 49.0 / (512.0 * PI * degree * degree),
    )
}

/// Formal second-order frozen contraction.
///
/// `kappa_d = 1 - 1/d + 3/(4 d^2) + O(d^-3)`.
pub fn frozen_contraction_second_order(degree: usize) -> Option<f64> {
    let degree = nonzero_degree(degree)?;
    Some(1.0 - 1.0 / degree + 3.0 / (4.0 * degree * degree))
}

/// Compare the exact frozen contraction at one row with the formal large-degree
/// approximation.
pub fn frozen_contraction_asymptotic_diagnostic(
    block_size: usize,
    row: usize,
    parity: ProlateParity,
) -> Result<FrozenContractionAsymptoticDiagnostic, ResolventTraceError> {
    let degree = parity.degree(row);
    let fixed_point = frozen_row_cavity_fixed_point(block_size, row, parity, 0.0)?;
    let second_order_contraction = frozen_contraction_second_order(degree)
        .expect("an admissible frozen interior row has positive parity degree");

    Ok(FrozenContractionAsymptoticDiagnostic {
        degree,
        exact_contraction: fixed_point.contraction_factor(),
        second_order_contraction,
    })
}

#[inline]
fn nonzero_degree(degree: usize) -> Option<f64> {
    (degree > 0).then_some(degree as f64)
}
