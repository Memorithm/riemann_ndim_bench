//! Source-locked boundary transform used in the Weil positivity criterion.
//!
//! Connes--Consani, *Weil positivity and Trace formula, the archimedean
//! place*, Lemma 3.3, identifies the ideal of compactly supported
//! multiplicative test functions satisfying the two boundary conditions
//!
//! `integral f(rho) rho^(+1/2) d*rho = 0`
//! `integral f(rho) rho^(-1/2) d*rho = 0`
//!
//! as the range of
//!
//! `Q = -(rho d/drho)^2 + 1/4`.
//!
//! This module only encodes that source-level differential/moment contract.
//! It does not assert Weil positivity and does not prove RH.

use std::fmt;

use crate::quadrature::{GaussLegendreUnit, QuadratureError};

/// Compact support interval in the multiplicative variable `rho > 0`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultiplicativeSupport {
    lower: f64,
    upper: f64,
    log_lower: f64,
    log_upper: f64,
}

impl MultiplicativeSupport {
    /// Construct a finite non-empty multiplicative support interval.
    pub fn new(lower: f64, upper: f64) -> Result<Self, WeilBoundaryError> {
        if !lower.is_finite()
            || !upper.is_finite()
            || lower <= 0.0
            || upper <= lower
        {
            return Err(WeilBoundaryError::InvalidSupport { lower, upper });
        }
        Ok(Self {
            lower,
            upper,
            log_lower: lower.ln(),
            log_upper: upper.ln(),
        })
    }

    #[inline]
    pub fn lower(self) -> f64 {
        self.lower
    }

    #[inline]
    pub fn upper(self) -> f64 {
        self.upper
    }

    #[inline]
    pub fn log_lower(self) -> f64 {
        self.log_lower
    }

    #[inline]
    pub fn log_upper(self) -> f64 {
        self.log_upper
    }

    #[inline]
    pub fn contains(self, rho: f64) -> bool {
        rho >= self.lower && rho <= self.upper
    }
}

/// Error returned by the source-locked boundary-transform helpers.
#[derive(Debug)]
pub enum WeilBoundaryError {
    InvalidSupport { lower: f64, upper: f64 },
    InvalidRho { rho: f64 },
    InvalidExponent { exponent: f64 },
    Quadrature(QuadratureError),
}

impl fmt::Display for WeilBoundaryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidSupport { lower, upper } => {
                write!(f, "invalid multiplicative support ({lower}, {upper})")
            }
            Self::InvalidRho { rho } => write!(f, "rho must be finite and positive: {rho}"),
            Self::InvalidExponent { exponent } => {
                write!(f, "Mellin power exponent must be finite: {exponent}")
            }
            Self::Quadrature(error) => write!(f, "quadrature construction failed: {error:?}"),
        }
    }
}

impl std::error::Error for WeilBoundaryError {}

impl From<QuadratureError> for WeilBoundaryError {
    fn from(value: QuadratureError) -> Self {
        Self::Quadrature(value)
    }
}

/// The multiplier of `Q` on the multiplicative character `rho^exponent`.
///
/// In logarithmic coordinate `x = log rho`, `Q = -d_x^2 + 1/4`, so
/// `Q exp(exponent x) = (1/4 - exponent^2) exp(exponent x)`.
#[inline]
pub fn character_multiplier(exponent: f64) -> f64 {
    0.25 - exponent * exponent
}

/// Apply `Q` from a value and its second derivative in `x = log rho`.
#[inline]
pub fn q_from_log_second_derivative(value: f64, second_log_derivative: f64) -> f64 {
    -second_log_derivative + 0.25 * value
}

/// Apply `Q` from ordinary derivatives with respect to `rho`.
///
/// Since `(rho d_rho)^2 f = rho^2 f'' + rho f'`, this is exactly the same
/// source operator expressed in the multiplicative coordinate.
pub fn q_from_rho_derivatives(
    rho: f64,
    value: f64,
    first_derivative: f64,
    second_derivative: f64,
) -> Result<f64, WeilBoundaryError> {
    checked_rho(rho)?;
    let second_log_derivative = rho * rho * second_derivative + rho * first_derivative;
    Ok(q_from_log_second_derivative(
        value,
        second_log_derivative,
    ))
}

/// Apply `Q` on a declared compact support interval, returning exact zero
/// outside the interval.
///
/// The caller supplies the value and second logarithmic derivative only for
/// points inside the support. This mirrors Lemma 3.3(ii): `Q` does not enlarge
/// support.
pub fn q_on_support(
    support: MultiplicativeSupport,
    rho: f64,
    mut value_and_second_log_derivative: impl FnMut(f64) -> (f64, f64),
) -> Result<f64, WeilBoundaryError> {
    checked_rho(rho)?;
    if !support.contains(rho) {
        return Ok(0.0);
    }
    let (value, second_log_derivative) = value_and_second_log_derivative(rho);
    Ok(q_from_log_second_derivative(
        value,
        second_log_derivative,
    ))
}

/// Multiplicative Mellin-power moment
///
/// `integral f(rho) rho^exponent d*rho`, where `d*rho = d rho / rho`.
///
/// The implementation integrates in logarithmic coordinate, where Haar
/// measure is ordinary Lebesgue measure.
pub fn mellin_power_moment(
    support: MultiplicativeSupport,
    quadrature_order: usize,
    exponent: f64,
    mut value: impl FnMut(f64) -> f64,
) -> Result<f64, WeilBoundaryError> {
    if !exponent.is_finite() {
        return Err(WeilBoundaryError::InvalidExponent { exponent });
    }
    let quadrature = GaussLegendreUnit::new(quadrature_order)?;
    let span = support.log_upper() - support.log_lower();
    let total = quadrature
        .nodes()
        .iter()
        .zip(quadrature.weights().iter())
        .map(|(&unit_x, &weight)| {
            let log_rho = support.log_lower() + span * unit_x;
            let rho = log_rho.exp();
            weight * value(rho) * (exponent * log_rho).exp()
        })
        .sum::<f64>();
    Ok(span * total)
}

/// The two source boundary moments at powers `+1/2` and `-1/2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WeilBoundaryMoments {
    pub plus_half: f64,
    pub minus_half: f64,
}

/// Evaluate the two boundary moments used by the Weil criterion.
pub fn critical_boundary_moments(
    support: MultiplicativeSupport,
    quadrature_order: usize,
    mut value: impl FnMut(f64) -> f64,
) -> Result<WeilBoundaryMoments, WeilBoundaryError> {
    let plus_half = mellin_power_moment(support, quadrature_order, 0.5, &mut value)?;
    let minus_half = mellin_power_moment(support, quadrature_order, -0.5, value)?;
    Ok(WeilBoundaryMoments {
        plus_half,
        minus_half,
    })
}

fn checked_rho(rho: f64) -> Result<(), WeilBoundaryError> {
    if rho.is_finite() && rho > 0.0 {
        Ok(())
    } else {
        Err(WeilBoundaryError::InvalidRho { rho })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boundary_characters_are_exactly_annihilated() {
        assert_eq!(character_multiplier(0.5), 0.0);
        assert_eq!(character_multiplier(-0.5), 0.0);
        assert_eq!(character_multiplier(0.0), 0.25);
    }

    #[test]
    fn rho_and_log_coordinate_forms_agree() {
        let rho = 1.75_f64;
        let value = rho.powi(3);
        let first = 3.0 * rho.powi(2);
        let second = 6.0 * rho;
        let from_rho = q_from_rho_derivatives(rho, value, first, second).unwrap();
        let from_log = q_from_log_second_derivative(value, 9.0 * value);
        assert_eq!(from_rho, from_log);
    }
}
