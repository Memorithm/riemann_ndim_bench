//! Source-locked Fourier multiplier attached to the semilocal Poisson map.
//!
//! Connes--Consani, *The Scaling Hamiltonian*, Section 4.1, explains that the
//! semilocal Poisson map
//!
//! `E(f)(x) = |x|^(1/2) sum_{m in M_S} f(m x)`
//!
//! is formally read in multiplicative Fourier variables through
//!
//! `sum_{m in M_S} m^(-1/2+it)
//!    = zeta(1/2-it) product_{p in S_f} (1-p^(-1/2+it))`.
//!
//! The source explicitly presents this discussion at a formal/heuristic level
//! on the critical line. This module therefore separates two statuses:
//!
//! - the finite Euler-deletion factor is an exact finite algebraic object for
//!   any finite complex exponent;
//! - the Dirichlet-series identity is numerically certified here only in the
//!   absolutely convergent half-plane `Re(s) > 1`.
//!
//! Nothing in this module promotes the critical-line formal identity to a
//! convergence theorem or proves a Weil-positivity statement.

use std::fmt;

use num_complex::Complex64;

use crate::semilocal_poisson::SemilocalPoissonMonoid;
use crate::semilocal_trace_contract::FinitePlaceSet;

/// Error returned by the semilocal Fourier/Dirichlet helpers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SemilocalFourierMultiplierError {
    /// A complex exponent must have finite real and imaginary parts.
    InvalidExponent { re: f64, im: f64 },
    /// Direct Dirichlet summation is accepted only for a real exponent > 1.
    NonConvergentRealExponent { sigma: f64 },
    /// A finite prefix must contain at least the first positive integer.
    EmptyPrefix,
}

impl fmt::Display for SemilocalFourierMultiplierError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidExponent { re, im } => {
                write!(f, "Dirichlet exponent must be finite: re={re}, im={im}")
            }
            Self::NonConvergentRealExponent { sigma } => write!(
                f,
                "direct semilocal Dirichlet certification requires sigma > 1: {sigma}"
            ),
            Self::EmptyPrefix => write!(f, "Dirichlet prefix bound must be at least 1"),
        }
    }
}

impl std::error::Error for SemilocalFourierMultiplierError {}

/// Exponent used by the source's critical-line notation.
///
/// The source writes `m^(-1/2+it) = m^(-s)` with
/// `s = 1/2 - i t`.
pub fn critical_line_dirichlet_exponent(
    t: f64,
) -> Result<Complex64, SemilocalFourierMultiplierError> {
    if !t.is_finite() {
        return Err(SemilocalFourierMultiplierError::InvalidExponent {
            re: 0.5,
            im: -t,
        });
    }
    Ok(Complex64::new(0.5, -t))
}

/// Compute `p^(-s)` without relying on a branch of complex logarithm.
fn positive_integer_negative_power(p: u64, s: Complex64) -> Complex64 {
    let log_p = (p as f64).ln();
    let magnitude = (-s.re * log_p).exp();
    let angle = -s.im * log_p;
    Complex64::new(magnitude * angle.cos(), magnitude * angle.sin())
}

fn checked_exponent(
    exponent: Complex64,
) -> Result<Complex64, SemilocalFourierMultiplierError> {
    if exponent.re.is_finite() && exponent.im.is_finite() {
        Ok(exponent)
    } else {
        Err(SemilocalFourierMultiplierError::InvalidExponent {
            re: exponent.re,
            im: exponent.im,
        })
    }
}

/// Exact finite Euler factor deleting the finite places in `S`:
///
/// `product_{p in S_f} (1 - p^(-s))`.
///
/// This is a finite product and therefore does not depend on any convergence
/// claim for the zeta Dirichlet series.
pub fn finite_euler_deletion_factor(
    places: &FinitePlaceSet,
    exponent: Complex64,
) -> Result<Complex64, SemilocalFourierMultiplierError> {
    let exponent = checked_exponent(exponent)?;
    let mut product = Complex64::new(1.0, 0.0);
    for &prime in places.finite_primes() {
        product *= Complex64::new(1.0, 0.0)
            - positive_integer_negative_power(prime, exponent);
    }
    Ok(product)
}

/// Multiply a caller-supplied `zeta(s)` value by the exact finite Euler
/// deletion factor associated with the semilocal place set.
///
/// Supplying `zeta(s)` is deliberately separate from this module. In
/// particular, calling this function on the critical line does not assert that
/// the defining Dirichlet series converges there.
pub fn semilocal_multiplier_from_zeta(
    zeta_value: Complex64,
    places: &FinitePlaceSet,
    exponent: Complex64,
) -> Result<Complex64, SemilocalFourierMultiplierError> {
    let factor = finite_euler_deletion_factor(places, exponent)?;
    Ok(zeta_value * factor)
}

/// Auditable finite prefix of the absolutely convergent semilocal Dirichlet
/// series `sum_{m in M_S} m^(-sigma)` for `sigma > 1`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConvergentDirichletPrefix {
    value: f64,
    max_m: u64,
    term_count: usize,
    tail_upper_bound: f64,
}

impl ConvergentDirichletPrefix {
    #[inline]
    pub fn value(self) -> f64 {
        self.value
    }

    #[inline]
    pub fn max_m(self) -> u64 {
        self.max_m
    }

    #[inline]
    pub fn term_count(self) -> usize {
        self.term_count
    }

    /// Rigorous elementary upper bound obtained by dropping the `M_S`
    /// restriction and applying the integral test to the remaining positive
    /// integers.
    #[inline]
    pub fn tail_upper_bound(self) -> f64 {
        self.tail_upper_bound
    }
}

/// Sum the semilocal Dirichlet series through `max_m` in its absolutely
/// convergent real half-plane.
///
/// The omitted tail is bounded by
///
/// `sum_{n>max_m} n^(-sigma) <= max_m^(1-sigma)/(sigma-1)`.
pub fn convergent_dirichlet_prefix(
    monoid: &SemilocalPoissonMonoid,
    sigma: f64,
    max_m: u64,
) -> Result<ConvergentDirichletPrefix, SemilocalFourierMultiplierError> {
    if !sigma.is_finite() || sigma <= 1.0 {
        return Err(SemilocalFourierMultiplierError::NonConvergentRealExponent {
            sigma,
        });
    }
    if max_m == 0 {
        return Err(SemilocalFourierMultiplierError::EmptyPrefix);
    }

    let mut value = 0.0;
    let mut term_count = 0_usize;
    for m in 1..=max_m {
        if monoid.contains(m) {
            value += (m as f64).powf(-sigma);
            term_count += 1;
        }
    }

    let tail_upper_bound =
        (max_m as f64).powf(1.0 - sigma) / (sigma - 1.0);

    Ok(ConvergentDirichletPrefix {
        value,
        max_m,
        term_count,
        tail_upper_bound,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_critical_exponent_has_expected_sign() {
        let exponent = critical_line_dirichlet_exponent(3.0).unwrap();
        assert_eq!(exponent, Complex64::new(0.5, -3.0));
    }

    #[test]
    fn finite_euler_factor_is_exact_at_real_sigma_two() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let factor = finite_euler_deletion_factor(
            &places,
            Complex64::new(2.0, 0.0),
        )
        .unwrap();
        assert!((factor.re - 2.0 / 3.0).abs() < 2.0e-15);
        assert!(factor.im.abs() < 2.0e-15);
    }

    #[test]
    fn direct_series_refuses_the_critical_line() {
        let monoid = SemilocalPoissonMonoid::new(FinitePlaceSet::new(vec![2]).unwrap());
        assert_eq!(
            convergent_dirichlet_prefix(&monoid, 0.5, 100).unwrap_err(),
            SemilocalFourierMultiplierError::NonConvergentRealExponent { sigma: 0.5 }
        );
    }
}
