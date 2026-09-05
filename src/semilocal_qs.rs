//! Exact arithmetic for the semilocal ring `Q_S` and its unit/monoid split.
//!
//! Connes--Consani, *The Scaling Hamiltonian*, equations (2.13)--(2.14), define
//! `Q_S` as the subring of rational numbers whose denominator uses only the
//! finite primes in `S`, and
//!
//! `Q_S^* = { +/- product_j p_j^(n_j) : n_j in Z }`.
//!
//! Immediately before equation (4.6), the source uses the unique decomposition
//!
//! `q = u m`,
//!
//! with `u in Q_S^*` and `m in M_S`, where `M_S` is the monoid of positive
//! integers prime to all finite places in `S`.
//!
//! This module makes that finite arithmetic statement executable. It does not
//! model the adele ring `A_S`, the quotient `X_S`, or the semilocal Fourier
//! transform.

use std::fmt;

use crate::semilocal_trace_contract::FinitePlaceSet;

/// Canonical rational element of `Q_S`.
///
/// The denominator is always positive and the numerator magnitude and
/// denominator are coprime. Zero is represented canonically as `0/1`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct QsRational {
    sign: i8,
    numerator_magnitude: u64,
    denominator: u64,
}

impl QsRational {
    /// Construct and reduce a rational, rejecting denominators with prime
    /// factors outside the declared finite place set.
    pub fn new(
        numerator: i64,
        denominator: u64,
        places: &FinitePlaceSet,
    ) -> Result<Self, QsArithmeticError> {
        if denominator == 0 {
            return Err(QsArithmeticError::ZeroDenominator);
        }
        if numerator == 0 {
            return Ok(Self {
                sign: 0,
                numerator_magnitude: 0,
                denominator: 1,
            });
        }

        let sign = if numerator < 0 { -1 } else { 1 };
        let magnitude = numerator.unsigned_abs();
        let divisor = gcd(magnitude, denominator);
        let reduced_numerator = magnitude / divisor;
        let reduced_denominator = denominator / divisor;
        validate_denominator(reduced_denominator, places)?;

        Ok(Self {
            sign,
            numerator_magnitude: reduced_numerator,
            denominator: reduced_denominator,
        })
    }

    #[inline]
    pub fn sign(self) -> i8 {
        self.sign
    }

    #[inline]
    pub fn numerator_magnitude(self) -> u64 {
        self.numerator_magnitude
    }

    #[inline]
    pub fn denominator(self) -> u64 {
        self.denominator
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.sign == 0
    }

    /// Decompose a non-zero element uniquely as `q = u m`, with
    /// `u in Q_S^*` and `m in M_S`.
    pub fn unit_monoid_decomposition(
        self,
        places: &FinitePlaceSet,
    ) -> Result<QsUnitMonoidDecomposition, QsArithmeticError> {
        if self.is_zero() {
            return Err(QsArithmeticError::ZeroHasNoUnitMonoidDecomposition);
        }
        validate_denominator(self.denominator, places)?;

        let mut numerator = self.numerator_magnitude;
        let mut denominator = self.denominator;
        let mut unit_exponents = Vec::with_capacity(places.finite_primes().len());

        for &prime in places.finite_primes() {
            let (numerator_valuation, stripped_numerator) = strip_prime(numerator, prime);
            let (denominator_valuation, stripped_denominator) = strip_prime(denominator, prime);
            numerator = stripped_numerator;
            denominator = stripped_denominator;
            unit_exponents.push((
                prime,
                numerator_valuation as i32 - denominator_valuation as i32,
            ));
        }

        debug_assert_eq!(denominator, 1);
        debug_assert!(places
            .finite_primes()
            .iter()
            .all(|&prime| !numerator.is_multiple_of(prime)));

        Ok(QsUnitMonoidDecomposition {
            unit_sign: self.sign,
            unit_exponents,
            monoid_element: numerator,
        })
    }
}

/// Canonical `Q_S^* x M_S` decomposition of a non-zero element of `Q_S`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QsUnitMonoidDecomposition {
    unit_sign: i8,
    unit_exponents: Vec<(u64, i32)>,
    monoid_element: u64,
}

impl QsUnitMonoidDecomposition {
    /// Sign of the unit factor (`+1` or `-1`).
    #[inline]
    pub fn unit_sign(&self) -> i8 {
        self.unit_sign
    }

    /// Prime/exponent signature of the `Q_S^*` unit in increasing prime order.
    #[inline]
    pub fn unit_exponents(&self) -> &[(u64, i32)] {
        &self.unit_exponents
    }

    /// Positive `M_S` factor.
    #[inline]
    pub fn monoid_element(&self) -> u64 {
        self.monoid_element
    }

    /// Exponent of a finite-place prime in the unit factor.
    pub fn exponent_for(&self, prime: u64) -> Option<i32> {
        self.unit_exponents
            .binary_search_by_key(&prime, |&(candidate, _)| candidate)
            .ok()
            .map(|index| self.unit_exponents[index].1)
    }

    /// Recompose the exact reduced rational represented by this decomposition.
    ///
    /// This method is primarily an audit hook. The decomposition can only be
    /// constructed from an existing `QsRational`, so the products are bounded
    /// by the original `u64` numerator and denominator.
    pub fn recompose(&self) -> QsRational {
        let mut numerator = self.monoid_element;
        let mut denominator = 1_u64;

        for &(prime, exponent) in &self.unit_exponents {
            if exponent > 0 {
                for _ in 0..exponent {
                    numerator = numerator
                        .checked_mul(prime)
                        .expect("Q_S decomposition numerator must recompose without overflow");
                }
            } else {
                for _ in 0..exponent.unsigned_abs() {
                    denominator = denominator
                        .checked_mul(prime)
                        .expect("Q_S decomposition denominator must recompose without overflow");
                }
            }
        }

        QsRational {
            sign: self.unit_sign,
            numerator_magnitude: numerator,
            denominator,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum QsArithmeticError {
    ZeroDenominator,
    DenominatorOutsidePlaceSet { residual: u64 },
    ZeroHasNoUnitMonoidDecomposition,
}

impl fmt::Display for QsArithmeticError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroDenominator => write!(f, "Q_S rational denominator cannot be zero"),
            Self::DenominatorOutsidePlaceSet { residual } => write!(
                f,
                "Q_S denominator contains a factor outside the finite place set: residual={residual}"
            ),
            Self::ZeroHasNoUnitMonoidDecomposition => write!(
                f,
                "zero cannot be written as a unit of Q_S times a positive M_S element"
            ),
        }
    }
}

impl std::error::Error for QsArithmeticError {}

fn validate_denominator(
    denominator: u64,
    places: &FinitePlaceSet,
) -> Result<(), QsArithmeticError> {
    let mut residual = denominator;
    for &prime in places.finite_primes() {
        while residual.is_multiple_of(prime) {
            residual /= prime;
        }
    }
    if residual == 1 {
        Ok(())
    } else {
        Err(QsArithmeticError::DenominatorOutsidePlaceSet { residual })
    }
}

fn strip_prime(mut value: u64, prime: u64) -> (u32, u64) {
    let mut valuation = 0_u32;
    while value.is_multiple_of(prime) {
        value /= prime;
        valuation += 1;
    }
    (valuation, value)
}

fn gcd(mut left: u64, mut right: u64) -> u64 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructor_reduces_before_checking_localized_denominator() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let rational = QsRational::new(-150, 72, &places).unwrap();
        assert_eq!(rational.sign(), -1);
        assert_eq!(rational.numerator_magnitude(), 25);
        assert_eq!(rational.denominator(), 12);
    }

    #[test]
    fn denominator_outside_s_is_rejected() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        assert_eq!(
            QsRational::new(1, 10, &places).unwrap_err(),
            QsArithmeticError::DenominatorOutsidePlaceSet { residual: 5 }
        );
    }

    #[test]
    fn unit_monoid_decomposition_recomposes_exactly() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let rational = QsRational::new(45, 8, &places).unwrap();
        let decomposition = rational.unit_monoid_decomposition(&places).unwrap();

        assert_eq!(decomposition.unit_sign(), 1);
        assert_eq!(decomposition.exponent_for(2), Some(-3));
        assert_eq!(decomposition.exponent_for(3), Some(2));
        assert_eq!(decomposition.monoid_element(), 5);
        assert_eq!(decomposition.recompose(), rational);
    }

    #[test]
    fn zero_has_no_source_unit_monoid_split() {
        let places = FinitePlaceSet::new(vec![2]).unwrap();
        let zero = QsRational::new(0, 16, &places).unwrap();
        assert_eq!(zero.denominator(), 1);
        assert_eq!(
            zero.unit_monoid_decomposition(&places).unwrap_err(),
            QsArithmeticError::ZeroHasNoUnitMonoidDecomposition
        );
    }
}
