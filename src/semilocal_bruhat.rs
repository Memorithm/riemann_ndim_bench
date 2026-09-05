//! Elementary semilocal Bruhat--Schwartz data on the diagonal copy of `Q_S`.
//!
//! This module deliberately models only the finite local factor
//!
//! `prod_{p in S_f} 1_{Z_p}`
//!
//! evaluated on the diagonal embedding of `Q_S`.  It does **not** model a
//! general adele in `A_S`, the quotient `X_S`, or the semilocal Fourier
//! transform.  The purpose is to make the finite-place content of manufactured
//! Poisson fixtures explicit instead of hiding it behind a scalar callback.

use crate::semilocal_qs::QsRational;
use crate::semilocal_trace_contract::FinitePlaceSet;

/// Standard elementary finite Bruhat--Schwartz factor
/// `prod_{p in S_f} 1_{Z_p}` on the diagonal copy of `Q_S`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ElementaryFiniteBruhatFactor {
    places: FinitePlaceSet,
}

impl ElementaryFiniteBruhatFactor {
    /// Build the standard finite factor for the declared semilocal place set.
    pub fn new(places: FinitePlaceSet) -> Self {
        Self { places }
    }

    /// Finite places carried by this factor.
    #[inline]
    pub fn places(&self) -> &FinitePlaceSet {
        &self.places
    }

    /// Exact p-adic valuation of a diagonal `Q_S` rational.
    ///
    /// Zero is treated separately because its p-adic valuation is `+infinity`.
    /// For non-zero `q=a/b` in reduced form this returns `v_p(a)-v_p(b)`.
    pub fn diagonal_valuation(&self, q: QsRational, prime: u64) -> Option<i32> {
        if !self.places.contains_prime(prime) || q.is_zero() {
            return None;
        }
        Some(
            valuation_u64(q.numerator_magnitude(), prime) as i32
                - valuation_u64(q.denominator(), prime) as i32,
        )
    }

    /// Evaluate `1_{Z_p}` on the diagonal sample at one finite place.
    ///
    /// Returns `None` when `prime` is not one of the declared finite places.
    pub fn local_unit_ball_indicator(&self, q: QsRational, prime: u64) -> Option<u8> {
        if !self.places.contains_prime(prime) {
            return None;
        }
        if q.is_zero() {
            return Some(1);
        }
        Some((self.diagonal_valuation(q, prime)? >= 0) as u8)
    }

    /// Evaluate the product finite factor `prod_p 1_{Z_p}`.
    pub fn evaluate_diagonal(&self, q: QsRational) -> u8 {
        if q.is_zero() {
            return 1;
        }
        self.places
            .finite_primes()
            .iter()
            .all(|&prime| self.diagonal_valuation(q, prime).is_some_and(|v| v >= 0)) as u8
    }

    /// Combine an externally supplied archimedean value with the exact finite
    /// factor on a diagonal `Q_S` sample.
    ///
    /// This is a restricted factorizable fixture, not an evaluator on general
    /// `A_S`: the archimedean coordinate is supplied independently by the caller
    /// while the finite coordinates are the diagonal image of `q`.
    #[inline]
    pub fn evaluate_factorizable_diagonal(&self, q: QsRational, archimedean_value: f64) -> f64 {
        f64::from(self.evaluate_diagonal(q)) * archimedean_value
    }
}

fn valuation_u64(mut value: u64, prime: u64) -> u32 {
    debug_assert!(prime >= 2);
    let mut valuation = 0_u32;
    while value != 0 && value.is_multiple_of(prime) {
        value /= prime;
        valuation += 1;
    }
    valuation
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonal_valuations_are_exact() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let factor = ElementaryFiniteBruhatFactor::new(places.clone());
        let q = QsRational::new(45, 8, &places).unwrap();

        assert_eq!(factor.diagonal_valuation(q, 2), Some(-3));
        assert_eq!(factor.diagonal_valuation(q, 3), Some(2));
        assert_eq!(factor.diagonal_valuation(q, 5), None);
    }

    #[test]
    fn zero_belongs_to_every_local_unit_ball() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let factor = ElementaryFiniteBruhatFactor::new(places.clone());
        let zero = QsRational::new(0, 1, &places).unwrap();

        assert_eq!(factor.local_unit_ball_indicator(zero, 2), Some(1));
        assert_eq!(factor.local_unit_ball_indicator(zero, 3), Some(1));
        assert_eq!(factor.evaluate_diagonal(zero), 1);
    }

    #[test]
    fn product_indicator_detects_negative_finite_valuation() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let factor = ElementaryFiniteBruhatFactor::new(places.clone());

        let integral = QsRational::new(18, 1, &places).unwrap();
        let non_integral = QsRational::new(9, 2, &places).unwrap();

        assert_eq!(factor.evaluate_diagonal(integral), 1);
        assert_eq!(factor.evaluate_diagonal(non_integral), 0);
    }

    #[test]
    fn multiplication_by_monoid_elements_preserves_local_valuations() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let factor = ElementaryFiniteBruhatFactor::new(places.clone());

        // 5 is prime to every finite place in S, so multiplying by 5 leaves
        // all p-adic valuations at p=2,3 unchanged.
        let q = QsRational::new(3, 4, &places).unwrap();
        let five_q = QsRational::new(15, 4, &places).unwrap();

        assert_eq!(
            factor.diagonal_valuation(q, 2),
            factor.diagonal_valuation(five_q, 2)
        );
        assert_eq!(
            factor.diagonal_valuation(q, 3),
            factor.diagonal_valuation(five_q, 3)
        );
        assert_eq!(
            factor.evaluate_diagonal(q),
            factor.evaluate_diagonal(five_q)
        );
    }
}
