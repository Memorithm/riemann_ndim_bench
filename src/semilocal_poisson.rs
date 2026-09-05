//! Source arithmetic behind the semilocal Poisson map.
//!
//! Section 4.1 of Connes--Consani, *The Scaling Hamiltonian*, writes every
//! element of `Q_S` uniquely as `u m`, with `u in Q_S^*` and `m` in the monoid
//! `M_S` of positive integers prime to every finite prime in `S`. The source map
//!
//! `E(f)(x) = |x|^(1/2) sum_{m in M_S} f(m x)`
//!
//! is one of the bridges through which Poisson summation produces zeta and the
//! omitted local Euler factors.
//!
//! RiemannBench does not yet represent a general adele `x in A_S`. This module
//! therefore exposes the exact monoid and a finite callback-driven sum. A caller
//! must provide a bound known to contain every non-zero term. No Poisson identity
//! or Weil-positivity statement is asserted here.

use std::fmt;

use crate::semilocal_trace_contract::{FinitePlaceSet, SemilocalSpaceContract};

/// Positive-integer monoid `M_S` attached to a finite semilocal place set.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SemilocalPoissonMonoid {
    places: FinitePlaceSet,
}

impl SemilocalPoissonMonoid {
    /// Build `M_S` from a source-locked semilocal space contract.
    pub fn from_space(space: &SemilocalSpaceContract) -> Self {
        Self {
            places: space.places().clone(),
        }
    }

    /// Build directly from a finite place set.
    pub fn new(places: FinitePlaceSet) -> Self {
        Self { places }
    }

    #[inline]
    pub fn places(&self) -> &FinitePlaceSet {
        &self.places
    }

    /// Membership in `M_S`: positive and not divisible by any finite place.
    pub fn contains(&self, value: u64) -> bool {
        value >= 1
            && self
                .places
                .finite_primes()
                .iter()
                .all(|&prime| !value.is_multiple_of(prime))
    }

    /// Enumerate the finite prefix `M_S intersect [1,max]`.
    pub fn elements_through(&self, max: u64) -> Vec<u64> {
        (1..=max).filter(|&value| self.contains(value)).collect()
    }

    /// Evaluate a finite source `E` sum.
    ///
    /// The callback receives `m in M_S` and must return the already-evaluated
    /// value `f(m x)` for the caller's semilocal point `x`. `modulus` is
    /// `|x|_S`. `max_m` must be large enough that omitted callback values are
    /// known to vanish, e.g. from compact support.
    pub fn finite_e_sum(
        &self,
        modulus: f64,
        max_m: u64,
        mut f_at_mx: impl FnMut(u64) -> f64,
    ) -> Result<FiniteESum, SemilocalPoissonError> {
        if !modulus.is_finite() || modulus <= 0.0 {
            return Err(SemilocalPoissonError::InvalidModulus { modulus });
        }

        let mut raw_sum = 0.0;
        let mut term_count = 0_usize;
        for m in 1..=max_m {
            if self.contains(m) {
                raw_sum += f_at_mx(m);
                term_count += 1;
            }
        }

        Ok(FiniteESum {
            value: modulus.sqrt() * raw_sum,
            raw_sum,
            term_count,
            max_m,
        })
    }
}

/// Auditable result of a finite `E` sum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteESum {
    value: f64,
    raw_sum: f64,
    term_count: usize,
    max_m: u64,
}

impl FiniteESum {
    #[inline]
    pub fn value(self) -> f64 {
        self.value
    }

    #[inline]
    pub fn raw_sum(self) -> f64 {
        self.raw_sum
    }

    #[inline]
    pub fn term_count(self) -> usize {
        self.term_count
    }

    #[inline]
    pub fn max_m(self) -> u64 {
        self.max_m
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SemilocalPoissonError {
    InvalidModulus { modulus: f64 },
}

impl fmt::Display for SemilocalPoissonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidModulus { modulus } => {
                write!(
                    f,
                    "semilocal modulus must be finite and positive: {modulus}"
                )
            }
        }
    }
}

impl std::error::Error for SemilocalPoissonError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn monoid_excludes_finite_place_factors() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let monoid = SemilocalPoissonMonoid::new(places);
        assert_eq!(monoid.elements_through(12), vec![1, 5, 7, 11]);
    }

    #[test]
    fn archimedean_monoid_contains_all_positive_integers() {
        let places = FinitePlaceSet::new(vec![]).unwrap();
        let monoid = SemilocalPoissonMonoid::new(places);
        assert_eq!(monoid.elements_through(5), vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn finite_e_sum_applies_source_half_density_factor() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let monoid = SemilocalPoissonMonoid::new(places);
        let result = monoid.finite_e_sum(4.0, 10, |_| 1.0).unwrap();
        assert_eq!(result.term_count(), 3);
        assert_eq!(result.raw_sum(), 3.0);
        assert_eq!(result.value(), 6.0);
    }
}
