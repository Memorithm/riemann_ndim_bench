//! Compact-support bookkeeping for the Riemann--Weil explicit formula.
//!
//! Fact 3.2 of Connes--Consani, *The Scaling Hamiltonian*, states that a
//! compactly supported Riemann--Weil test function involves only finitely many
//! places because the contribution at a finite prime samples non-zero powers of
//! that prime.
//!
//! This module makes that finite-support reduction executable. It deliberately
//! does **not** assert Conjecture 4.1, which is the stronger claim that the
//! semilocal framework for those places suffices to prove Weil positivity on a
//! whole support window.

use std::fmt;

use crate::semilocal_trace_contract::FinitePlaceSet;
use crate::weil_boundary::{MultiplicativeSupport, WeilBoundaryError};

/// Open support window from Conjecture 4.1.
///
/// The window for `h_1` is `(q^(-1/2), q^(1/2))`. The multiplicative
/// convolution square then has open support envelope `(q^(-1), q)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WeilSupportWindow {
    q: f64,
}

impl WeilSupportWindow {
    pub fn new(q: f64) -> Result<Self, WeilSupportError> {
        if !q.is_finite() || q <= 1.0 {
            return Err(WeilSupportError::InvalidWindowParameter { q });
        }
        Ok(Self { q })
    }

    #[inline]
    pub fn q(self) -> f64 {
        self.q
    }

    #[inline]
    pub fn h1_lower(self) -> f64 {
        self.q.sqrt().recip()
    }

    #[inline]
    pub fn h1_upper(self) -> f64 {
        self.q.sqrt()
    }

    #[inline]
    pub fn convolution_lower(self) -> f64 {
        self.q.recip()
    }

    #[inline]
    pub fn convolution_upper(self) -> f64 {
        self.q
    }

    /// Test whether an actual compact support is strictly inside the source
    /// `h_1` window.
    pub fn contains_h1_support(self, support: MultiplicativeSupport) -> bool {
        support.lower() > self.h1_lower() && support.upper() < self.h1_upper()
    }

    /// Test whether an actual compact convolution support is strictly inside
    /// `(q^-1, q)`.
    pub fn contains_convolution_support(self, support: MultiplicativeSupport) -> bool {
        support.lower() > self.convolution_lower() && support.upper() < self.convolution_upper()
    }

    /// Membership in the finite-place set appearing in Conjecture 4.1:
    /// `S(q) = {infinity} union {p : p < q}`.
    ///
    /// This method is only set bookkeeping; it does not assert the conjectured
    /// positivity sufficiency of `S(q)`.
    #[inline]
    pub fn source_set_contains_prime(self, prime: u64) -> bool {
        (prime as f64) < self.q
    }
}

/// Support envelope of the multiplicative inverse `f*(rho)=f(rho^-1)`.
pub fn inverse_support(
    support: MultiplicativeSupport,
) -> Result<MultiplicativeSupport, WeilBoundaryError> {
    MultiplicativeSupport::new(support.upper().recip(), support.lower().recip())
}

/// Conservative support envelope for multiplicative convolution.
///
/// If `supp(f) subset [a,b]` and `supp(g) subset [c,d]`, then
/// `supp(f*g) subset [ac,bd]`.
pub fn convolution_support_envelope(
    left: MultiplicativeSupport,
    right: MultiplicativeSupport,
) -> Result<MultiplicativeSupport, WeilBoundaryError> {
    MultiplicativeSupport::new(left.lower() * right.lower(), left.upper() * right.upper())
}

/// Support envelope for `h_1 * h_1^*`, where
/// `h_1^*(rho)=overline(h_1(rho^-1))`.
pub fn convolution_square_support_envelope(
    h1_support: MultiplicativeSupport,
) -> Result<MultiplicativeSupport, WeilBoundaryError> {
    convolution_support_envelope(h1_support, inverse_support(h1_support)?)
}

/// A finite prime is certainly absent from the explicit formula when none of
/// its non-zero integral powers can meet the declared support.
///
/// The sufficient test below is exact for the exclusion direction: if
/// `p > upper` and `p^-1 < lower`, then every positive power lies above the
/// support and every negative power lies below it.
pub fn prime_is_excluded_by_support(prime: u64, support: MultiplicativeSupport) -> bool {
    if prime < 2 {
        return false;
    }
    let p = prime as f64;
    p > support.upper() && p.recip() < support.lower()
}

/// Build the exact finite source set `{infinity} union {p : p < q}` for an
/// integer bound `q`.
///
/// This helper exists for deterministic finite regressions. It does not encode
/// Conjecture 4.1's sufficiency claim.
pub fn source_place_set_below_integer_bound(q: u64) -> Result<FinitePlaceSet, WeilSupportError> {
    if q < 2 {
        return Err(WeilSupportError::InvalidIntegerPlaceBound { q });
    }
    let primes = (2..q).filter(|&value| is_prime(value)).collect();
    Ok(FinitePlaceSet::new(primes).expect("generated values are prime"))
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WeilSupportError {
    InvalidWindowParameter { q: f64 },
    InvalidIntegerPlaceBound { q: u64 },
}

impl fmt::Display for WeilSupportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidWindowParameter { q } => {
                write!(f, "support-window parameter q must be finite and > 1: {q}")
            }
            Self::InvalidIntegerPlaceBound { q } => {
                write!(f, "integer finite-place bound must be >= 2: {q}")
            }
        }
    }
}

impl std::error::Error for WeilSupportError {}

fn is_prime(value: u64) -> bool {
    if value < 2 {
        return false;
    }
    if value == 2 {
        return true;
    }
    if value.is_multiple_of(2) {
        return false;
    }
    let mut divisor = 3_u64;
    while divisor <= value / divisor {
        if value.is_multiple_of(divisor) {
            return false;
        }
        divisor += 2;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_window_scales_exactly_under_convolution_envelope() {
        let window = WeilSupportWindow::new(10.0).unwrap();
        assert!((window.h1_lower() - 10.0_f64.sqrt().recip()).abs() < f64::EPSILON);
        assert!((window.h1_upper() - 10.0_f64.sqrt()).abs() < f64::EPSILON);
        assert_eq!(window.convolution_lower(), 0.1);
        assert_eq!(window.convolution_upper(), 10.0);
    }

    #[test]
    fn integer_source_place_set_uses_strict_prime_bound() {
        let places = source_place_set_below_integer_bound(10).unwrap();
        assert_eq!(places.finite_primes(), &[2, 3, 5, 7]);
        let places = source_place_set_below_integer_bound(2).unwrap();
        assert!(places.finite_primes().is_empty());
    }
}
