//! Source-level data contract for the semilocal trace formula.
//!
//! This module mirrors the objects and normalizations in Connes--Consani,
//! *The Scaling Hamiltonian*, Section 2. It deliberately does not construct
//! `L^2(X_S)`, the Fourier transform, or the cutoff operators numerically.
//!
//! The purpose is to prevent later Riemann-specific code from silently changing
//! the finite place set, the self-dual additive-character normalization, the
//! cutoff convention, or the leading term in Theorem 2.5.

use std::fmt;

/// Finite set of non-archimedean places, with the archimedean place implicit.
///
/// The source always assumes `infinity in S`; this type enforces that convention
/// by construction and stores only the finite prime places.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FinitePlaceSet {
    finite_primes: Vec<u64>,
}

impl FinitePlaceSet {
    /// Construct `S = {infinity} union finite_primes`.
    ///
    /// Inputs are sorted and deduplicated. Every finite place must be prime.
    pub fn new(mut finite_primes: Vec<u64>) -> Result<Self, SemilocalTraceContractError> {
        if let Some(&value) = finite_primes.iter().find(|&&value| !is_prime(value)) {
            return Err(SemilocalTraceContractError::InvalidPrimePlace { value });
        }
        finite_primes.sort_unstable();
        finite_primes.dedup();
        Ok(Self { finite_primes })
    }

    /// Finite prime places in increasing order.
    #[inline]
    pub fn finite_primes(&self) -> &[u64] {
        &self.finite_primes
    }

    /// Number of places, including the mandatory archimedean place.
    #[inline]
    pub fn place_count(&self) -> usize {
        self.finite_primes.len() + 1
    }

    /// Whether the finite place `p` belongs to `S`.
    #[inline]
    pub fn contains_prime(&self, p: u64) -> bool {
        self.finite_primes.binary_search(&p).is_ok()
    }

    /// Prime generators occurring in
    /// `Q_S^* = {+/- product p_j^n_j : n_j in Z}`.
    #[inline]
    pub fn qs_unit_prime_generators(&self) -> &[u64] {
        &self.finite_primes
    }
}

/// Additive-character normalization required by the semilocal source.
///
/// For each place, the additive Haar measure is normalized self-dually and the
/// product character is chosen so that `Q_S` is a self-dual lattice. Different
/// admissible choices differ by `Q_S^*`-scaling and therefore disappear on
/// `X_S = A_S / Q_S^*`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BasicCharacterNormalization {
    QsSelfDual,
}

/// Source contract for the semilocal Hilbert-space geometry.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SemilocalSpaceContract {
    places: FinitePlaceSet,
    character_normalization: BasicCharacterNormalization,
}

impl SemilocalSpaceContract {
    /// Construct the only normalization currently accepted by RiemannBench.
    pub fn qs_self_dual(places: FinitePlaceSet) -> Self {
        Self {
            places,
            character_normalization: BasicCharacterNormalization::QsSelfDual,
        }
    }

    #[inline]
    pub fn places(&self) -> &FinitePlaceSet {
        &self.places
    }

    #[inline]
    pub fn character_normalization(&self) -> BasicCharacterNormalization {
        self.character_normalization
    }
}

/// Shared infrared/ultraviolet cutoff parameter from equations (2.24)-(2.25).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SemilocalCutoff {
    lambda: f64,
}

impl SemilocalCutoff {
    pub fn new(lambda: f64) -> Result<Self, SemilocalTraceContractError> {
        if !lambda.is_finite() || lambda <= 0.0 {
            return Err(SemilocalTraceContractError::InvalidCutoff { lambda });
        }
        Ok(Self { lambda })
    }

    #[inline]
    pub fn lambda(self) -> f64 {
        self.lambda
    }

    #[inline]
    pub fn log_lambda(self) -> f64 {
        self.lambda.ln()
    }

    /// Upper endpoint `2 log Lambda / (2 pi) = log Lambda / pi` in Lemma 2.4.
    #[inline]
    pub fn quantized_band_endpoint(self) -> f64 {
        self.lambda.ln() / std::f64::consts::PI
    }

    /// Leading divergent term `2 f(1) log Lambda` of Theorem 2.5.
    #[inline]
    pub fn theorem_2_5_leading_term(self, f_at_one: f64) -> f64 {
        2.0 * f_at_one * self.lambda.ln()
    }
}

/// Complete source metadata needed before a numerical representation of
/// Theorem 2.5 can be attempted.
#[derive(Clone, Debug, PartialEq)]
pub struct SemilocalTraceContract {
    space: SemilocalSpaceContract,
    cutoff: SemilocalCutoff,
}

impl SemilocalTraceContract {
    pub fn new(space: SemilocalSpaceContract, cutoff: SemilocalCutoff) -> Self {
        Self { space, cutoff }
    }

    #[inline]
    pub fn space(&self) -> &SemilocalSpaceContract {
        &self.space
    }

    #[inline]
    pub fn cutoff(&self) -> SemilocalCutoff {
        self.cutoff
    }
}

/// Convert the source test-function convention
/// `h(lambda) = |lambda|^(1/2) f(lambda)` from equations (2.29)-(2.30).
pub fn symmetric_test_value(
    modulus: f64,
    f_value: f64,
) -> Result<f64, SemilocalTraceContractError> {
    if !modulus.is_finite() || modulus <= 0.0 {
        return Err(SemilocalTraceContractError::InvalidModulus { modulus });
    }
    Ok(modulus.sqrt() * f_value)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SemilocalTraceContractError {
    InvalidPrimePlace { value: u64 },
    InvalidCutoff { lambda: f64 },
    InvalidModulus { modulus: f64 },
}

impl fmt::Display for SemilocalTraceContractError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPrimePlace { value } => {
                write!(f, "finite semilocal place must be prime: {value}")
            }
            Self::InvalidCutoff { lambda } => {
                write!(
                    f,
                    "semilocal cutoff Lambda must be finite and positive: {lambda}"
                )
            }
            Self::InvalidModulus { modulus } => {
                write!(
                    f,
                    "idele-class modulus must be finite and positive: {modulus}"
                )
            }
        }
    }
}

impl std::error::Error for SemilocalTraceContractError {}

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
    fn place_set_is_sorted_deduplicated_and_archimedean_by_contract() {
        let places = FinitePlaceSet::new(vec![5, 2, 3, 2]).unwrap();
        assert_eq!(places.finite_primes(), &[2, 3, 5]);
        assert_eq!(places.place_count(), 4);
        assert!(places.contains_prime(3));
        assert!(!places.contains_prime(7));
        assert_eq!(places.qs_unit_prime_generators(), &[2, 3, 5]);
    }

    #[test]
    fn composite_place_is_rejected() {
        assert_eq!(
            FinitePlaceSet::new(vec![2, 9]).unwrap_err(),
            SemilocalTraceContractError::InvalidPrimePlace { value: 9 }
        );
    }

    #[test]
    fn cutoff_matches_source_band_and_leading_trace_normalization() {
        let lambda = (3.0 * std::f64::consts::PI).exp();
        let cutoff = SemilocalCutoff::new(lambda).unwrap();
        assert!((cutoff.quantized_band_endpoint() - 3.0).abs() < 4.0e-15);

        let lambda = 2.0_f64.exp();
        let cutoff = SemilocalCutoff::new(lambda).unwrap();
        assert!((cutoff.theorem_2_5_leading_term(1.5) - 6.0).abs() < 2.0e-15);
    }

    #[test]
    fn only_qs_self_dual_character_convention_is_constructed() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let space = SemilocalSpaceContract::qs_self_dual(places);
        assert_eq!(
            space.character_normalization(),
            BasicCharacterNormalization::QsSelfDual
        );
    }

    #[test]
    fn symmetric_test_normalization_matches_equation_2_30() {
        assert_eq!(symmetric_test_value(4.0, 3.0).unwrap(), 6.0);
        assert!(symmetric_test_value(0.0, 3.0).is_err());
    }
}
