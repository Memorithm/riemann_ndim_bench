//! Certified factorizable semilocal Poisson fixture.
//!
//! For a finite place set `S_f` and elementary local factors
//!
//! `prod_p 1_{p^{k_p} Z_p}`,
//!
//! the diagonal lattice `Q_S = Z[S_f^{-1}]` surviving all finite-place
//! indicators is exactly
//!
//! `(prod_p p^{k_p}) Z`.
//!
//! Under the source-locked self-dual local Fourier normalization,
//!
//! `F_p[1_{p^k Z_p}] = p^{-k} 1_{p^{-k} Z_p}`,
//!
//! so the transformed finite factor contributes the reciprocal scale and the
//! reciprocal diagonal lattice. Combining those exact local identities with the
//! analytic archimedean Gaussian Fourier pair yields a certified manufactured
//! instance of additive semilocal Poisson summation.
//!
//! This module is deliberately much narrower than a representation of `A_S`,
//! `X_S`, or a general Bruhat--Schwartz function. It does not implement the
//! quotient Poisson map `E`, Conjecture 4.1, Weil positivity, or RH.

use std::fmt;

use crate::archimedean_poisson::{
    ArchimedeanPoissonError, CertifiedLatticeSum, certified_source_fixture_fourier_lattice_sum,
    certified_source_fixture_lattice_sum,
};
use crate::semilocal_padic_fourier::{PadicBall, PadicFourierError};
use crate::semilocal_trace_contract::{FinitePlaceSet, SemilocalTraceContractError};

/// One elementary local factor `1_{p^k Z_p}`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LocalBallSpec {
    prime: u64,
    exponent: i32,
}

impl LocalBallSpec {
    #[inline]
    pub const fn new(prime: u64, exponent: i32) -> Self {
        Self { prime, exponent }
    }

    #[inline]
    pub const fn prime(self) -> u64 {
        self.prime
    }

    #[inline]
    pub const fn exponent(self) -> i32 {
        self.exponent
    }
}

/// Certified comparison of a manufactured factorizable semilocal Poisson pair.
#[derive(Clone, Debug, PartialEq)]
pub struct FactorizablePoissonComparison {
    local_balls: Vec<LocalBallSpec>,
    diagonal_lattice_scale: f64,
    fourier_local_scale: f64,
    original_sum: CertifiedLatticeSum,
    transformed_archimedean_sum: CertifiedLatticeSum,
}

impl FactorizablePoissonComparison {
    #[inline]
    pub fn local_balls(&self) -> &[LocalBallSpec] {
        &self.local_balls
    }

    /// Positive real number `A = prod_p p^{k_p}` such that the surviving
    /// diagonal lattice is `A Z`.
    #[inline]
    pub fn diagonal_lattice_scale(&self) -> f64 {
        self.diagonal_lattice_scale
    }

    /// Product of exact local Fourier scales, numerically evaluated as `A^-1`.
    #[inline]
    pub fn fourier_local_scale(&self) -> f64 {
        self.fourier_local_scale
    }

    #[inline]
    pub fn original_sum(&self) -> CertifiedLatticeSum {
        self.original_sum
    }

    #[inline]
    pub fn transformed_archimedean_sum(&self) -> CertifiedLatticeSum {
        self.transformed_archimedean_sum
    }

    /// Left side of the manufactured semilocal Poisson identity.
    #[inline]
    pub fn left_value(&self) -> f64 {
        self.original_sum.value()
    }

    /// Right side, including the product of finite-place Fourier scales.
    #[inline]
    pub fn right_value(&self) -> f64 {
        self.fourier_local_scale * self.transformed_archimedean_sum.value()
    }

    #[inline]
    pub fn residual(&self) -> f64 {
        self.left_value() - self.right_value()
    }

    /// Rigorous bound for the two omitted Gaussian lattice tails. Floating
    /// point roundoff is intentionally not included in this quantity.
    #[inline]
    pub fn combined_tail_bound(&self) -> f64 {
        self.original_sum.absolute_tail_bound()
            + self.fourier_local_scale
                * self.transformed_archimedean_sum.absolute_tail_bound()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FactorizablePoissonError {
    DuplicatePrime { prime: u64 },
    ScaleOutOfRange,
    TraceContract(SemilocalTraceContractError),
    PadicFourier(PadicFourierError),
    Archimedean(ArchimedeanPoissonError),
}

impl fmt::Display for FactorizablePoissonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicatePrime { prime } => {
                write!(f, "factorizable local ball list repeats finite place p={prime}")
            }
            Self::ScaleOutOfRange => write!(
                f,
                "factorizable diagonal lattice scale is outside finite positive f64 range"
            ),
            Self::TraceContract(error) => write!(f, "invalid semilocal place set: {error}"),
            Self::PadicFourier(error) => write!(f, "p-adic Fourier fixture failed: {error}"),
            Self::Archimedean(error) => write!(f, "archimedean Poisson fixture failed: {error}"),
        }
    }
}

impl std::error::Error for FactorizablePoissonError {}

impl From<SemilocalTraceContractError> for FactorizablePoissonError {
    fn from(value: SemilocalTraceContractError) -> Self {
        Self::TraceContract(value)
    }
}

impl From<PadicFourierError> for FactorizablePoissonError {
    fn from(value: PadicFourierError) -> Self {
        Self::PadicFourier(value)
    }
}

impl From<ArchimedeanPoissonError> for FactorizablePoissonError {
    fn from(value: ArchimedeanPoissonError) -> Self {
        Self::Archimedean(value)
    }
}

/// Compare the two sides of additive semilocal Poisson summation for the
/// manufactured factorizable class
///
/// `f_infinity tensor prod_p 1_{p^{k_p} Z_p}`.
///
/// `max_abs_n` controls the certified bilateral archimedean lattice prefixes.
pub fn compare_factorizable_ball_poisson(
    local_balls: &[LocalBallSpec],
    max_abs_n: u64,
) -> Result<FactorizablePoissonComparison, FactorizablePoissonError> {
    let mut sorted = local_balls.to_vec();
    sorted.sort_unstable_by_key(|spec| spec.prime);
    for pair in sorted.windows(2) {
        if pair[0].prime == pair[1].prime {
            return Err(FactorizablePoissonError::DuplicatePrime {
                prime: pair[0].prime,
            });
        }
    }

    let places = FinitePlaceSet::new(sorted.iter().map(|spec| spec.prime).collect())?;

    let mut log_scale = 0.0_f64;
    for spec in &sorted {
        let ball = PadicBall::new(spec.prime, spec.exponent, &places)?;
        let transformed = ball.fourier_transform();

        debug_assert_eq!(transformed.ball().prime(), spec.prime);
        debug_assert_eq!(transformed.ball().exponent(), -i64::from(spec.exponent));
        debug_assert_eq!(transformed.scale().prime(), spec.prime);
        debug_assert_eq!(transformed.scale().exponent(), -i64::from(spec.exponent));

        log_scale += f64::from(spec.exponent) * (spec.prime as f64).ln();
    }

    let diagonal_lattice_scale = log_scale.exp();
    let fourier_local_scale = (-log_scale).exp();
    if !diagonal_lattice_scale.is_finite()
        || diagonal_lattice_scale <= 0.0
        || !fourier_local_scale.is_finite()
        || fourier_local_scale <= 0.0
    {
        return Err(FactorizablePoissonError::ScaleOutOfRange);
    }

    let original_sum =
        certified_source_fixture_lattice_sum(diagonal_lattice_scale, max_abs_n)?;
    let transformed_archimedean_sum =
        certified_source_fixture_fourier_lattice_sum(fourier_local_scale, max_abs_n)?;

    Ok(FactorizablePoissonComparison {
        local_balls: sorted,
        diagonal_lattice_scale,
        fourier_local_scale,
        original_sum,
        transformed_archimedean_sum,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_certified(comparison: &FactorizablePoissonComparison) {
        let scale = comparison
            .left_value()
            .abs()
            .max(comparison.right_value().abs())
            .max(1.0);
        let roundoff = 5.0e-13 * scale;
        assert!(
            comparison.residual().abs() <= comparison.combined_tail_bound() + roundoff,
            "residual={:.3e} tail={:.3e} roundoff={:.3e}",
            comparison.residual(),
            comparison.combined_tail_bound(),
            roundoff
        );
    }

    #[test]
    fn self_dual_unit_ball_recovers_a_certified_product_poisson_fixture() {
        let comparison =
            compare_factorizable_ball_poisson(&[LocalBallSpec::new(2, 0)], 128).unwrap();
        assert_eq!(comparison.diagonal_lattice_scale(), 1.0);
        assert_eq!(comparison.fourier_local_scale(), 1.0);
        assert_certified(&comparison);
    }

    #[test]
    fn nontrivial_local_ball_and_dual_scale_are_both_exercised() {
        for exponent in [-2, -1, 1, 2] {
            let comparison =
                compare_factorizable_ball_poisson(&[LocalBallSpec::new(3, exponent)], 256)
                    .unwrap();
            let expected = 3.0_f64.powi(exponent);
            assert!((comparison.diagonal_lattice_scale() - expected).abs() < 2.0e-15 * expected);
            assert!(
                (comparison.fourier_local_scale() - expected.recip()).abs()
                    < 2.0e-15 * expected.recip()
            );
            assert_certified(&comparison);
        }
    }

    #[test]
    fn multiple_finite_places_compose_multiplicatively() {
        let comparison = compare_factorizable_ball_poisson(
            &[LocalBallSpec::new(2, 1), LocalBallSpec::new(3, -1)],
            256,
        )
        .unwrap();
        assert!((comparison.diagonal_lattice_scale() - 2.0 / 3.0).abs() < 2.0e-15);
        assert!((comparison.fourier_local_scale() - 1.5).abs() < 2.0e-15);
        assert_certified(&comparison);
    }

    #[test]
    fn duplicate_and_composite_places_are_rejected() {
        assert_eq!(
            compare_factorizable_ball_poisson(
                &[LocalBallSpec::new(2, 0), LocalBallSpec::new(2, 1)],
                32,
            )
            .unwrap_err(),
            FactorizablePoissonError::DuplicatePrime { prime: 2 }
        );
        assert!(matches!(
            compare_factorizable_ball_poisson(&[LocalBallSpec::new(9, 0)], 32),
            Err(FactorizablePoissonError::TraceContract(
                SemilocalTraceContractError::InvalidPrimePlace { value: 9 }
            ))
        ));
    }
}
