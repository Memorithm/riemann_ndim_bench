//! Manufactured covariance of the factorizable additive Poisson fixture under
//! the exact `Q_S^*` unit action.
//!
//! PR #45 supplied a certified factorizable additive Poisson fixture. PR #46
//! supplies the exact local transport `k_p -> k_p - n_p` for a unit
//! `u = +/- prod p^{n_p}`. This module combines the two facts and checks the
//! source statement that the additive `Q_S` summation is invariant under the
//! unit action without erasing the finite local coordinates.
//!
//! The exact proof-relevant identity remains the exponent compensation in
//! `semilocal_unit_orbit`; the floating-point lattice comparison here is a
//! manufactured regression of that identity, not a proof of the quotient
//! Poisson theorem on `L^2(X_S)`.

use std::fmt;

use crate::archimedean_poisson::{
    ArchimedeanPoissonError, CertifiedLatticeSum, certified_source_fixture_lattice_sum,
};
use crate::semilocal_factorizable_poisson::{
    FactorizablePoissonComparison, FactorizablePoissonError, LocalBallSpec,
    compare_factorizable_ball_poisson,
};
use crate::semilocal_qs::{QsArithmeticError, QsRational};
use crate::semilocal_trace_contract::FinitePlaceSet;
use crate::semilocal_unit_orbit::{SemilocalUnitOrbitTransport, UnitOrbitError};

/// Auditable manufactured comparison before and after one `Q_S^*` unit action.
#[derive(Clone, Debug, PartialEq)]
pub struct UnitCovarianceComparison {
    monoid_representative: u64,
    original: FactorizablePoissonComparison,
    transported_balls: Vec<LocalBallSpec>,
    archimedean_unit_scale: f64,
    transported_local_lattice_scale: f64,
    effective_transported_step: f64,
    rescaled_transported_sum: CertifiedLatticeSum,
}

impl UnitCovarianceComparison {
    /// Unique `m in M_S` in the source decomposition `q = u m`.
    #[inline]
    pub fn monoid_representative(&self) -> u64 {
        self.monoid_representative
    }

    #[inline]
    pub fn original(&self) -> &FactorizablePoissonComparison {
        &self.original
    }

    #[inline]
    pub fn transported_balls(&self) -> &[LocalBallSpec] {
        &self.transported_balls
    }

    #[inline]
    pub fn archimedean_unit_scale(&self) -> f64 {
        self.archimedean_unit_scale
    }

    #[inline]
    pub fn transported_local_lattice_scale(&self) -> f64 {
        self.transported_local_lattice_scale
    }

    /// Product of the archimedean unit scale with the transported finite local
    /// lattice scale. Exact exponent bookkeeping says this equals the original
    /// lattice scale; this `f64` value is only an audit representation.
    #[inline]
    pub fn effective_transported_step(&self) -> f64 {
        self.effective_transported_step
    }

    #[inline]
    pub fn rescaled_transported_sum(&self) -> CertifiedLatticeSum {
        self.rescaled_transported_sum
    }

    #[inline]
    pub fn residual(&self) -> f64 {
        self.original.original_sum().value() - self.rescaled_transported_sum.value()
    }

    /// Sum of the rigorous Gaussian truncation bounds on the independently
    /// evaluated original and unit-transported lattice sums.
    #[inline]
    pub fn combined_tail_bound(&self) -> f64 {
        self.original.original_sum().absolute_tail_bound()
            + self.rescaled_transported_sum.absolute_tail_bound()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum UnitCovarianceError {
    Qs(QsArithmeticError),
    UnitOrbit(UnitOrbitError),
    Factorizable(FactorizablePoissonError),
    Archimedean(ArchimedeanPoissonError),
    EffectiveStepOutOfRange,
}

impl fmt::Display for UnitCovarianceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Qs(error) => write!(f, "Q_S decomposition failed: {error}"),
            Self::UnitOrbit(error) => write!(f, "Q_S unit transport failed: {error}"),
            Self::Factorizable(error) => write!(f, "factorizable Poisson fixture failed: {error}"),
            Self::Archimedean(error) => write!(f, "archimedean lattice fixture failed: {error}"),
            Self::EffectiveStepOutOfRange => write!(
                f,
                "unit-transported effective archimedean lattice step is not finite positive"
            ),
        }
    }
}

impl std::error::Error for UnitCovarianceError {}

impl From<QsArithmeticError> for UnitCovarianceError {
    fn from(value: QsArithmeticError) -> Self {
        Self::Qs(value)
    }
}

impl From<UnitOrbitError> for UnitCovarianceError {
    fn from(value: UnitOrbitError) -> Self {
        Self::UnitOrbit(value)
    }
}

impl From<FactorizablePoissonError> for UnitCovarianceError {
    fn from(value: FactorizablePoissonError) -> Self {
        Self::Factorizable(value)
    }
}

impl From<ArchimedeanPoissonError> for UnitCovarianceError {
    fn from(value: ArchimedeanPoissonError) -> Self {
        Self::Archimedean(value)
    }
}

/// Compare the original factorizable additive lattice sum with the same fixture
/// after transporting all finite local balls by the unit part of `q = u m` and
/// compensating by the archimedean absolute value of `u`.
pub fn compare_unit_covariance(
    q: QsRational,
    places: &FinitePlaceSet,
    original_balls: &[LocalBallSpec],
    max_abs_n: u64,
) -> Result<UnitCovarianceComparison, UnitCovarianceError> {
    let decomposition = q.unit_monoid_decomposition(places)?;
    let action = SemilocalUnitOrbitTransport::from_decomposition(&decomposition, places)?;
    let transported_balls = action.transport_complete_product(original_balls)?;

    let original = compare_factorizable_ball_poisson(original_balls, max_abs_n)?;

    let transported_fixture =
        compare_factorizable_ball_poisson(&transported_balls, max_abs_n)?;
    let transported_local_lattice_scale = transported_fixture.diagonal_lattice_scale();
    let archimedean_unit_scale = action.archimedean_absolute_scale()?;
    let effective_transported_step = archimedean_unit_scale * transported_local_lattice_scale;
    if !effective_transported_step.is_finite() || effective_transported_step <= 0.0 {
        return Err(UnitCovarianceError::EffectiveStepOutOfRange);
    }

    let rescaled_transported_sum =
        certified_source_fixture_lattice_sum(effective_transported_step, max_abs_n)?;

    Ok(UnitCovarianceComparison {
        monoid_representative: decomposition.monoid_element(),
        original,
        transported_balls,
        archimedean_unit_scale,
        transported_local_lattice_scale,
        effective_transported_step,
        rescaled_transported_sum,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_covariant(comparison: &UnitCovarianceComparison, relative_roundoff: f64) {
        let scale = comparison
            .original()
            .original_sum()
            .value()
            .abs()
            .max(comparison.rescaled_transported_sum().value().abs())
            .max(1.0);
        let roundoff = relative_roundoff * scale;
        assert!(
            comparison.residual().abs() <= comparison.combined_tail_bound() + roundoff,
            "residual={:.3e} tails={:.3e} roundoff={:.3e}",
            comparison.residual(),
            comparison.combined_tail_bound(),
            roundoff
        );
    }

    #[test]
    fn dyadic_unit_covariance_is_exactly_representable() {
        let places = FinitePlaceSet::new(vec![2]).unwrap();
        let q = QsRational::new(3, 8, &places).unwrap();
        let comparison = compare_unit_covariance(
            q,
            &places,
            &[LocalBallSpec::new(2, 1)],
            256,
        )
        .unwrap();

        assert_eq!(comparison.monoid_representative(), 3);
        assert_eq!(comparison.transported_balls(), &[LocalBallSpec::new(2, 4)]);
        assert_eq!(comparison.archimedean_unit_scale(), 1.0 / 8.0);
        assert_eq!(comparison.transported_local_lattice_scale(), 16.0);
        assert_eq!(comparison.effective_transported_step(), 2.0);
        assert_eq!(comparison.original().diagonal_lattice_scale(), 2.0);
        assert_covariant(&comparison, 0.0);
    }

    #[test]
    fn multi_prime_unit_covariance_survives_float_audit() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let q = QsRational::new(45, 8, &places).unwrap();
        let comparison = compare_unit_covariance(
            q,
            &places,
            &[LocalBallSpec::new(2, 1), LocalBallSpec::new(3, -2)],
            256,
        )
        .unwrap();

        assert_eq!(comparison.monoid_representative(), 5);
        assert_eq!(
            comparison.transported_balls(),
            &[LocalBallSpec::new(2, 4), LocalBallSpec::new(3, -4)]
        );
        assert_covariant(&comparison, 5.0e-13);
    }
}
