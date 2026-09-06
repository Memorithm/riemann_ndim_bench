//! Quadrature-refinement diagnostics for the finite Gram-normalized Weil audit.
//!
//! The same support window and basis dimension are re-evaluated at a declared
//! sequence of quadrature orders. The resulting variation is reported as an
//! observed refinement spread and consecutive deltas. These quantities are
//! empirical resolution diagnostics, not certified numerical error bounds.

use std::fmt;

use crate::semilocal_compact_archimedean::{CompactArchimedeanBump, CompactArchimedeanError};
use crate::weil_generalized_spectrum::{
    FiniteWeilGeneralizedSpectrumError, audit_finite_weil_generalized_spectrum,
};
use crate::weil_support_sweep::WeilSupportWindow;

/// One declared set of quadrature orders.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WeilQuadratureLevel {
    correlation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
    gram_order: usize,
}

impl WeilQuadratureLevel {
    #[inline]
    pub const fn new(
        correlation_order: usize,
        archimedean_order: usize,
        boundary_order: usize,
        gram_order: usize,
    ) -> Self {
        Self {
            correlation_order,
            archimedean_order,
            boundary_order,
            gram_order,
        }
    }

    #[inline]
    pub const fn correlation_order(self) -> usize {
        self.correlation_order
    }

    #[inline]
    pub const fn archimedean_order(self) -> usize {
        self.archimedean_order
    }

    #[inline]
    pub const fn boundary_order(self) -> usize {
        self.boundary_order
    }

    #[inline]
    pub const fn gram_order(self) -> usize {
        self.gram_order
    }
}

/// One finite spectral evaluation at a declared quadrature level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WeilQuadratureSample {
    level: WeilQuadratureLevel,
    raw_minimum_eigenvalue: f64,
    generalized_minimum_eigenvalue: f64,
    gram_condition_number: f64,
    max_boundary_residual: f64,
    max_pairing_asymmetry: f64,
    max_whitened_asymmetry: f64,
}

impl WeilQuadratureSample {
    #[inline]
    pub const fn level(self) -> WeilQuadratureLevel {
        self.level
    }

    #[inline]
    pub const fn raw_minimum_eigenvalue(self) -> f64 {
        self.raw_minimum_eigenvalue
    }

    #[inline]
    pub const fn generalized_minimum_eigenvalue(self) -> f64 {
        self.generalized_minimum_eigenvalue
    }

    #[inline]
    pub const fn gram_condition_number(self) -> f64 {
        self.gram_condition_number
    }

    #[inline]
    pub const fn max_boundary_residual(self) -> f64 {
        self.max_boundary_residual
    }

    #[inline]
    pub const fn max_pairing_asymmetry(self) -> f64 {
        self.max_pairing_asymmetry
    }

    #[inline]
    pub const fn max_whitened_asymmetry(self) -> f64 {
        self.max_whitened_asymmetry
    }
}

/// Refinement audit for one exact support window and one fixed basis dimension.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilQuadratureRefinementAudit {
    window: WeilSupportWindow,
    dimension: usize,
    samples: Vec<WeilQuadratureSample>,
    raw_observed_minimum: f64,
    raw_observed_maximum: f64,
    generalized_observed_minimum: f64,
    generalized_observed_maximum: f64,
}

impl FiniteWeilQuadratureRefinementAudit {
    #[inline]
    pub const fn window(&self) -> WeilSupportWindow {
        self.window
    }

    #[inline]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    #[inline]
    pub fn samples(&self) -> &[WeilQuadratureSample] {
        &self.samples
    }

    #[inline]
    pub const fn raw_observed_interval(&self) -> (f64, f64) {
        (self.raw_observed_minimum, self.raw_observed_maximum)
    }

    #[inline]
    pub const fn generalized_observed_interval(&self) -> (f64, f64) {
        (
            self.generalized_observed_minimum,
            self.generalized_observed_maximum,
        )
    }

    #[inline]
    pub const fn raw_observed_span(&self) -> f64 {
        self.raw_observed_maximum - self.raw_observed_minimum
    }

    #[inline]
    pub const fn generalized_observed_span(&self) -> f64 {
        self.generalized_observed_maximum - self.generalized_observed_minimum
    }

    pub fn last_raw_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some((last.raw_minimum_eigenvalue - previous.raw_minimum_eigenvalue).abs())
    }

    pub fn last_generalized_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some(
            (last.generalized_minimum_eigenvalue - previous.generalized_minimum_eigenvalue).abs(),
        )
    }
}

#[derive(Debug)]
pub enum FiniteWeilRefinementError {
    EmptyLevelSet,
    CompactSupport(CompactArchimedeanError),
    Generalized(FiniteWeilGeneralizedSpectrumError),
}

impl fmt::Display for FiniteWeilRefinementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyLevelSet => write!(f, "finite Weil refinement audit requires at least one level"),
            Self::CompactSupport(error) => write!(f, "invalid compact support window: {error}"),
            Self::Generalized(error) => write!(f, "generalized finite Weil audit failed: {error}"),
        }
    }
}

impl std::error::Error for FiniteWeilRefinementError {}

impl From<CompactArchimedeanError> for FiniteWeilRefinementError {
    fn from(value: CompactArchimedeanError) -> Self {
        Self::CompactSupport(value)
    }
}

impl From<FiniteWeilGeneralizedSpectrumError> for FiniteWeilRefinementError {
    fn from(value: FiniteWeilGeneralizedSpectrumError) -> Self {
        Self::Generalized(value)
    }
}

/// Re-evaluate one support/dimension point at all declared quadrature levels.
pub fn audit_finite_weil_quadrature_refinement(
    window: WeilSupportWindow,
    dimension: usize,
    levels: &[WeilQuadratureLevel],
) -> Result<FiniteWeilQuadratureRefinementAudit, FiniteWeilRefinementError> {
    if levels.is_empty() {
        return Err(FiniteWeilRefinementError::EmptyLevelSet);
    }

    let bump = CompactArchimedeanBump::new(window.lower(), window.upper())?;
    let mut samples = Vec::with_capacity(levels.len());
    for &level in levels {
        let audit = audit_finite_weil_generalized_spectrum(
            bump,
            dimension,
            level.correlation_order(),
            level.archimedean_order(),
            level.boundary_order(),
            level.gram_order(),
        )?;
        samples.push(WeilQuadratureSample {
            level,
            raw_minimum_eigenvalue: audit.pairing().minimum_eigenvalue(),
            generalized_minimum_eigenvalue: audit.minimum_generalized_eigenvalue(),
            gram_condition_number: audit.gram_condition_number(),
            max_boundary_residual: audit.pairing().max_boundary_residual(),
            max_pairing_asymmetry: audit.pairing().max_raw_pairing_asymmetry(),
            max_whitened_asymmetry: audit.max_whitened_asymmetry(),
        });
    }

    let raw_observed_minimum = samples
        .iter()
        .map(|sample| sample.raw_minimum_eigenvalue)
        .min_by(f64::total_cmp)
        .expect("non-empty refinement sample set");
    let raw_observed_maximum = samples
        .iter()
        .map(|sample| sample.raw_minimum_eigenvalue)
        .max_by(f64::total_cmp)
        .expect("non-empty refinement sample set");
    let generalized_observed_minimum = samples
        .iter()
        .map(|sample| sample.generalized_minimum_eigenvalue)
        .min_by(f64::total_cmp)
        .expect("non-empty refinement sample set");
    let generalized_observed_maximum = samples
        .iter()
        .map(|sample| sample.generalized_minimum_eigenvalue)
        .max_by(f64::total_cmp)
        .expect("non-empty refinement sample set");

    Ok(FiniteWeilQuadratureRefinementAudit {
        window,
        dimension,
        samples,
        raw_observed_minimum,
        raw_observed_maximum,
        generalized_observed_minimum,
        generalized_observed_maximum,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::PositiveRational;

    fn rational(numerator: u64, denominator: u64) -> PositiveRational {
        PositiveRational::new(numerator, denominator).unwrap()
    }

    fn window() -> WeilSupportWindow {
        WeilSupportWindow::new(rational(1, 2), rational(7, 2))
    }

    #[test]
    fn empty_refinement_level_set_is_rejected() {
        let error = audit_finite_weil_quadrature_refinement(window(), 1, &[]).unwrap_err();
        assert!(matches!(error, FiniteWeilRefinementError::EmptyLevelSet));
    }

    #[test]
    fn refinement_reports_observed_spreads_without_sign_assumptions() {
        let levels = [
            WeilQuadratureLevel::new(24, 24, 32, 32),
            WeilQuadratureLevel::new(32, 32, 48, 48),
        ];
        let audit = audit_finite_weil_quadrature_refinement(window(), 1, &levels).unwrap();

        assert_eq!(audit.dimension(), 1);
        assert_eq!(audit.samples().len(), 2);
        assert_eq!(audit.samples()[0].level(), levels[0]);
        assert_eq!(audit.samples()[1].level(), levels[1]);
        assert!(audit.raw_observed_span().is_finite());
        assert!(audit.raw_observed_span() >= 0.0);
        assert!(audit.generalized_observed_span().is_finite());
        assert!(audit.generalized_observed_span() >= 0.0);
        assert!(audit.last_raw_delta().unwrap().is_finite());
        assert!(audit.last_generalized_delta().unwrap().is_finite());
        for sample in audit.samples() {
            assert!(sample.raw_minimum_eigenvalue().is_finite());
            assert!(sample.generalized_minimum_eigenvalue().is_finite());
            assert!(sample.gram_condition_number().is_finite());
            assert!(sample.max_boundary_residual().is_finite());
            assert!(sample.max_pairing_asymmetry().is_finite());
            assert!(sample.max_whitened_asymmetry().is_finite());
        }
    }
}
