//! Support-by-refinement evidence grid for one selected Bernstein subspace.
//!
//! The localized Bernstein subspace from `weil_bernstein_subspace` is evaluated
//! over exact rational support windows and declared quadrature levels.  Every
//! sample is paired with the dimension-matched leading-Legendre control extracted
//! from the same parent `(A,G)` computation.
//!
//! The resulting spreads and deltas are empirical sensitivity diagnostics.  They
//! are not certified error bounds, significance scores, or evidence of
//! complete-space Weil positivity.  Finite agreement across support, resolution,
//! and this one subspace family does not establish density/completeness,
//! semilocal `L^2(X_S)`, Conjecture 4.1, or RH.

use std::fmt;

use crate::semilocal_compact_archimedean::{CompactArchimedeanBump, CompactArchimedeanError};
use crate::weil_bernstein_subspace::{
    BernsteinIndexSubspace, FiniteWeilBernsteinSubspaceError, audit_finite_weil_bernstein_subspace,
};
use crate::weil_refinement::WeilQuadratureLevel;
use crate::weil_support_sweep::WeilSupportWindow;

/// One localized-subspace evaluation at one declared quadrature level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WeilBernsteinEvidenceSample {
    level: WeilQuadratureLevel,
    bernstein_raw_minimum_eigenvalue: f64,
    bernstein_generalized_minimum_eigenvalue: f64,
    leading_legendre_generalized_minimum_eigenvalue: f64,
    generalized_family_delta: f64,
    bernstein_gram_condition_number: f64,
    leading_legendre_gram_condition_number: f64,
    max_boundary_residual: f64,
    max_pairing_asymmetry: f64,
    max_whitened_asymmetry: f64,
}

impl WeilBernsteinEvidenceSample {
    #[inline]
    pub const fn level(self) -> WeilQuadratureLevel {
        self.level
    }

    #[inline]
    pub const fn bernstein_raw_minimum_eigenvalue(self) -> f64 {
        self.bernstein_raw_minimum_eigenvalue
    }

    #[inline]
    pub const fn bernstein_generalized_minimum_eigenvalue(self) -> f64 {
        self.bernstein_generalized_minimum_eigenvalue
    }

    #[inline]
    pub const fn leading_legendre_generalized_minimum_eigenvalue(self) -> f64 {
        self.leading_legendre_generalized_minimum_eigenvalue
    }

    /// Signed finite-dimensional difference
    /// `lambda_min(Bernstein) - lambda_min(leading Legendre)`.
    #[inline]
    pub const fn generalized_family_delta(self) -> f64 {
        self.generalized_family_delta
    }

    #[inline]
    pub const fn bernstein_gram_condition_number(self) -> f64 {
        self.bernstein_gram_condition_number
    }

    #[inline]
    pub const fn leading_legendre_gram_condition_number(self) -> f64 {
        self.leading_legendre_gram_condition_number
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

/// Refinement evidence for one exact support window.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilBernsteinEvidenceCell {
    window: WeilSupportWindow,
    samples: Vec<WeilBernsteinEvidenceSample>,
    bernstein_observed_minimum: f64,
    bernstein_observed_maximum: f64,
    leading_observed_minimum: f64,
    leading_observed_maximum: f64,
    family_delta_observed_minimum: f64,
    family_delta_observed_maximum: f64,
}

impl FiniteWeilBernsteinEvidenceCell {
    #[inline]
    pub const fn window(&self) -> WeilSupportWindow {
        self.window
    }

    #[inline]
    pub fn samples(&self) -> &[WeilBernsteinEvidenceSample] {
        &self.samples
    }

    #[inline]
    pub const fn bernstein_observed_interval(&self) -> (f64, f64) {
        (
            self.bernstein_observed_minimum,
            self.bernstein_observed_maximum,
        )
    }

    #[inline]
    pub const fn leading_observed_interval(&self) -> (f64, f64) {
        (self.leading_observed_minimum, self.leading_observed_maximum)
    }

    #[inline]
    pub const fn family_delta_observed_interval(&self) -> (f64, f64) {
        (
            self.family_delta_observed_minimum,
            self.family_delta_observed_maximum,
        )
    }

    #[inline]
    pub const fn bernstein_observed_span(&self) -> f64 {
        self.bernstein_observed_maximum - self.bernstein_observed_minimum
    }

    #[inline]
    pub const fn leading_observed_span(&self) -> f64 {
        self.leading_observed_maximum - self.leading_observed_minimum
    }

    #[inline]
    pub const fn family_delta_observed_span(&self) -> f64 {
        self.family_delta_observed_maximum - self.family_delta_observed_minimum
    }

    pub fn last_bernstein_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some(
            (last.bernstein_generalized_minimum_eigenvalue
                - previous.bernstein_generalized_minimum_eigenvalue)
                .abs(),
        )
    }

    pub fn last_leading_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some(
            (last.leading_legendre_generalized_minimum_eigenvalue
                - previous.leading_legendre_generalized_minimum_eigenvalue)
                .abs(),
        )
    }

    pub fn last_family_delta_change(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some((last.generalized_family_delta - previous.generalized_family_delta).abs())
    }

    pub fn maximum_bernstein_gram_condition_number(&self) -> f64 {
        self.samples
            .iter()
            .map(|sample| sample.bernstein_gram_condition_number)
            .fold(0.0_f64, f64::max)
    }

    pub fn maximum_leading_gram_condition_number(&self) -> f64 {
        self.samples
            .iter()
            .map(|sample| sample.leading_legendre_gram_condition_number)
            .fold(0.0_f64, f64::max)
    }
}

/// One fixed Bernstein subspace crossed with support windows and refinement levels.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilBernsteinEvidenceGrid {
    subspace: BernsteinIndexSubspace,
    windows: Vec<WeilSupportWindow>,
    levels: Vec<WeilQuadratureLevel>,
    cells: Vec<FiniteWeilBernsteinEvidenceCell>,
}

impl FiniteWeilBernsteinEvidenceGrid {
    #[inline]
    pub fn subspace(&self) -> &BernsteinIndexSubspace {
        &self.subspace
    }

    #[inline]
    pub fn windows(&self) -> &[WeilSupportWindow] {
        &self.windows
    }

    #[inline]
    pub fn levels(&self) -> &[WeilQuadratureLevel] {
        &self.levels
    }

    #[inline]
    pub fn cells(&self) -> &[FiniteWeilBernsteinEvidenceCell] {
        &self.cells
    }

    #[inline]
    pub fn cell(&self, window_index: usize) -> Option<&FiniteWeilBernsteinEvidenceCell> {
        self.cells.get(window_index)
    }
}

#[derive(Debug)]
pub enum FiniteWeilBernsteinEvidenceGridError {
    EmptyWindowSet,
    EmptyLevelSet,
    CompactSupport(CompactArchimedeanError),
    BernsteinSubspace(FiniteWeilBernsteinSubspaceError),
}

impl fmt::Display for FiniteWeilBernsteinEvidenceGridError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyWindowSet => write!(
                f,
                "Bernstein finite Weil evidence grid requires at least one support window"
            ),
            Self::EmptyLevelSet => write!(
                f,
                "Bernstein finite Weil evidence grid requires at least one quadrature level"
            ),
            Self::CompactSupport(error) => write!(f, "invalid compact support window: {error}"),
            Self::BernsteinSubspace(error) => {
                write!(f, "Bernstein finite Weil subspace audit failed: {error}")
            }
        }
    }
}

impl std::error::Error for FiniteWeilBernsteinEvidenceGridError {}

impl From<CompactArchimedeanError> for FiniteWeilBernsteinEvidenceGridError {
    fn from(value: CompactArchimedeanError) -> Self {
        Self::CompactSupport(value)
    }
}

impl From<FiniteWeilBernsteinSubspaceError> for FiniteWeilBernsteinEvidenceGridError {
    fn from(value: FiniteWeilBernsteinSubspaceError) -> Self {
        Self::BernsteinSubspace(value)
    }
}

/// Evaluate one fixed localized Bernstein subspace against its dimension-matched
/// leading-Legendre control over all declared support/refinement cells.
pub fn audit_finite_weil_bernstein_evidence_grid(
    windows: &[WeilSupportWindow],
    subspace: &BernsteinIndexSubspace,
    levels: &[WeilQuadratureLevel],
) -> Result<FiniteWeilBernsteinEvidenceGrid, FiniteWeilBernsteinEvidenceGridError> {
    if windows.is_empty() {
        return Err(FiniteWeilBernsteinEvidenceGridError::EmptyWindowSet);
    }
    if levels.is_empty() {
        return Err(FiniteWeilBernsteinEvidenceGridError::EmptyLevelSet);
    }

    let mut cells = Vec::with_capacity(windows.len());
    for &window in windows {
        let bump = CompactArchimedeanBump::new(window.lower(), window.upper())?;
        let mut samples = Vec::with_capacity(levels.len());

        for &level in levels {
            let audit = audit_finite_weil_bernstein_subspace(
                bump,
                subspace,
                level.correlation_order(),
                level.archimedean_order(),
                level.boundary_order(),
                level.gram_order(),
            )?;
            let bernstein_generalized_minimum_eigenvalue = audit.minimum_generalized_eigenvalue();
            let leading_legendre_generalized_minimum_eigenvalue =
                audit.leading_legendre_generalized_minimum();
            let generalized_family_delta = bernstein_generalized_minimum_eigenvalue
                - leading_legendre_generalized_minimum_eigenvalue;

            samples.push(WeilBernsteinEvidenceSample {
                level,
                bernstein_raw_minimum_eigenvalue: audit.minimum_raw_eigenvalue(),
                bernstein_generalized_minimum_eigenvalue,
                leading_legendre_generalized_minimum_eigenvalue,
                generalized_family_delta,
                bernstein_gram_condition_number: audit.gram_condition_number(),
                leading_legendre_gram_condition_number: audit
                    .leading_legendre_gram_condition_number(),
                max_boundary_residual: audit.max_boundary_residual(),
                max_pairing_asymmetry: audit.parent_max_raw_pairing_asymmetry(),
                max_whitened_asymmetry: audit.max_whitened_asymmetry(),
            });
        }

        let bernstein_observed_minimum = samples
            .iter()
            .map(|sample| sample.bernstein_generalized_minimum_eigenvalue)
            .fold(f64::INFINITY, f64::min);
        let bernstein_observed_maximum = samples
            .iter()
            .map(|sample| sample.bernstein_generalized_minimum_eigenvalue)
            .fold(f64::NEG_INFINITY, f64::max);
        let leading_observed_minimum = samples
            .iter()
            .map(|sample| sample.leading_legendre_generalized_minimum_eigenvalue)
            .fold(f64::INFINITY, f64::min);
        let leading_observed_maximum = samples
            .iter()
            .map(|sample| sample.leading_legendre_generalized_minimum_eigenvalue)
            .fold(f64::NEG_INFINITY, f64::max);
        let family_delta_observed_minimum = samples
            .iter()
            .map(|sample| sample.generalized_family_delta)
            .fold(f64::INFINITY, f64::min);
        let family_delta_observed_maximum = samples
            .iter()
            .map(|sample| sample.generalized_family_delta)
            .fold(f64::NEG_INFINITY, f64::max);

        cells.push(FiniteWeilBernsteinEvidenceCell {
            window,
            samples,
            bernstein_observed_minimum,
            bernstein_observed_maximum,
            leading_observed_minimum,
            leading_observed_maximum,
            family_delta_observed_minimum,
            family_delta_observed_maximum,
        });
    }

    Ok(FiniteWeilBernsteinEvidenceGrid {
        subspace: subspace.clone(),
        windows: windows.to_vec(),
        levels: levels.to_vec(),
        cells,
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
    fn grid_rejects_empty_axes() {
        let subspace = BernsteinIndexSubspace::new(2, vec![0, 2]).unwrap();
        let level = WeilQuadratureLevel::new(16, 16, 24, 24);
        assert!(matches!(
            audit_finite_weil_bernstein_evidence_grid(&[], &subspace, &[level]),
            Err(FiniteWeilBernsteinEvidenceGridError::EmptyWindowSet)
        ));
        assert!(matches!(
            audit_finite_weil_bernstein_evidence_grid(&[window()], &subspace, &[]),
            Err(FiniteWeilBernsteinEvidenceGridError::EmptyLevelSet)
        ));
    }

    #[test]
    fn grid_crosses_one_localized_subspace_with_refinement_without_sign_assumption() {
        let subspace = BernsteinIndexSubspace::new(2, vec![0, 2]).unwrap();
        let levels = [
            WeilQuadratureLevel::new(16, 16, 24, 24),
            WeilQuadratureLevel::new(20, 20, 28, 28),
        ];
        let grid =
            audit_finite_weil_bernstein_evidence_grid(&[window()], &subspace, &levels).unwrap();

        assert_eq!(grid.subspace(), &subspace);
        assert_eq!(grid.windows().len(), 1);
        assert_eq!(grid.levels(), &levels);
        assert_eq!(grid.cells().len(), 1);

        let cell = grid.cell(0).unwrap();
        assert_eq!(cell.window(), window());
        assert_eq!(cell.samples().len(), 2);
        assert!(cell.bernstein_observed_span().is_finite());
        assert!(cell.bernstein_observed_span() >= 0.0);
        assert!(cell.leading_observed_span().is_finite());
        assert!(cell.leading_observed_span() >= 0.0);
        assert!(cell.family_delta_observed_span().is_finite());
        assert!(cell.family_delta_observed_span() >= 0.0);
        assert!(cell.last_bernstein_delta().unwrap().is_finite());
        assert!(cell.last_leading_delta().unwrap().is_finite());
        assert!(cell.last_family_delta_change().unwrap().is_finite());
        assert!(cell.maximum_bernstein_gram_condition_number().is_finite());
        assert!(cell.maximum_leading_gram_condition_number().is_finite());

        for sample in cell.samples() {
            assert!(sample.bernstein_raw_minimum_eigenvalue().is_finite());
            assert!(
                sample
                    .bernstein_generalized_minimum_eigenvalue()
                    .is_finite()
            );
            assert!(
                sample
                    .leading_legendre_generalized_minimum_eigenvalue()
                    .is_finite()
            );
            assert!(sample.generalized_family_delta().is_finite());
            assert!(sample.bernstein_gram_condition_number().is_finite());
            assert!(sample.leading_legendre_gram_condition_number().is_finite());
            assert!(sample.max_boundary_residual().is_finite());
            assert!(sample.max_pairing_asymmetry().is_finite());
            assert!(sample.max_whitened_asymmetry().is_finite());
        }
    }
}
