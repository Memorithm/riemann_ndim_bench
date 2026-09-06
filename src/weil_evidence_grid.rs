//! Consolidated support x dimension x quadrature grid for the finite Weil audit.
//!
//! For each support window and quadrature level, one maximum-dimensional
//! Gram-normalized Weil audit is computed. All smaller dimensions are extracted
//! as leading principal subproblems. Results are then regrouped by
//! `(support, dimension)` across refinement levels.
//!
//! This module records numerical evidence without assigning a statistical or
//! theorem-level significance score. Observed spreads and last-step deltas are
//! empirical resolution diagnostics only.

use std::fmt;

use crate::semilocal_compact_archimedean::{CompactArchimedeanBump, CompactArchimedeanError};
use crate::weil_generalized_spectrum::{
    FiniteWeilGeneralizedSpectrumError, audit_finite_weil_generalized_spectrum,
};
use crate::weil_refinement::WeilQuadratureLevel;
use crate::weil_support_sweep::WeilSupportWindow;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WeilEvidenceSample {
    level: WeilQuadratureLevel,
    raw_minimum_eigenvalue: f64,
    generalized_minimum_eigenvalue: f64,
    gram_condition_number: f64,
    max_boundary_residual: f64,
    max_pairing_asymmetry: f64,
    max_whitened_asymmetry: f64,
}

impl WeilEvidenceSample {
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

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilEvidenceCell {
    window: WeilSupportWindow,
    dimension: usize,
    samples: Vec<WeilEvidenceSample>,
    raw_observed_minimum: f64,
    raw_observed_maximum: f64,
    generalized_observed_minimum: f64,
    generalized_observed_maximum: f64,
}

impl FiniteWeilEvidenceCell {
    #[inline]
    pub const fn window(&self) -> WeilSupportWindow {
        self.window
    }

    #[inline]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    #[inline]
    pub fn samples(&self) -> &[WeilEvidenceSample] {
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
        Some((last.generalized_minimum_eigenvalue - previous.generalized_minimum_eigenvalue).abs())
    }

    pub fn maximum_gram_condition_number(&self) -> f64 {
        self.samples
            .iter()
            .map(|sample| sample.gram_condition_number)
            .fold(0.0_f64, f64::max)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilEvidenceGrid {
    max_dimension: usize,
    windows: Vec<WeilSupportWindow>,
    levels: Vec<WeilQuadratureLevel>,
    cells: Vec<FiniteWeilEvidenceCell>,
}

impl FiniteWeilEvidenceGrid {
    #[inline]
    pub const fn max_dimension(&self) -> usize {
        self.max_dimension
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
    pub fn cells(&self) -> &[FiniteWeilEvidenceCell] {
        &self.cells
    }

    pub fn cell(&self, window_index: usize, dimension: usize) -> Option<&FiniteWeilEvidenceCell> {
        if window_index >= self.windows.len() || dimension == 0 || dimension > self.max_dimension {
            return None;
        }
        self.cells
            .get(window_index * self.max_dimension + (dimension - 1))
    }
}

#[derive(Debug)]
pub enum FiniteWeilEvidenceGridError {
    EmptyWindowSet,
    EmptyLevelSet,
    CompactSupport(CompactArchimedeanError),
    Generalized(FiniteWeilGeneralizedSpectrumError),
}

impl fmt::Display for FiniteWeilEvidenceGridError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyWindowSet => write!(
                f,
                "finite Weil evidence grid requires at least one support window"
            ),
            Self::EmptyLevelSet => write!(
                f,
                "finite Weil evidence grid requires at least one quadrature level"
            ),
            Self::CompactSupport(error) => write!(f, "invalid compact support window: {error}"),
            Self::Generalized(error) => write!(f, "generalized finite Weil audit failed: {error}"),
        }
    }
}

impl std::error::Error for FiniteWeilEvidenceGridError {}

impl From<CompactArchimedeanError> for FiniteWeilEvidenceGridError {
    fn from(value: CompactArchimedeanError) -> Self {
        Self::CompactSupport(value)
    }
}

impl From<FiniteWeilGeneralizedSpectrumError> for FiniteWeilEvidenceGridError {
    fn from(value: FiniteWeilGeneralizedSpectrumError) -> Self {
        Self::Generalized(value)
    }
}

pub fn audit_finite_weil_evidence_grid(
    windows: &[WeilSupportWindow],
    max_dimension: usize,
    levels: &[WeilQuadratureLevel],
) -> Result<FiniteWeilEvidenceGrid, FiniteWeilEvidenceGridError> {
    if windows.is_empty() {
        return Err(FiniteWeilEvidenceGridError::EmptyWindowSet);
    }
    if levels.is_empty() {
        return Err(FiniteWeilEvidenceGridError::EmptyLevelSet);
    }

    let mut cells = Vec::with_capacity(windows.len() * max_dimension);
    for &window in windows {
        let bump = CompactArchimedeanBump::new(window.lower(), window.upper())?;
        let mut samples_by_dimension = vec![Vec::with_capacity(levels.len()); max_dimension];

        for &level in levels {
            let audit = audit_finite_weil_generalized_spectrum(
                bump,
                max_dimension,
                level.correlation_order(),
                level.archimedean_order(),
                level.boundary_order(),
                level.gram_order(),
            )?;
            let rows = audit.principal_sweep()?;
            let max_boundary_residual = audit.pairing().max_boundary_residual();
            let max_pairing_asymmetry = audit.pairing().max_raw_pairing_asymmetry();
            let max_whitened_asymmetry = audit.max_whitened_asymmetry();

            for row in rows {
                samples_by_dimension[row.dimension() - 1].push(WeilEvidenceSample {
                    level,
                    raw_minimum_eigenvalue: row.raw_minimum_eigenvalue(),
                    generalized_minimum_eigenvalue: row.generalized_minimum_eigenvalue(),
                    gram_condition_number: row.gram_condition_number(),
                    max_boundary_residual,
                    max_pairing_asymmetry,
                    max_whitened_asymmetry,
                });
            }
        }

        for (offset, samples) in samples_by_dimension.into_iter().enumerate() {
            let raw_observed_minimum = samples
                .iter()
                .map(|sample| sample.raw_minimum_eigenvalue)
                .min_by(f64::total_cmp)
                .expect("non-empty quadrature level set");
            let raw_observed_maximum = samples
                .iter()
                .map(|sample| sample.raw_minimum_eigenvalue)
                .max_by(f64::total_cmp)
                .expect("non-empty quadrature level set");
            let generalized_observed_minimum = samples
                .iter()
                .map(|sample| sample.generalized_minimum_eigenvalue)
                .min_by(f64::total_cmp)
                .expect("non-empty quadrature level set");
            let generalized_observed_maximum = samples
                .iter()
                .map(|sample| sample.generalized_minimum_eigenvalue)
                .max_by(f64::total_cmp)
                .expect("non-empty quadrature level set");

            cells.push(FiniteWeilEvidenceCell {
                window,
                dimension: offset + 1,
                samples,
                raw_observed_minimum,
                raw_observed_maximum,
                generalized_observed_minimum,
                generalized_observed_maximum,
            });
        }
    }

    Ok(FiniteWeilEvidenceGrid {
        max_dimension,
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

    #[test]
    fn evidence_grid_rejects_empty_axes() {
        let window = WeilSupportWindow::new(rational(1, 2), rational(7, 2));
        let level = WeilQuadratureLevel::new(16, 16, 24, 24);
        assert!(matches!(
            audit_finite_weil_evidence_grid(&[], 1, &[level]).unwrap_err(),
            FiniteWeilEvidenceGridError::EmptyWindowSet
        ));
        assert!(matches!(
            audit_finite_weil_evidence_grid(&[window], 1, &[]).unwrap_err(),
            FiniteWeilEvidenceGridError::EmptyLevelSet
        ));
    }

    #[test]
    fn evidence_grid_reuses_max_dimension_across_levels() {
        let windows = [WeilSupportWindow::new(rational(1, 2), rational(7, 2))];
        let levels = [
            WeilQuadratureLevel::new(20, 20, 28, 28),
            WeilQuadratureLevel::new(24, 24, 32, 32),
        ];
        let grid = audit_finite_weil_evidence_grid(&windows, 2, &levels).unwrap();

        assert_eq!(grid.max_dimension(), 2);
        assert_eq!(grid.windows(), windows);
        assert_eq!(grid.levels(), levels);
        assert_eq!(grid.cells().len(), 2);
        for dimension in 1..=2 {
            let cell = grid.cell(0, dimension).unwrap();
            assert_eq!(cell.dimension(), dimension);
            assert_eq!(cell.samples().len(), 2);
            assert!(cell.raw_observed_span().is_finite());
            assert!(cell.raw_observed_span() >= 0.0);
            assert!(cell.generalized_observed_span().is_finite());
            assert!(cell.generalized_observed_span() >= 0.0);
            assert!(cell.last_raw_delta().unwrap().is_finite());
            assert!(cell.last_generalized_delta().unwrap().is_finite());
            assert!(cell.maximum_gram_condition_number().is_finite());
            for sample in cell.samples() {
                assert!(sample.raw_minimum_eigenvalue().is_finite());
                assert!(sample.generalized_minimum_eigenvalue().is_finite());
                assert!(sample.gram_condition_number().is_finite());
                assert!(sample.max_boundary_residual().is_finite());
                assert!(sample.max_pairing_asymmetry().is_finite());
                assert!(sample.max_whitened_asymmetry().is_finite());
            }
        }
    }
}
