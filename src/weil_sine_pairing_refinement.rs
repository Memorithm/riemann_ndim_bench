//! Quadrature-refinement grid for the full direct-vs-projected sine Weil pairing.
//!
//! The direct sine pairing audit already varies the Legendre parent dimension.
//! This module crosses that approximation axis with an explicit sequence of
//! coefficient and Weil-pairing quadrature levels. The resulting spreads and
//! consecutive deltas are empirical resolution diagnostics only; they are not
//! certified error bounds, confidence intervals, or evidence of complete-space
//! Weil positivity.

use std::fmt;

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_compact_pairing::CompactWeilPairingConfig;
use crate::weil_sine_pairing::{
    DirectSinePairingAuditConfig, DirectSinePairingSample, FiniteWeilDirectSinePairingError,
    audit_finite_weil_direct_sine_pairing,
};
use crate::weil_sine_truncation::SineModeSet;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SinePairingRefinementLevel {
    coefficient_quadrature_order: usize,
    correlation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
}

impl SinePairingRefinementLevel {
    #[inline]
    pub const fn new(
        coefficient_quadrature_order: usize,
        correlation_order: usize,
        archimedean_order: usize,
        boundary_order: usize,
    ) -> Self {
        Self {
            coefficient_quadrature_order,
            correlation_order,
            archimedean_order,
            boundary_order,
        }
    }

    #[inline]
    pub const fn coefficient_quadrature_order(self) -> usize {
        self.coefficient_quadrature_order
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
    fn audit_config(self) -> DirectSinePairingAuditConfig {
        DirectSinePairingAuditConfig::new(
            self.coefficient_quadrature_order,
            CompactWeilPairingConfig::new(
                self.correlation_order,
                self.archimedean_order,
                self.boundary_order,
            ),
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SinePairingRefinementSample {
    level: SinePairingRefinementLevel,
    max_pairing_amplitude: f64,
    max_direct_projected_pairing_residual: f64,
    max_normalized_pairing_residual: f64,
    max_parent_matrix_projection_residual: f64,
    max_pole_term_residual: f64,
    max_archimedean_term_residual: f64,
    max_prime_total_residual: f64,
    max_direct_pairing_asymmetry: f64,
    max_projected_pairing_asymmetry: f64,
    max_direct_boundary_residual: f64,
    max_projected_boundary_residual: f64,
}

impl SinePairingRefinementSample {
    #[inline]
    pub const fn level(self) -> SinePairingRefinementLevel {
        self.level
    }

    #[inline]
    pub const fn max_pairing_amplitude(self) -> f64 {
        self.max_pairing_amplitude
    }

    #[inline]
    pub const fn max_direct_projected_pairing_residual(self) -> f64 {
        self.max_direct_projected_pairing_residual
    }

    #[inline]
    pub const fn max_normalized_pairing_residual(self) -> f64 {
        self.max_normalized_pairing_residual
    }

    #[inline]
    pub const fn max_parent_matrix_projection_residual(self) -> f64 {
        self.max_parent_matrix_projection_residual
    }

    #[inline]
    pub const fn max_pole_term_residual(self) -> f64 {
        self.max_pole_term_residual
    }

    #[inline]
    pub const fn max_archimedean_term_residual(self) -> f64 {
        self.max_archimedean_term_residual
    }

    #[inline]
    pub const fn max_prime_total_residual(self) -> f64 {
        self.max_prime_total_residual
    }

    #[inline]
    pub const fn max_direct_pairing_asymmetry(self) -> f64 {
        self.max_direct_pairing_asymmetry
    }

    #[inline]
    pub const fn max_projected_pairing_asymmetry(self) -> f64 {
        self.max_projected_pairing_asymmetry
    }

    #[inline]
    pub const fn max_direct_boundary_residual(self) -> f64 {
        self.max_direct_boundary_residual
    }

    #[inline]
    pub const fn max_projected_boundary_residual(self) -> f64 {
        self.max_projected_boundary_residual
    }
}

impl From<(SinePairingRefinementLevel, DirectSinePairingSample)> for SinePairingRefinementSample {
    fn from((level, sample): (SinePairingRefinementLevel, DirectSinePairingSample)) -> Self {
        Self {
            level,
            max_pairing_amplitude: sample.max_pairing_amplitude(),
            max_direct_projected_pairing_residual: sample.max_direct_projected_pairing_residual(),
            max_normalized_pairing_residual: sample.max_normalized_pairing_residual(),
            max_parent_matrix_projection_residual: sample.max_parent_matrix_projection_residual(),
            max_pole_term_residual: sample.max_pole_term_residual(),
            max_archimedean_term_residual: sample.max_archimedean_term_residual(),
            max_prime_total_residual: sample.max_prime_total_residual(),
            max_direct_pairing_asymmetry: sample.max_direct_pairing_asymmetry(),
            max_projected_pairing_asymmetry: sample.max_projected_pairing_asymmetry(),
            max_direct_boundary_residual: sample.max_direct_boundary_residual(),
            max_projected_boundary_residual: sample.max_projected_boundary_residual(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilSinePairingRefinementCell {
    parent_dimension: usize,
    samples: Vec<SinePairingRefinementSample>,
    pairing_residual_observed_minimum: f64,
    pairing_residual_observed_maximum: f64,
    normalized_residual_observed_minimum: f64,
    normalized_residual_observed_maximum: f64,
}

impl FiniteWeilSinePairingRefinementCell {
    #[inline]
    pub const fn parent_dimension(&self) -> usize {
        self.parent_dimension
    }

    #[inline]
    pub fn samples(&self) -> &[SinePairingRefinementSample] {
        &self.samples
    }

    #[inline]
    pub const fn pairing_residual_observed_interval(&self) -> (f64, f64) {
        (
            self.pairing_residual_observed_minimum,
            self.pairing_residual_observed_maximum,
        )
    }

    #[inline]
    pub const fn pairing_residual_observed_span(&self) -> f64 {
        self.pairing_residual_observed_maximum - self.pairing_residual_observed_minimum
    }

    #[inline]
    pub const fn normalized_residual_observed_interval(&self) -> (f64, f64) {
        (
            self.normalized_residual_observed_minimum,
            self.normalized_residual_observed_maximum,
        )
    }

    #[inline]
    pub const fn normalized_residual_observed_span(&self) -> f64 {
        self.normalized_residual_observed_maximum - self.normalized_residual_observed_minimum
    }

    pub fn last_pairing_residual_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some(
            (last.max_direct_projected_pairing_residual
                - previous.max_direct_projected_pairing_residual)
                .abs(),
        )
    }

    pub fn last_normalized_residual_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some(
            (last.max_normalized_pairing_residual - previous.max_normalized_pairing_residual).abs(),
        )
    }

    pub fn max_parent_matrix_projection_residual(&self) -> f64 {
        self.samples
            .iter()
            .map(|sample| sample.max_parent_matrix_projection_residual)
            .fold(0.0_f64, f64::max)
    }

    pub fn max_pairing_asymmetry(&self) -> f64 {
        self.samples
            .iter()
            .map(|sample| {
                sample
                    .max_direct_pairing_asymmetry
                    .max(sample.max_projected_pairing_asymmetry)
            })
            .fold(0.0_f64, f64::max)
    }

    pub fn max_boundary_residual(&self) -> f64 {
        self.samples
            .iter()
            .map(|sample| {
                sample
                    .max_direct_boundary_residual
                    .max(sample.max_projected_boundary_residual)
            })
            .fold(0.0_f64, f64::max)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilSinePairingRefinementGrid {
    modes: SineModeSet,
    parent_dimensions: Vec<usize>,
    levels: Vec<SinePairingRefinementLevel>,
    cells: Vec<FiniteWeilSinePairingRefinementCell>,
}

impl FiniteWeilSinePairingRefinementGrid {
    #[inline]
    pub fn modes(&self) -> &SineModeSet {
        &self.modes
    }

    #[inline]
    pub fn parent_dimensions(&self) -> &[usize] {
        &self.parent_dimensions
    }

    #[inline]
    pub fn levels(&self) -> &[SinePairingRefinementLevel] {
        &self.levels
    }

    #[inline]
    pub fn cells(&self) -> &[FiniteWeilSinePairingRefinementCell] {
        &self.cells
    }

    pub fn cell(&self, parent_dimension: usize) -> Option<&FiniteWeilSinePairingRefinementCell> {
        self.cells
            .iter()
            .find(|cell| cell.parent_dimension == parent_dimension)
    }
}

#[derive(Debug)]
pub enum FiniteWeilSinePairingRefinementError {
    EmptyLevelSet,
    DirectPairing(FiniteWeilDirectSinePairingError),
}

impl fmt::Display for FiniteWeilSinePairingRefinementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyLevelSet => write!(
                f,
                "finite sine pairing refinement grid requires at least one quadrature level"
            ),
            Self::DirectPairing(error) => {
                write!(
                    f,
                    "finite sine pairing refinement evaluation failed: {error}"
                )
            }
        }
    }
}

impl std::error::Error for FiniteWeilSinePairingRefinementError {}

impl From<FiniteWeilDirectSinePairingError> for FiniteWeilSinePairingRefinementError {
    fn from(value: FiniteWeilDirectSinePairingError) -> Self {
        Self::DirectPairing(value)
    }
}

pub fn audit_finite_weil_sine_pairing_refinement(
    bump: CompactArchimedeanBump,
    modes: &SineModeSet,
    parent_dimensions: &[usize],
    levels: &[SinePairingRefinementLevel],
) -> Result<FiniteWeilSinePairingRefinementGrid, FiniteWeilSinePairingRefinementError> {
    if levels.is_empty() {
        return Err(FiniteWeilSinePairingRefinementError::EmptyLevelSet);
    }

    let mut samples_by_dimension = vec![Vec::with_capacity(levels.len()); parent_dimensions.len()];
    for &level in levels {
        let audit = audit_finite_weil_direct_sine_pairing(
            bump,
            modes,
            parent_dimensions,
            level.audit_config(),
        )?;
        for (dimension_index, &sample) in audit.samples().iter().enumerate() {
            samples_by_dimension[dimension_index].push((level, sample).into());
        }
    }

    let cells = parent_dimensions
        .iter()
        .copied()
        .zip(samples_by_dimension)
        .map(|(parent_dimension, samples)| build_cell(parent_dimension, samples))
        .collect::<Vec<_>>();

    Ok(FiniteWeilSinePairingRefinementGrid {
        modes: modes.clone(),
        parent_dimensions: parent_dimensions.to_vec(),
        levels: levels.to_vec(),
        cells,
    })
}

fn build_cell(
    parent_dimension: usize,
    samples: Vec<SinePairingRefinementSample>,
) -> FiniteWeilSinePairingRefinementCell {
    let pairing_residual_observed_minimum = samples
        .iter()
        .map(|sample| sample.max_direct_projected_pairing_residual)
        .min_by(f64::total_cmp)
        .expect("level set is validated non-empty");
    let pairing_residual_observed_maximum = samples
        .iter()
        .map(|sample| sample.max_direct_projected_pairing_residual)
        .max_by(f64::total_cmp)
        .expect("level set is validated non-empty");
    let normalized_residual_observed_minimum = samples
        .iter()
        .map(|sample| sample.max_normalized_pairing_residual)
        .min_by(f64::total_cmp)
        .expect("level set is validated non-empty");
    let normalized_residual_observed_maximum = samples
        .iter()
        .map(|sample| sample.max_normalized_pairing_residual)
        .max_by(f64::total_cmp)
        .expect("level set is validated non-empty");

    FiniteWeilSinePairingRefinementCell {
        parent_dimension,
        samples,
        pairing_residual_observed_minimum,
        pairing_residual_observed_maximum,
        normalized_residual_observed_minimum,
        normalized_residual_observed_maximum,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::PositiveRational;

    fn bump() -> CompactArchimedeanBump {
        CompactArchimedeanBump::new(
            PositiveRational::new(1, 2).unwrap(),
            PositiveRational::new(7, 2).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn refinement_grid_crosses_parent_dimension_and_quadrature_without_sign_assumption() {
        let modes = SineModeSet::new(vec![1, 2]).unwrap();
        let levels = [
            SinePairingRefinementLevel::new(32, 16, 16, 24),
            SinePairingRefinementLevel::new(40, 20, 20, 28),
        ];
        let grid =
            audit_finite_weil_sine_pairing_refinement(bump(), &modes, &[3, 5], &levels).unwrap();

        assert_eq!(grid.modes(), &modes);
        assert_eq!(grid.parent_dimensions(), &[3, 5]);
        assert_eq!(grid.levels(), &levels);
        assert_eq!(grid.cells().len(), 2);
        for cell in grid.cells() {
            assert_eq!(cell.samples().len(), 2);
            assert!(cell.pairing_residual_observed_span().is_finite());
            assert!(cell.pairing_residual_observed_span() >= 0.0);
            assert!(cell.normalized_residual_observed_span().is_finite());
            assert!(cell.normalized_residual_observed_span() >= 0.0);
            assert!(cell.last_pairing_residual_delta().unwrap().is_finite());
            assert!(cell.last_normalized_residual_delta().unwrap().is_finite());
            assert!(cell.max_parent_matrix_projection_residual().is_finite());
            assert!(cell.max_pairing_asymmetry().is_finite());
            assert!(cell.max_boundary_residual().is_finite());
            for sample in cell.samples().iter().copied() {
                assert!(sample.max_pairing_amplitude().is_finite());
                assert!(sample.max_direct_projected_pairing_residual().is_finite());
                assert!(sample.max_normalized_pairing_residual().is_finite());
                assert!(sample.max_parent_matrix_projection_residual().is_finite());
                assert!(sample.max_pole_term_residual().is_finite());
                assert!(sample.max_archimedean_term_residual().is_finite());
                assert!(sample.max_prime_total_residual().is_finite());
                assert!(sample.max_direct_pairing_asymmetry().is_finite());
                assert!(sample.max_projected_pairing_asymmetry().is_finite());
                assert!(sample.max_direct_boundary_residual().is_finite());
                assert!(sample.max_projected_boundary_residual().is_finite());
            }
        }
    }

    #[test]
    fn refinement_grid_rejects_empty_level_axis() {
        let modes = SineModeSet::new(vec![1]).unwrap();
        assert!(matches!(
            audit_finite_weil_sine_pairing_refinement(bump(), &modes, &[2], &[]),
            Err(FiniteWeilSinePairingRefinementError::EmptyLevelSet)
        ));
    }
}
