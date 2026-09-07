//! Legendre-truncation audit for non-polynomial sine-enriched compact generators.
//!
//! The target enrichments are
//!
//! `g_m(rho) = bump(rho) sin(m pi t)`, `m >= 1`,
//!
//! with the same affine support coordinate `t` used by the existing compact
//! Legendre family.  Rather than introducing a second Riemann--Weil pairing
//! implementation, each sine enrichment is projected onto a finite shifted-
//! Legendre parent basis.  The parent dimension is then increased explicitly.
//!
//! Every finite truncation still lies inside a finite polynomial parent space.
//! Convergence of reconstruction and spectral diagnostics across increasing
//! parent dimensions is numerical evidence about the approximation only; it is
//! not an exact evaluation of the infinite Legendre series, a density theorem,
//! complete-space Weil positivity, Conjecture 4.1, or RH.

use std::f64::consts::PI;
use std::fmt;

use crate::quadrature::{GaussLegendreUnit, QuadratureError};
use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_coefficient_subspace::{
    FiniteWeilCoefficientSubspaceError, audit_finite_weil_coefficient_subspace,
};
use crate::weil_generalized_spectrum::{
    FiniteWeilGeneralizedSpectrumError, audit_finite_weil_generalized_spectrum,
};
use crate::weil_refinement::WeilQuadratureLevel;

const RECONSTRUCTION_INTERVALS: usize = 256;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SineModeSet {
    modes: Vec<usize>,
}

impl SineModeSet {
    pub fn new(modes: Vec<usize>) -> Result<Self, FiniteWeilSineTruncationError> {
        if modes.is_empty() {
            return Err(FiniteWeilSineTruncationError::EmptyModeSet);
        }
        for &mode in &modes {
            if mode == 0 {
                return Err(FiniteWeilSineTruncationError::ZeroMode);
            }
        }
        for pair in modes.windows(2) {
            if pair[0] >= pair[1] {
                return Err(FiniteWeilSineTruncationError::ModesNotStrictlyIncreasing {
                    previous: pair[0],
                    next: pair[1],
                });
            }
        }
        Ok(Self { modes })
    }

    #[inline]
    pub fn modes(&self) -> &[usize] {
        &self.modes
    }

    #[inline]
    pub fn dimension(&self) -> usize {
        self.modes.len()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SineTruncationSample {
    parent_dimension: usize,
    max_reconstruction_residual: f64,
    max_coefficient_l1_norm: f64,
    raw_minimum_eigenvalue: f64,
    generalized_minimum_eigenvalue: f64,
    leading_legendre_generalized_minimum_eigenvalue: f64,
    generalized_family_delta: f64,
    gram_condition_number: f64,
    leading_legendre_gram_condition_number: f64,
    max_boundary_residual: f64,
    max_pairing_asymmetry: f64,
    max_whitened_asymmetry: f64,
}

impl SineTruncationSample {
    #[inline]
    pub const fn parent_dimension(self) -> usize {
        self.parent_dimension
    }

    #[inline]
    pub const fn max_reconstruction_residual(self) -> f64 {
        self.max_reconstruction_residual
    }

    #[inline]
    pub const fn max_coefficient_l1_norm(self) -> f64 {
        self.max_coefficient_l1_norm
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
    pub const fn leading_legendre_generalized_minimum_eigenvalue(self) -> f64 {
        self.leading_legendre_generalized_minimum_eigenvalue
    }

    #[inline]
    pub const fn generalized_family_delta(self) -> f64 {
        self.generalized_family_delta
    }

    #[inline]
    pub const fn gram_condition_number(self) -> f64 {
        self.gram_condition_number
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

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilSineTruncationAudit {
    modes: SineModeSet,
    parent_dimensions: Vec<usize>,
    coefficient_quadrature_order: usize,
    weil_level: WeilQuadratureLevel,
    samples: Vec<SineTruncationSample>,
    generalized_observed_minimum: f64,
    generalized_observed_maximum: f64,
    family_delta_observed_minimum: f64,
    family_delta_observed_maximum: f64,
}

impl FiniteWeilSineTruncationAudit {
    #[inline]
    pub fn modes(&self) -> &SineModeSet {
        &self.modes
    }

    #[inline]
    pub fn parent_dimensions(&self) -> &[usize] {
        &self.parent_dimensions
    }

    #[inline]
    pub const fn coefficient_quadrature_order(&self) -> usize {
        self.coefficient_quadrature_order
    }

    #[inline]
    pub const fn weil_level(&self) -> WeilQuadratureLevel {
        self.weil_level
    }

    #[inline]
    pub fn samples(&self) -> &[SineTruncationSample] {
        &self.samples
    }

    #[inline]
    pub const fn generalized_observed_interval(&self) -> (f64, f64) {
        (
            self.generalized_observed_minimum,
            self.generalized_observed_maximum,
        )
    }

    #[inline]
    pub const fn generalized_observed_span(&self) -> f64 {
        self.generalized_observed_maximum - self.generalized_observed_minimum
    }

    #[inline]
    pub const fn family_delta_observed_interval(&self) -> (f64, f64) {
        (
            self.family_delta_observed_minimum,
            self.family_delta_observed_maximum,
        )
    }

    #[inline]
    pub const fn family_delta_observed_span(&self) -> f64 {
        self.family_delta_observed_maximum - self.family_delta_observed_minimum
    }

    pub fn last_generalized_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some((last.generalized_minimum_eigenvalue - previous.generalized_minimum_eigenvalue).abs())
    }

    pub fn last_reconstruction_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some((last.max_reconstruction_residual - previous.max_reconstruction_residual).abs())
    }

    pub fn last_family_delta_change(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some((last.generalized_family_delta - previous.generalized_family_delta).abs())
    }
}

#[derive(Debug)]
pub enum FiniteWeilSineTruncationError {
    EmptyModeSet,
    ZeroMode,
    ModesNotStrictlyIncreasing {
        previous: usize,
        next: usize,
    },
    EmptyParentDimensionSet,
    ParentDimensionsNotStrictlyIncreasing {
        previous: usize,
        next: usize,
    },
    ParentDimensionTooSmall {
        parent_dimension: usize,
        modes: usize,
    },
    Quadrature(QuadratureError),
    Generalized(FiniteWeilGeneralizedSpectrumError),
    CoefficientSubspace(FiniteWeilCoefficientSubspaceError),
    NonFiniteEvaluation {
        stage: &'static str,
        value: f64,
    },
}

impl fmt::Display for FiniteWeilSineTruncationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyModeSet => write!(f, "sine-enriched finite Weil family must contain a mode"),
            Self::ZeroMode => write!(
                f,
                "sine-enriched finite Weil modes must be positive integers"
            ),
            Self::ModesNotStrictlyIncreasing { previous, next } => write!(
                f,
                "sine modes must be strictly increasing: previous={previous}, next={next}"
            ),
            Self::EmptyParentDimensionSet => write!(
                f,
                "sine truncation audit requires at least one parent dimension"
            ),
            Self::ParentDimensionsNotStrictlyIncreasing { previous, next } => write!(
                f,
                "sine parent dimensions must be strictly increasing: previous={previous}, next={next}"
            ),
            Self::ParentDimensionTooSmall {
                parent_dimension,
                modes,
            } => write!(
                f,
                "sine parent dimension {parent_dimension} is smaller than subspace dimension {modes}"
            ),
            Self::Quadrature(error) => write!(f, "sine coefficient quadrature failed: {error:?}"),
            Self::Generalized(error) => write!(f, "sine parent finite Weil audit failed: {error}"),
            Self::CoefficientSubspace(error) => {
                write!(f, "sine coefficient-subspace audit failed: {error}")
            }
            Self::NonFiniteEvaluation { stage, value } => {
                write!(f, "non-finite sine truncation value at {stage}: {value}")
            }
        }
    }
}

impl std::error::Error for FiniteWeilSineTruncationError {}

impl From<QuadratureError> for FiniteWeilSineTruncationError {
    fn from(value: QuadratureError) -> Self {
        Self::Quadrature(value)
    }
}

impl From<FiniteWeilGeneralizedSpectrumError> for FiniteWeilSineTruncationError {
    fn from(value: FiniteWeilGeneralizedSpectrumError) -> Self {
        Self::Generalized(value)
    }
}

impl From<FiniteWeilCoefficientSubspaceError> for FiniteWeilSineTruncationError {
    fn from(value: FiniteWeilCoefficientSubspaceError) -> Self {
        Self::CoefficientSubspace(value)
    }
}

/// Project the enrichment `sin(mode*pi*t)` onto shifted Legendre degrees
/// `0..parent_dimension` by numerical Gauss--Legendre integration.
pub fn sine_legendre_coefficients(
    mode: usize,
    parent_dimension: usize,
    quadrature_order: usize,
) -> Result<Vec<f64>, FiniteWeilSineTruncationError> {
    if mode == 0 {
        return Err(FiniteWeilSineTruncationError::ZeroMode);
    }
    let quadrature = GaussLegendreUnit::new(quadrature_order)?;
    let mut coefficients = Vec::with_capacity(parent_dimension);
    for degree in 0..parent_dimension {
        let mut integral = 0.0_f64;
        for (&node, &weight) in quadrature.nodes().iter().zip(quadrature.weights().iter()) {
            integral +=
                weight * (mode as f64 * PI * node).sin() * shifted_legendre_value(degree, node);
        }
        let coefficient = (2 * degree + 1) as f64 * integral;
        checked_finite("sine Legendre coefficient", coefficient)?;
        coefficients.push(coefficient);
    }
    Ok(coefficients)
}

/// Re-evaluate the same sine-enriched finite subspace through increasing
/// Legendre parent dimensions while holding the declared Weil quadrature level fixed.
pub fn audit_finite_weil_sine_truncation(
    bump: CompactArchimedeanBump,
    modes: &SineModeSet,
    parent_dimensions: &[usize],
    coefficient_quadrature_order: usize,
    weil_level: WeilQuadratureLevel,
) -> Result<FiniteWeilSineTruncationAudit, FiniteWeilSineTruncationError> {
    validate_parent_dimensions(parent_dimensions, modes.dimension())?;

    let mut samples = Vec::with_capacity(parent_dimensions.len());
    for &parent_dimension in parent_dimensions {
        let parent = audit_finite_weil_generalized_spectrum(
            bump,
            parent_dimension,
            weil_level.correlation_order(),
            weil_level.archimedean_order(),
            weil_level.boundary_order(),
            weil_level.gram_order(),
        )?;

        let coefficients = modes
            .modes()
            .iter()
            .map(|&mode| {
                sine_legendre_coefficients(mode, parent_dimension, coefficient_quadrature_order)
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut max_reconstruction_residual = 0.0_f64;
        let mut max_coefficient_l1_norm = 0.0_f64;
        for (&mode, coefficient_vector) in modes.modes().iter().zip(coefficients.iter()) {
            let residual = sine_reconstruction_residual(mode, coefficient_vector);
            checked_finite("sine reconstruction residual", residual)?;
            max_reconstruction_residual = max_reconstruction_residual.max(residual);
            let l1 = coefficient_vector
                .iter()
                .map(|value| value.abs())
                .sum::<f64>();
            checked_finite("sine coefficient L1 norm", l1)?;
            max_coefficient_l1_norm = max_coefficient_l1_norm.max(l1);
        }

        let subspace = audit_finite_weil_coefficient_subspace(
            bump,
            &parent,
            &coefficients,
            weil_level.boundary_order(),
        )?;
        let generalized_minimum_eigenvalue = subspace.minimum_generalized_eigenvalue();
        let leading_legendre_generalized_minimum_eigenvalue =
            subspace.leading_legendre_generalized_minimum();
        let generalized_family_delta =
            generalized_minimum_eigenvalue - leading_legendre_generalized_minimum_eigenvalue;

        samples.push(SineTruncationSample {
            parent_dimension,
            max_reconstruction_residual,
            max_coefficient_l1_norm,
            raw_minimum_eigenvalue: subspace.minimum_raw_eigenvalue(),
            generalized_minimum_eigenvalue,
            leading_legendre_generalized_minimum_eigenvalue,
            generalized_family_delta,
            gram_condition_number: subspace.gram_condition_number(),
            leading_legendre_gram_condition_number: subspace
                .leading_legendre_gram_condition_number(),
            max_boundary_residual: subspace.max_boundary_residual(),
            max_pairing_asymmetry: subspace.parent_max_raw_pairing_asymmetry(),
            max_whitened_asymmetry: subspace.max_whitened_asymmetry(),
        });
    }

    let generalized_observed_minimum = samples
        .iter()
        .map(|sample| sample.generalized_minimum_eigenvalue)
        .fold(f64::INFINITY, f64::min);
    let generalized_observed_maximum = samples
        .iter()
        .map(|sample| sample.generalized_minimum_eigenvalue)
        .fold(f64::NEG_INFINITY, f64::max);
    let family_delta_observed_minimum = samples
        .iter()
        .map(|sample| sample.generalized_family_delta)
        .fold(f64::INFINITY, f64::min);
    let family_delta_observed_maximum = samples
        .iter()
        .map(|sample| sample.generalized_family_delta)
        .fold(f64::NEG_INFINITY, f64::max);

    Ok(FiniteWeilSineTruncationAudit {
        modes: modes.clone(),
        parent_dimensions: parent_dimensions.to_vec(),
        coefficient_quadrature_order,
        weil_level,
        samples,
        generalized_observed_minimum,
        generalized_observed_maximum,
        family_delta_observed_minimum,
        family_delta_observed_maximum,
    })
}

fn validate_parent_dimensions(
    parent_dimensions: &[usize],
    modes: usize,
) -> Result<(), FiniteWeilSineTruncationError> {
    if parent_dimensions.is_empty() {
        return Err(FiniteWeilSineTruncationError::EmptyParentDimensionSet);
    }
    for &parent_dimension in parent_dimensions {
        if parent_dimension < modes {
            return Err(FiniteWeilSineTruncationError::ParentDimensionTooSmall {
                parent_dimension,
                modes,
            });
        }
    }
    for pair in parent_dimensions.windows(2) {
        if pair[0] >= pair[1] {
            return Err(
                FiniteWeilSineTruncationError::ParentDimensionsNotStrictlyIncreasing {
                    previous: pair[0],
                    next: pair[1],
                },
            );
        }
    }
    Ok(())
}

fn sine_reconstruction_residual(mode: usize, coefficients: &[f64]) -> f64 {
    let mut maximum = 0.0_f64;
    for sample in 0..=RECONSTRUCTION_INTERVALS {
        let t = sample as f64 / RECONSTRUCTION_INTERVALS as f64;
        let direct = (mode as f64 * PI * t).sin();
        let expanded = coefficients
            .iter()
            .enumerate()
            .map(|(degree, &coefficient)| coefficient * shifted_legendre_value(degree, t))
            .sum::<f64>();
        maximum = maximum.max((direct - expanded).abs());
    }
    maximum
}

fn shifted_legendre_value(degree: usize, t: f64) -> f64 {
    let x = 2.0 * t - 1.0;
    if degree == 0 {
        return 1.0;
    }
    if degree == 1 {
        return x;
    }

    let mut p_nm2 = 1.0_f64;
    let mut p_nm1 = x;
    for n in 2..=degree {
        let p = ((2 * n - 1) as f64 * x * p_nm1 - (n - 1) as f64 * p_nm2) / n as f64;
        p_nm2 = p_nm1;
        p_nm1 = p;
    }
    p_nm1
}

fn checked_finite(stage: &'static str, value: f64) -> Result<(), FiniteWeilSineTruncationError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilSineTruncationError::NonFiniteEvaluation { stage, value })
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
    fn mode_set_rejects_zero_and_non_increasing_modes() {
        assert!(matches!(
            SineModeSet::new(vec![0, 1]),
            Err(FiniteWeilSineTruncationError::ZeroMode)
        ));
        assert!(matches!(
            SineModeSet::new(vec![1, 2, 2]),
            Err(FiniteWeilSineTruncationError::ModesNotStrictlyIncreasing { .. })
        ));
    }

    #[test]
    fn first_sine_mode_has_numerically_even_shifted_legendre_parity() {
        let coefficients = sine_legendre_coefficients(1, 6, 64).unwrap();
        assert_eq!(coefficients.len(), 6);
        for degree in [1_usize, 3, 5] {
            assert!(coefficients[degree].abs() <= 2.0e-13);
        }
    }

    #[test]
    fn truncation_sweep_records_finite_diagnostics_without_sign_assumption() {
        let modes = SineModeSet::new(vec![1, 2]).unwrap();
        let level = WeilQuadratureLevel::new(20, 20, 28, 28);
        let audit =
            audit_finite_weil_sine_truncation(bump(), &modes, &[3, 5], 48, level).unwrap();

        assert_eq!(audit.modes(), &modes);
        assert_eq!(audit.parent_dimensions(), &[3, 5]);
        assert_eq!(audit.weil_level(), level);
        assert_eq!(audit.samples().len(), 2);
        assert!(audit.generalized_observed_span().is_finite());
        assert!(audit.generalized_observed_span() >= 0.0);
        assert!(audit.family_delta_observed_span().is_finite());
        assert!(audit.family_delta_observed_span() >= 0.0);
        assert!(audit.last_generalized_delta().unwrap().is_finite());
        assert!(audit.last_reconstruction_delta().unwrap().is_finite());
        assert!(audit.last_family_delta_change().unwrap().is_finite());

        for sample in audit.samples() {
            assert!(sample.max_reconstruction_residual().is_finite());
            assert!(sample.max_coefficient_l1_norm().is_finite());
            assert!(sample.raw_minimum_eigenvalue().is_finite());
            assert!(sample.generalized_minimum_eigenvalue().is_finite());
            assert!(
                sample
                    .leading_legendre_generalized_minimum_eigenvalue()
                    .is_finite()
            );
            assert!(sample.generalized_family_delta().is_finite());
            assert!(sample.gram_condition_number().is_finite());
            assert!(sample.leading_legendre_gram_condition_number().is_finite());
            assert!(sample.max_boundary_residual().is_finite());
            assert!(sample.max_pairing_asymmetry().is_finite());
            assert!(sample.max_whitened_asymmetry().is_finite());
        }
    }
}
