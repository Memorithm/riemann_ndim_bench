//! Direct-Q cross-check for sine-enriched compact Weil generators.
//!
//! The truncation audit projects `sin(m*pi*t)` into the shifted-Legendre parent
//! before applying the already validated finite Weil machinery.  This module
//! adds a distinct local check: evaluate
//!
//! `g_m(rho) = bump(rho) sin(m*pi*t)`
//!
//! under `Q = -(rho d/drho)^2 + 1/4` by the product rule, then compare that
//! direct `Qg_m` with the linear combination of the existing Legendre images
//! `Q(bump P_j)` using the projection coefficients.
//!
//! The comparison is a numerical consistency audit.  Its sampled residuals are
//! not certified uniform bounds, do not bound the full Weil functional, and do
//! not establish convergence in the topology required by the Weil criterion or
//! RH.

use std::f64::consts::PI;
use std::fmt;

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::semilocal_compact_weil::CompactWeilTestFunction;
use crate::weil_boundary::{WeilBoundaryError, WeilBoundaryMoments, critical_boundary_moments};
use crate::weil_quadratic_matrix::CompactWeilBasisFunction;
use crate::weil_sine_truncation::{
    FiniteWeilSineTruncationError, SineModeSet, sine_legendre_coefficients,
};

/// One exact target generator `bump(rho) sin(mode*pi*t)` with direct product-rule Q.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactSineWeilFunction {
    bump: CompactArchimedeanBump,
    mode: usize,
}

impl CompactSineWeilFunction {
    pub fn new(
        bump: CompactArchimedeanBump,
        mode: usize,
    ) -> Result<Self, FiniteWeilDirectSineError> {
        if mode == 0 {
            return Err(FiniteWeilDirectSineError::ZeroMode);
        }
        Ok(Self { bump, mode })
    }

    #[inline]
    pub const fn bump(self) -> CompactArchimedeanBump {
        self.bump
    }

    #[inline]
    pub const fn mode(self) -> usize {
        self.mode
    }

    /// Evaluate the target generator before Q.
    pub fn generator_value(self, rho: f64) -> Result<f64, WeilBoundaryError> {
        let base = CompactWeilTestFunction::new(self.bump);
        let bump = base.generator_value(rho)?;
        if bump == 0.0 {
            return Ok(0.0);
        }
        let t = self.affine_coordinate(rho);
        Ok(bump * (self.omega() * t).sin())
    }

    /// Evaluate `Q(bump(rho) sin(mode*pi*t))` directly by the product rule.
    ///
    /// Writing `D=rho d/drho`,
    ///
    /// `Q(bs) = s Qb - 2(Db)(Ds) - b D^2 s`.
    ///
    /// `Qb` and `b` are taken from the existing compact Weil implementation;
    /// only the first logarithmic derivative of the bump and the analytic sine
    /// derivatives are assembled here.
    pub fn q_value(self, rho: f64) -> Result<f64, WeilBoundaryError> {
        if !rho.is_finite() || rho <= 0.0 {
            return Err(WeilBoundaryError::InvalidRho { rho });
        }
        if !self.bump.support().contains(rho) {
            return Ok(0.0);
        }
        self.q_value_inside(rho)
    }

    /// Numerically audit the two critical Mellin moments of the directly
    /// evaluated `Qg_m`.
    pub fn boundary_moments(
        self,
        quadrature_order: usize,
    ) -> Result<WeilBoundaryMoments, WeilBoundaryError> {
        critical_boundary_moments(self.bump.support(), quadrature_order, |rho| {
            self.q_value_inside(rho)
                .expect("boundary quadrature supplies finite positive interior rho")
        })
    }

    fn q_value_inside(self, rho: f64) -> Result<f64, WeilBoundaryError> {
        let base = CompactWeilTestFunction::new(self.bump);
        let bump = base.generator_value(rho)?;
        let q_bump = base.q_value(rho)?;
        if bump == 0.0 && q_bump == 0.0 {
            return Ok(0.0);
        }

        let support = self.bump.support();
        let width = support.upper() - support.lower();
        let t = (rho - support.lower()) / width;
        if !(0.0 < t && t < 1.0) {
            return Ok(0.0);
        }

        let d = t * (1.0 - t);
        let first_log_profile = (1.0 - 2.0 * t) / d.powi(2);
        let scale = rho / width;
        let d_bump = scale * bump * first_log_profile;

        let omega = self.omega();
        let sine = (omega * t).sin();
        let sine_t = omega * (omega * t).cos();
        let sine_tt = -omega * omega * sine;
        let d_sine = scale * sine_t;
        let d2_sine = scale * sine_t + scale * scale * sine_tt;

        let value = sine * q_bump - 2.0 * d_bump * d_sine - bump * d2_sine;
        debug_assert!(value.is_finite());
        Ok(value)
    }

    #[inline]
    fn affine_coordinate(self, rho: f64) -> f64 {
        let support = self.bump.support();
        (rho - support.lower()) / (support.upper() - support.lower())
    }

    #[inline]
    fn omega(self) -> f64 {
        self.mode as f64 * PI
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DirectSineQAuditConfig {
    coefficient_quadrature_order: usize,
    boundary_order: usize,
    comparison_intervals: usize,
}

impl DirectSineQAuditConfig {
    #[inline]
    pub const fn new(
        coefficient_quadrature_order: usize,
        boundary_order: usize,
        comparison_intervals: usize,
    ) -> Self {
        Self {
            coefficient_quadrature_order,
            boundary_order,
            comparison_intervals,
        }
    }

    #[inline]
    pub const fn coefficient_quadrature_order(self) -> usize {
        self.coefficient_quadrature_order
    }

    #[inline]
    pub const fn boundary_order(self) -> usize {
        self.boundary_order
    }

    #[inline]
    pub const fn comparison_intervals(self) -> usize {
        self.comparison_intervals
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DirectSineQComparisonSample {
    parent_dimension: usize,
    max_generator_residual: f64,
    max_direct_q_amplitude: f64,
    max_projected_q_amplitude: f64,
    max_q_absolute_residual: f64,
    normalized_q_residual: f64,
    max_direct_boundary_residual: f64,
    max_projected_boundary_residual: f64,
}

impl DirectSineQComparisonSample {
    #[inline]
    pub const fn parent_dimension(self) -> usize {
        self.parent_dimension
    }

    #[inline]
    pub const fn max_generator_residual(self) -> f64 {
        self.max_generator_residual
    }

    #[inline]
    pub const fn max_direct_q_amplitude(self) -> f64 {
        self.max_direct_q_amplitude
    }

    #[inline]
    pub const fn max_projected_q_amplitude(self) -> f64 {
        self.max_projected_q_amplitude
    }

    #[inline]
    pub const fn max_q_absolute_residual(self) -> f64 {
        self.max_q_absolute_residual
    }

    #[inline]
    pub const fn normalized_q_residual(self) -> f64 {
        self.normalized_q_residual
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

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilDirectSineQAudit {
    modes: SineModeSet,
    parent_dimensions: Vec<usize>,
    config: DirectSineQAuditConfig,
    samples: Vec<DirectSineQComparisonSample>,
}

impl FiniteWeilDirectSineQAudit {
    #[inline]
    pub fn modes(&self) -> &SineModeSet {
        &self.modes
    }

    #[inline]
    pub fn parent_dimensions(&self) -> &[usize] {
        &self.parent_dimensions
    }

    #[inline]
    pub const fn config(&self) -> DirectSineQAuditConfig {
        self.config
    }

    #[inline]
    pub fn samples(&self) -> &[DirectSineQComparisonSample] {
        &self.samples
    }

    pub fn last_q_residual_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some((last.max_q_absolute_residual - previous.max_q_absolute_residual).abs())
    }

    pub fn last_normalized_q_residual_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some((last.normalized_q_residual - previous.normalized_q_residual).abs())
    }
}

#[derive(Debug)]
pub enum FiniteWeilDirectSineError {
    ZeroMode,
    EmptyParentDimensionSet,
    ZeroParentDimension,
    ParentDimensionsNotStrictlyIncreasing { previous: usize, next: usize },
    ZeroComparisonIntervals,
    SineProjection(FiniteWeilSineTruncationError),
    Boundary(WeilBoundaryError),
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilDirectSineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroMode => write!(f, "direct sine Weil mode must be positive"),
            Self::EmptyParentDimensionSet => {
                write!(
                    f,
                    "direct sine Q audit requires at least one parent dimension"
                )
            }
            Self::ZeroParentDimension => {
                write!(f, "direct sine Q parent dimension must be positive")
            }
            Self::ParentDimensionsNotStrictlyIncreasing { previous, next } => write!(
                f,
                "direct sine Q parent dimensions must be strictly increasing: previous={previous}, next={next}"
            ),
            Self::ZeroComparisonIntervals => {
                write!(
                    f,
                    "direct sine Q audit requires at least one comparison interval"
                )
            }
            Self::SineProjection(error) => write!(f, "sine Legendre projection failed: {error}"),
            Self::Boundary(error) => write!(f, "direct sine Q boundary evaluation failed: {error}"),
            Self::NonFiniteEvaluation { stage, value } => {
                write!(f, "non-finite direct sine Q diagnostic at {stage}: {value}")
            }
        }
    }
}

impl std::error::Error for FiniteWeilDirectSineError {}

impl From<FiniteWeilSineTruncationError> for FiniteWeilDirectSineError {
    fn from(value: FiniteWeilSineTruncationError) -> Self {
        Self::SineProjection(value)
    }
}

impl From<WeilBoundaryError> for FiniteWeilDirectSineError {
    fn from(value: WeilBoundaryError) -> Self {
        Self::Boundary(value)
    }
}

pub fn audit_finite_weil_direct_sine_q(
    bump: CompactArchimedeanBump,
    modes: &SineModeSet,
    parent_dimensions: &[usize],
    config: DirectSineQAuditConfig,
) -> Result<FiniteWeilDirectSineQAudit, FiniteWeilDirectSineError> {
    validate_parent_dimensions(parent_dimensions)?;
    if config.comparison_intervals() == 0 {
        return Err(FiniteWeilDirectSineError::ZeroComparisonIntervals);
    }

    let direct_functions = modes
        .modes()
        .iter()
        .map(|&mode| CompactSineWeilFunction::new(bump, mode))
        .collect::<Result<Vec<_>, _>>()?;

    let mut direct_boundary_residuals = Vec::with_capacity(direct_functions.len());
    for &direct in &direct_functions {
        let moments = direct.boundary_moments(config.boundary_order())?;
        direct_boundary_residuals.push(moments.plus_half.abs().max(moments.minus_half.abs()));
    }

    let support = bump.support();
    let width = support.upper() - support.lower();
    let mut samples = Vec::with_capacity(parent_dimensions.len());

    for &parent_dimension in parent_dimensions {
        let coefficients = modes
            .modes()
            .iter()
            .map(|&mode| {
                sine_legendre_coefficients(
                    mode,
                    parent_dimension,
                    config.coefficient_quadrature_order(),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;

        let parent_boundary_moments = (0..parent_dimension)
            .map(|degree| {
                CompactWeilBasisFunction::new(bump, degree)
                    .boundary_moments(config.boundary_order())
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut max_generator_residual = 0.0_f64;
        let mut max_direct_q_amplitude = 0.0_f64;
        let mut max_projected_q_amplitude = 0.0_f64;
        let mut max_q_absolute_residual = 0.0_f64;
        let mut max_projected_boundary_residual = 0.0_f64;

        for (mode_index, (&mode, coefficient_vector)) in
            modes.modes().iter().zip(coefficients.iter()).enumerate()
        {
            let direct = direct_functions[mode_index];

            let projected_plus = coefficient_vector
                .iter()
                .zip(parent_boundary_moments.iter())
                .map(|(&coefficient, moments)| coefficient * moments.plus_half)
                .sum::<f64>();
            let projected_minus = coefficient_vector
                .iter()
                .zip(parent_boundary_moments.iter())
                .map(|(&coefficient, moments)| coefficient * moments.minus_half)
                .sum::<f64>();
            max_projected_boundary_residual = max_projected_boundary_residual
                .max(projected_plus.abs().max(projected_minus.abs()));

            for sample_index in 0..=config.comparison_intervals() {
                let t = sample_index as f64 / config.comparison_intervals() as f64;
                let rho = support.lower() + width * t;
                let direct_generator = direct.generator_value(rho)?;
                let direct_q = direct.q_value(rho)?;

                let mut projected_profile = 0.0_f64;
                let mut projected_q = 0.0_f64;
                for (degree, &coefficient) in coefficient_vector.iter().enumerate() {
                    projected_profile += coefficient * shifted_legendre_value(degree, t);
                    projected_q +=
                        coefficient * CompactWeilBasisFunction::new(bump, degree).value(rho)?;
                }
                let bump_value = CompactWeilTestFunction::new(bump).generator_value(rho)?;
                let projected_generator = bump_value * projected_profile;

                max_generator_residual =
                    max_generator_residual.max((direct_generator - projected_generator).abs());
                max_direct_q_amplitude = max_direct_q_amplitude.max(direct_q.abs());
                max_projected_q_amplitude = max_projected_q_amplitude.max(projected_q.abs());
                max_q_absolute_residual =
                    max_q_absolute_residual.max((direct_q - projected_q).abs());
            }

            let direct_boundary = direct_boundary_residuals[mode_index];
            checked_finite("direct boundary residual", direct_boundary)?;
            checked_finite("mode", mode as f64)?;
        }

        let max_direct_boundary_residual = direct_boundary_residuals
            .iter()
            .copied()
            .fold(0.0_f64, f64::max);
        let normalized_q_residual = if max_direct_q_amplitude > 0.0 {
            max_q_absolute_residual / max_direct_q_amplitude
        } else {
            max_q_absolute_residual
        };

        for (stage, value) in [
            ("generator residual", max_generator_residual),
            ("direct Q amplitude", max_direct_q_amplitude),
            ("projected Q amplitude", max_projected_q_amplitude),
            ("Q absolute residual", max_q_absolute_residual),
            ("normalized Q residual", normalized_q_residual),
            ("direct boundary residual", max_direct_boundary_residual),
            (
                "projected boundary residual",
                max_projected_boundary_residual,
            ),
        ] {
            checked_finite(stage, value)?;
        }

        samples.push(DirectSineQComparisonSample {
            parent_dimension,
            max_generator_residual,
            max_direct_q_amplitude,
            max_projected_q_amplitude,
            max_q_absolute_residual,
            normalized_q_residual,
            max_direct_boundary_residual,
            max_projected_boundary_residual,
        });
    }

    Ok(FiniteWeilDirectSineQAudit {
        modes: modes.clone(),
        parent_dimensions: parent_dimensions.to_vec(),
        config,
        samples,
    })
}

fn validate_parent_dimensions(
    parent_dimensions: &[usize],
) -> Result<(), FiniteWeilDirectSineError> {
    if parent_dimensions.is_empty() {
        return Err(FiniteWeilDirectSineError::EmptyParentDimensionSet);
    }
    if parent_dimensions.contains(&0) {
        return Err(FiniteWeilDirectSineError::ZeroParentDimension);
    }
    for pair in parent_dimensions.windows(2) {
        if pair[0] >= pair[1] {
            return Err(
                FiniteWeilDirectSineError::ParentDimensionsNotStrictlyIncreasing {
                    previous: pair[0],
                    next: pair[1],
                },
            );
        }
    }
    Ok(())
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

fn checked_finite(stage: &'static str, value: f64) -> Result<(), FiniteWeilDirectSineError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilDirectSineError::NonFiniteEvaluation { stage, value })
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
    fn direct_sine_q_preserves_compact_support() {
        let function = CompactSineWeilFunction::new(bump(), 1).unwrap();
        let support = bump().support();
        assert_eq!(function.q_value(support.lower()).unwrap(), 0.0);
        assert_eq!(function.q_value(support.upper()).unwrap(), 0.0);
        assert_eq!(function.q_value(0.5 * support.lower()).unwrap(), 0.0);
        assert_eq!(function.q_value(2.0 * support.upper()).unwrap(), 0.0);
    }

    #[test]
    fn direct_sine_q_boundary_moments_are_finite() {
        let function = CompactSineWeilFunction::new(bump(), 1).unwrap();
        let moments = function.boundary_moments(64).unwrap();
        assert!(moments.plus_half.is_finite());
        assert!(moments.minus_half.is_finite());
    }

    #[test]
    fn direct_projection_audit_records_q_residuals_without_sign_claim() {
        let modes = SineModeSet::new(vec![1, 2]).unwrap();
        let config = DirectSineQAuditConfig::new(48, 40, 64);
        let audit = audit_finite_weil_direct_sine_q(bump(), &modes, &[3, 5], config).unwrap();
        assert_eq!(audit.modes(), &modes);
        assert_eq!(audit.parent_dimensions(), &[3, 5]);
        assert_eq!(audit.config(), config);
        assert_eq!(audit.samples().len(), 2);
        assert!(audit.last_q_residual_delta().unwrap().is_finite());
        assert!(
            audit
                .last_normalized_q_residual_delta()
                .unwrap()
                .is_finite()
        );
        for sample in audit.samples() {
            assert!(sample.max_generator_residual().is_finite());
            assert!(sample.max_direct_q_amplitude().is_finite());
            assert!(sample.max_projected_q_amplitude().is_finite());
            assert!(sample.max_q_absolute_residual().is_finite());
            assert!(sample.normalized_q_residual().is_finite());
            assert!(sample.max_direct_boundary_residual().is_finite());
            assert!(sample.max_projected_boundary_residual().is_finite());
        }
    }
}
