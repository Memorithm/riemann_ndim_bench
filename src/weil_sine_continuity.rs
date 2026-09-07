//! Empirical L2 continuity probe for direct-vs-projected sine Weil pairings.
//!
//! For each sine mode this module compares the direct compact function
//! `h = Q(bump * sin(m*pi*t))` with its finite shifted-Legendre projection `h_N`
//! in the multiplicative norm
//!
//! `||f||_2^2 = integral |f(rho)|^2 d rho / rho`.
//!
//! It then compares the full finite Riemann--Weil pairing residual with the
//! symmetric perturbation scale
//!
//! `e_i max(||h_j||, ||h_j,N||) + e_j max(||h_i||, ||h_i,N||)`,
//!
//! where `e_k = ||h_k-h_k,N||_2`.
//!
//! The resulting quotient is an observed finite-dimensional diagnostic only.
//! It is not a proved operator norm, continuity theorem, certified error bound,
//! density/completeness result, complete-space Weil positivity statement,
//! Conjecture 4.1, or RH.

use std::fmt;

use crate::quadrature::{GaussLegendreUnit, QuadratureError};
use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_boundary::{WeilBoundaryError, WeilBoundaryMoments};
use crate::weil_compact_pairing::{
    CompactWeilEvaluand, CompactWeilPairingConfig, CompactWeilPairingError,
    FiniteCompactWeilPairingAudit, audit_compact_weil_pairing,
};
use crate::weil_quadratic_matrix::CompactWeilBasisFunction;
use crate::weil_sine_direct_q::{CompactSineWeilFunction, FiniteWeilDirectSineError};
use crate::weil_sine_truncation::{
    FiniteWeilSineTruncationError, SineModeSet, sine_legendre_coefficients,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SineContinuityProbeConfig {
    coefficient_quadrature_order: usize,
    norm_quadrature_order: usize,
    pairing: CompactWeilPairingConfig,
}

impl SineContinuityProbeConfig {
    #[inline]
    pub const fn new(
        coefficient_quadrature_order: usize,
        norm_quadrature_order: usize,
        pairing: CompactWeilPairingConfig,
    ) -> Self {
        Self {
            coefficient_quadrature_order,
            norm_quadrature_order,
            pairing,
        }
    }

    #[inline]
    pub const fn coefficient_quadrature_order(self) -> usize {
        self.coefficient_quadrature_order
    }

    #[inline]
    pub const fn norm_quadrature_order(self) -> usize {
        self.norm_quadrature_order
    }

    #[inline]
    pub const fn pairing(self) -> CompactWeilPairingConfig {
        self.pairing
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SineModeL2Comparison {
    mode: usize,
    direct_l2_norm: f64,
    projected_l2_norm: f64,
    l2_error: f64,
    relative_l2_error: f64,
    observed_max_absolute_error: f64,
}

impl SineModeL2Comparison {
    #[inline]
    pub const fn mode(self) -> usize {
        self.mode
    }

    #[inline]
    pub const fn direct_l2_norm(self) -> f64 {
        self.direct_l2_norm
    }

    #[inline]
    pub const fn projected_l2_norm(self) -> f64 {
        self.projected_l2_norm
    }

    #[inline]
    pub const fn l2_error(self) -> f64 {
        self.l2_error
    }

    #[inline]
    pub const fn relative_l2_error(self) -> f64 {
        self.relative_l2_error
    }

    #[inline]
    pub const fn observed_max_absolute_error(self) -> f64 {
        self.observed_max_absolute_error
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SinePairingContinuityEntry {
    left_mode: usize,
    right_mode: usize,
    direct_pairing: f64,
    projected_pairing: f64,
    pairing_residual: f64,
    pole_residual: f64,
    archimedean_residual: f64,
    prime_residual: f64,
    l2_perturbation_scale: f64,
    observed_l2_continuity_quotient: f64,
    direct_pairing_asymmetry: f64,
    projected_pairing_asymmetry: f64,
}

impl SinePairingContinuityEntry {
    #[inline]
    pub const fn left_mode(self) -> usize {
        self.left_mode
    }

    #[inline]
    pub const fn right_mode(self) -> usize {
        self.right_mode
    }

    #[inline]
    pub const fn direct_pairing(self) -> f64 {
        self.direct_pairing
    }

    #[inline]
    pub const fn projected_pairing(self) -> f64 {
        self.projected_pairing
    }

    #[inline]
    pub const fn pairing_residual(self) -> f64 {
        self.pairing_residual
    }

    #[inline]
    pub const fn pole_residual(self) -> f64 {
        self.pole_residual
    }

    #[inline]
    pub const fn archimedean_residual(self) -> f64 {
        self.archimedean_residual
    }

    #[inline]
    pub const fn prime_residual(self) -> f64 {
        self.prime_residual
    }

    #[inline]
    pub const fn l2_perturbation_scale(self) -> f64 {
        self.l2_perturbation_scale
    }

    #[inline]
    pub const fn observed_l2_continuity_quotient(self) -> f64 {
        self.observed_l2_continuity_quotient
    }

    #[inline]
    pub const fn direct_pairing_asymmetry(self) -> f64 {
        self.direct_pairing_asymmetry
    }

    #[inline]
    pub const fn projected_pairing_asymmetry(self) -> f64 {
        self.projected_pairing_asymmetry
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SineContinuityProbeSample {
    parent_dimension: usize,
    mode_norms: Vec<SineModeL2Comparison>,
    pairings: Vec<SinePairingContinuityEntry>,
    max_l2_error: f64,
    max_relative_l2_error: f64,
    max_pairing_residual: f64,
    max_observed_l2_continuity_quotient: f64,
    max_pole_residual: f64,
    max_archimedean_residual: f64,
    max_prime_residual: f64,
    max_pairing_asymmetry: f64,
}

impl SineContinuityProbeSample {
    #[inline]
    pub const fn parent_dimension(&self) -> usize {
        self.parent_dimension
    }

    #[inline]
    pub fn mode_norms(&self) -> &[SineModeL2Comparison] {
        &self.mode_norms
    }

    #[inline]
    pub fn pairings(&self) -> &[SinePairingContinuityEntry] {
        &self.pairings
    }

    #[inline]
    pub const fn max_l2_error(&self) -> f64 {
        self.max_l2_error
    }

    #[inline]
    pub const fn max_relative_l2_error(&self) -> f64 {
        self.max_relative_l2_error
    }

    #[inline]
    pub const fn max_pairing_residual(&self) -> f64 {
        self.max_pairing_residual
    }

    #[inline]
    pub const fn max_observed_l2_continuity_quotient(&self) -> f64 {
        self.max_observed_l2_continuity_quotient
    }

    #[inline]
    pub const fn max_pole_residual(&self) -> f64 {
        self.max_pole_residual
    }

    #[inline]
    pub const fn max_archimedean_residual(&self) -> f64 {
        self.max_archimedean_residual
    }

    #[inline]
    pub const fn max_prime_residual(&self) -> f64 {
        self.max_prime_residual
    }

    #[inline]
    pub const fn max_pairing_asymmetry(&self) -> f64 {
        self.max_pairing_asymmetry
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilSineContinuityProbe {
    modes: SineModeSet,
    parent_dimensions: Vec<usize>,
    config: SineContinuityProbeConfig,
    samples: Vec<SineContinuityProbeSample>,
}

impl FiniteWeilSineContinuityProbe {
    #[inline]
    pub fn modes(&self) -> &SineModeSet {
        &self.modes
    }

    #[inline]
    pub fn parent_dimensions(&self) -> &[usize] {
        &self.parent_dimensions
    }

    #[inline]
    pub const fn config(&self) -> SineContinuityProbeConfig {
        self.config
    }

    #[inline]
    pub fn samples(&self) -> &[SineContinuityProbeSample] {
        &self.samples
    }

    pub fn last_l2_error_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some((last.max_l2_error - previous.max_l2_error).abs())
    }

    pub fn last_pairing_residual_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some((last.max_pairing_residual - previous.max_pairing_residual).abs())
    }

    pub fn last_continuity_quotient_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some(
            (last.max_observed_l2_continuity_quotient
                - previous.max_observed_l2_continuity_quotient)
                .abs(),
        )
    }
}

#[derive(Debug)]
pub enum FiniteWeilSineContinuityError {
    EmptyParentDimensionSet,
    ZeroParentDimension,
    ParentDimensionsNotStrictlyIncreasing { previous: usize, next: usize },
    Quadrature(QuadratureError),
    DirectSine(FiniteWeilDirectSineError),
    SineProjection(FiniteWeilSineTruncationError),
    Pairing(CompactWeilPairingError),
    Boundary(WeilBoundaryError),
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilSineContinuityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyParentDimensionSet => {
                write!(f, "sine continuity probe requires a parent dimension")
            }
            Self::ZeroParentDimension => {
                write!(f, "sine continuity parent dimension must be positive")
            }
            Self::ParentDimensionsNotStrictlyIncreasing { previous, next } => write!(
                f,
                "sine continuity parent dimensions must be strictly increasing: previous={previous}, next={next}"
            ),
            Self::Quadrature(error) => write!(f, "sine continuity quadrature failed: {error:?}"),
            Self::DirectSine(error) => write!(f, "direct sine construction failed: {error}"),
            Self::SineProjection(error) => write!(f, "sine projection failed: {error}"),
            Self::Pairing(error) => write!(f, "compact Weil pairing failed: {error}"),
            Self::Boundary(error) => write!(f, "compact Weil norm evaluation failed: {error}"),
            Self::NonFiniteEvaluation { stage, value } => {
                write!(
                    f,
                    "non-finite sine continuity diagnostic at {stage}: {value}"
                )
            }
        }
    }
}

impl std::error::Error for FiniteWeilSineContinuityError {}

impl From<QuadratureError> for FiniteWeilSineContinuityError {
    fn from(value: QuadratureError) -> Self {
        Self::Quadrature(value)
    }
}

impl From<FiniteWeilDirectSineError> for FiniteWeilSineContinuityError {
    fn from(value: FiniteWeilDirectSineError) -> Self {
        Self::DirectSine(value)
    }
}

impl From<FiniteWeilSineTruncationError> for FiniteWeilSineContinuityError {
    fn from(value: FiniteWeilSineTruncationError) -> Self {
        Self::SineProjection(value)
    }
}

impl From<CompactWeilPairingError> for FiniteWeilSineContinuityError {
    fn from(value: CompactWeilPairingError) -> Self {
        Self::Pairing(value)
    }
}

impl From<WeilBoundaryError> for FiniteWeilSineContinuityError {
    fn from(value: WeilBoundaryError) -> Self {
        Self::Boundary(value)
    }
}

#[derive(Clone, Debug)]
struct ProjectedSineWeilFunction {
    bump: CompactArchimedeanBump,
    coefficients: Vec<f64>,
}

impl CompactWeilEvaluand for ProjectedSineWeilFunction {
    #[inline]
    fn bump(&self) -> CompactArchimedeanBump {
        self.bump
    }

    fn value(&self, rho: f64) -> Result<f64, WeilBoundaryError> {
        let mut total = 0.0_f64;
        for (degree, &coefficient) in self.coefficients.iter().enumerate() {
            total += coefficient * CompactWeilBasisFunction::new(self.bump, degree).value(rho)?;
        }
        Ok(total)
    }

    fn boundary_moments(
        &self,
        quadrature_order: usize,
    ) -> Result<WeilBoundaryMoments, WeilBoundaryError> {
        let mut plus_half = 0.0_f64;
        let mut minus_half = 0.0_f64;
        for (degree, &coefficient) in self.coefficients.iter().enumerate() {
            let moments = CompactWeilBasisFunction::new(self.bump, degree)
                .boundary_moments(quadrature_order)?;
            plus_half += coefficient * moments.plus_half;
            minus_half += coefficient * moments.minus_half;
        }
        Ok(WeilBoundaryMoments {
            plus_half,
            minus_half,
        })
    }
}

#[derive(Clone, Copy)]
struct SymmetricPairing {
    value: f64,
    pole_term: f64,
    archimedean_term: f64,
    prime_total: f64,
    asymmetry: f64,
}

pub fn audit_finite_weil_sine_continuity(
    bump: CompactArchimedeanBump,
    modes: &SineModeSet,
    parent_dimensions: &[usize],
    config: SineContinuityProbeConfig,
) -> Result<FiniteWeilSineContinuityProbe, FiniteWeilSineContinuityError> {
    validate_parent_dimensions(parent_dimensions)?;

    let direct_functions = modes
        .modes()
        .iter()
        .map(|&mode| CompactSineWeilFunction::new(bump, mode))
        .collect::<Result<Vec<_>, _>>()?;
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
        let projected_functions = coefficients
            .iter()
            .cloned()
            .map(|coefficients| ProjectedSineWeilFunction { bump, coefficients })
            .collect::<Vec<_>>();

        let mut mode_norms = Vec::with_capacity(modes.dimension());
        for (mode_index, &mode) in modes.modes().iter().enumerate() {
            mode_norms.push(audit_l2_comparison(
                mode,
                &direct_functions[mode_index],
                &projected_functions[mode_index],
                config.norm_quadrature_order(),
            )?);
        }

        let mut pairings = Vec::with_capacity(modes.dimension() * (modes.dimension() + 1) / 2);
        for i in 0..modes.dimension() {
            for j in i..modes.dimension() {
                let direct = symmetric_pairing(
                    &direct_functions[i],
                    &direct_functions[j],
                    i == j,
                    config.pairing(),
                )?;
                let projected = symmetric_pairing(
                    &projected_functions[i],
                    &projected_functions[j],
                    i == j,
                    config.pairing(),
                )?;
                let pairing_residual = (direct.value - projected.value).abs();
                let left_error = mode_norms[i].l2_error;
                let right_error = mode_norms[j].l2_error;
                let left_scale = mode_norms[i]
                    .direct_l2_norm
                    .max(mode_norms[i].projected_l2_norm);
                let right_scale = mode_norms[j]
                    .direct_l2_norm
                    .max(mode_norms[j].projected_l2_norm);
                let l2_perturbation_scale = left_error * right_scale + right_error * left_scale;
                let observed_l2_continuity_quotient = if l2_perturbation_scale > 0.0 {
                    pairing_residual / l2_perturbation_scale
                } else {
                    0.0
                };
                let entry = SinePairingContinuityEntry {
                    left_mode: modes.modes()[i],
                    right_mode: modes.modes()[j],
                    direct_pairing: direct.value,
                    projected_pairing: projected.value,
                    pairing_residual,
                    pole_residual: (direct.pole_term - projected.pole_term).abs(),
                    archimedean_residual: (direct.archimedean_term - projected.archimedean_term)
                        .abs(),
                    prime_residual: (direct.prime_total - projected.prime_total).abs(),
                    l2_perturbation_scale,
                    observed_l2_continuity_quotient,
                    direct_pairing_asymmetry: direct.asymmetry,
                    projected_pairing_asymmetry: projected.asymmetry,
                };
                check_entry(entry)?;
                pairings.push(entry);
            }
        }

        let max_l2_error = mode_norms
            .iter()
            .map(|sample| sample.l2_error)
            .fold(0.0_f64, f64::max);
        let max_relative_l2_error = mode_norms
            .iter()
            .map(|sample| sample.relative_l2_error)
            .fold(0.0_f64, f64::max);
        let max_pairing_residual = pairings
            .iter()
            .map(|entry| entry.pairing_residual)
            .fold(0.0_f64, f64::max);
        let max_observed_l2_continuity_quotient = pairings
            .iter()
            .map(|entry| entry.observed_l2_continuity_quotient)
            .fold(0.0_f64, f64::max);
        let max_pole_residual = pairings
            .iter()
            .map(|entry| entry.pole_residual)
            .fold(0.0_f64, f64::max);
        let max_archimedean_residual = pairings
            .iter()
            .map(|entry| entry.archimedean_residual)
            .fold(0.0_f64, f64::max);
        let max_prime_residual = pairings
            .iter()
            .map(|entry| entry.prime_residual)
            .fold(0.0_f64, f64::max);
        let max_pairing_asymmetry = pairings
            .iter()
            .map(|entry| {
                entry
                    .direct_pairing_asymmetry
                    .max(entry.projected_pairing_asymmetry)
            })
            .fold(0.0_f64, f64::max);

        for (stage, value) in [
            ("maximum L2 error", max_l2_error),
            ("maximum relative L2 error", max_relative_l2_error),
            ("maximum pairing residual", max_pairing_residual),
            (
                "maximum observed L2 continuity quotient",
                max_observed_l2_continuity_quotient,
            ),
            ("maximum pole residual", max_pole_residual),
            ("maximum archimedean residual", max_archimedean_residual),
            ("maximum prime residual", max_prime_residual),
            ("maximum pairing asymmetry", max_pairing_asymmetry),
        ] {
            checked_finite(stage, value)?;
        }

        samples.push(SineContinuityProbeSample {
            parent_dimension,
            mode_norms,
            pairings,
            max_l2_error,
            max_relative_l2_error,
            max_pairing_residual,
            max_observed_l2_continuity_quotient,
            max_pole_residual,
            max_archimedean_residual,
            max_prime_residual,
            max_pairing_asymmetry,
        });
    }

    Ok(FiniteWeilSineContinuityProbe {
        modes: modes.clone(),
        parent_dimensions: parent_dimensions.to_vec(),
        config,
        samples,
    })
}

fn audit_l2_comparison<L, R>(
    mode: usize,
    direct: &L,
    projected: &R,
    quadrature_order: usize,
) -> Result<SineModeL2Comparison, FiniteWeilSineContinuityError>
where
    L: CompactWeilEvaluand + ?Sized,
    R: CompactWeilEvaluand + ?Sized,
{
    let support = direct.bump().support();
    let log_lower = support.log_lower();
    let log_upper = support.log_upper();
    let log_span = log_upper - log_lower;
    let quadrature = GaussLegendreUnit::new(quadrature_order)?;
    let mut direct_square = 0.0_f64;
    let mut projected_square = 0.0_f64;
    let mut error_square = 0.0_f64;
    let mut observed_max_absolute_error = 0.0_f64;

    for (&node, &weight) in quadrature.nodes().iter().zip(quadrature.weights().iter()) {
        let rho = (log_lower + log_span * node).exp();
        let direct_value = direct.value(rho)?;
        let projected_value = projected.value(rho)?;
        let error = direct_value - projected_value;
        direct_square += weight * direct_value * direct_value;
        projected_square += weight * projected_value * projected_value;
        error_square += weight * error * error;
        observed_max_absolute_error = observed_max_absolute_error.max(error.abs());
    }

    let direct_l2_norm = (log_span * direct_square).sqrt();
    let projected_l2_norm = (log_span * projected_square).sqrt();
    let l2_error = (log_span * error_square).sqrt();
    let reference = direct_l2_norm.max(projected_l2_norm);
    let relative_l2_error = if reference > 0.0 {
        l2_error / reference
    } else {
        0.0
    };
    for (stage, value) in [
        ("direct L2 norm", direct_l2_norm),
        ("projected L2 norm", projected_l2_norm),
        ("L2 error", l2_error),
        ("relative L2 error", relative_l2_error),
        (
            "observed maximum absolute error",
            observed_max_absolute_error,
        ),
    ] {
        checked_finite(stage, value)?;
    }
    Ok(SineModeL2Comparison {
        mode,
        direct_l2_norm,
        projected_l2_norm,
        l2_error,
        relative_l2_error,
        observed_max_absolute_error,
    })
}

fn symmetric_pairing<L, R>(
    left: &L,
    right: &R,
    diagonal: bool,
    config: CompactWeilPairingConfig,
) -> Result<SymmetricPairing, FiniteWeilSineContinuityError>
where
    L: CompactWeilEvaluand + ?Sized,
    R: CompactWeilEvaluand + ?Sized,
{
    let forward = audit_compact_weil_pairing(left, right, config)?;
    if diagonal {
        return Ok(from_single_pairing(forward));
    }
    let reverse = audit_compact_weil_pairing(right, left, config)?;
    Ok(SymmetricPairing {
        value: 0.5 * (forward.value() + reverse.value()),
        pole_term: 0.5 * (forward.pole_term() + reverse.pole_term()),
        archimedean_term: 0.5 * (forward.archimedean_term() + reverse.archimedean_term()),
        prime_total: 0.5 * (forward.prime_total() + reverse.prime_total()),
        asymmetry: (forward.value() - reverse.value()).abs(),
    })
}

fn from_single_pairing(pairing: FiniteCompactWeilPairingAudit) -> SymmetricPairing {
    SymmetricPairing {
        value: pairing.value(),
        pole_term: pairing.pole_term(),
        archimedean_term: pairing.archimedean_term(),
        prime_total: pairing.prime_total(),
        asymmetry: 0.0,
    }
}

fn validate_parent_dimensions(
    parent_dimensions: &[usize],
) -> Result<(), FiniteWeilSineContinuityError> {
    if parent_dimensions.is_empty() {
        return Err(FiniteWeilSineContinuityError::EmptyParentDimensionSet);
    }
    if parent_dimensions.contains(&0) {
        return Err(FiniteWeilSineContinuityError::ZeroParentDimension);
    }
    for pair in parent_dimensions.windows(2) {
        if pair[0] >= pair[1] {
            return Err(
                FiniteWeilSineContinuityError::ParentDimensionsNotStrictlyIncreasing {
                    previous: pair[0],
                    next: pair[1],
                },
            );
        }
    }
    Ok(())
}

fn check_entry(entry: SinePairingContinuityEntry) -> Result<(), FiniteWeilSineContinuityError> {
    for (stage, value) in [
        ("direct pairing", entry.direct_pairing),
        ("projected pairing", entry.projected_pairing),
        ("pairing residual", entry.pairing_residual),
        ("pole residual", entry.pole_residual),
        ("archimedean residual", entry.archimedean_residual),
        ("prime residual", entry.prime_residual),
        ("L2 perturbation scale", entry.l2_perturbation_scale),
        (
            "observed L2 continuity quotient",
            entry.observed_l2_continuity_quotient,
        ),
        ("direct pairing asymmetry", entry.direct_pairing_asymmetry),
        (
            "projected pairing asymmetry",
            entry.projected_pairing_asymmetry,
        ),
    ] {
        checked_finite(stage, value)?;
    }
    Ok(())
}

fn checked_finite(stage: &'static str, value: f64) -> Result<(), FiniteWeilSineContinuityError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilSineContinuityError::NonFiniteEvaluation { stage, value })
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
    fn continuity_probe_records_l2_and_pairing_residuals_without_theorem_claim() {
        let modes = SineModeSet::new(vec![1, 2]).unwrap();
        let config =
            SineContinuityProbeConfig::new(40, 32, CompactWeilPairingConfig::new(20, 20, 28));
        let audit = audit_finite_weil_sine_continuity(bump(), &modes, &[3, 5], config).unwrap();

        assert_eq!(audit.modes(), &modes);
        assert_eq!(audit.parent_dimensions(), &[3, 5]);
        assert_eq!(audit.samples().len(), 2);
        assert!(audit.last_l2_error_delta().unwrap().is_finite());
        assert!(audit.last_pairing_residual_delta().unwrap().is_finite());
        assert!(audit.last_continuity_quotient_delta().unwrap().is_finite());
        for sample in audit.samples() {
            assert_eq!(sample.mode_norms().len(), 2);
            assert_eq!(sample.pairings().len(), 3);
            assert!(sample.max_l2_error().is_finite());
            assert!(sample.max_relative_l2_error().is_finite());
            assert!(sample.max_pairing_residual().is_finite());
            assert!(sample.max_observed_l2_continuity_quotient().is_finite());
            assert!(sample.max_pole_residual().is_finite());
            assert!(sample.max_archimedean_residual().is_finite());
            assert!(sample.max_prime_residual().is_finite());
            assert!(sample.max_pairing_asymmetry().is_finite());
            for norm in sample.mode_norms().iter().copied() {
                assert!(norm.direct_l2_norm().is_finite());
                assert!(norm.projected_l2_norm().is_finite());
                assert!(norm.l2_error().is_finite());
                assert!(norm.relative_l2_error().is_finite());
                assert!(norm.observed_max_absolute_error().is_finite());
            }
            for entry in sample.pairings().iter().copied() {
                assert!(entry.pairing_residual().is_finite());
                assert!(entry.l2_perturbation_scale().is_finite());
                assert!(entry.observed_l2_continuity_quotient().is_finite());
            }
        }
    }

    #[test]
    fn continuity_probe_rejects_invalid_parent_axis() {
        let modes = SineModeSet::new(vec![1]).unwrap();
        let config =
            SineContinuityProbeConfig::new(16, 16, CompactWeilPairingConfig::new(12, 12, 16));
        assert!(matches!(
            audit_finite_weil_sine_continuity(bump(), &modes, &[], config),
            Err(FiniteWeilSineContinuityError::EmptyParentDimensionSet)
        ));
        assert!(matches!(
            audit_finite_weil_sine_continuity(bump(), &modes, &[2, 2], config),
            Err(FiniteWeilSineContinuityError::ParentDimensionsNotStrictlyIncreasing { .. })
        ));
    }
}
