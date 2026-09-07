//! Generic mixed Riemann--Weil pairing for compact admissible functions.
//!
//! This module evaluates the same finite compact-support source decomposition
//! used by the existing Legendre matrix, but the inputs are expressed through
//! [`CompactWeilEvaluand`] rather than a polynomial degree.  It is therefore
//! suitable for independent cross-checks with non-polynomial compact families.
//!
//! The result is still a finite numerical audit.  It does not establish
//! complete-space Weil positivity, density/completeness, Conjecture 4.1, or RH.

use std::f64::consts::PI;
use std::fmt;

use crate::quadrature::{GaussLegendreUnit, QuadratureError};
use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_boundary::{WeilBoundaryError, WeilBoundaryMoments};
use crate::weil_quadratic_matrix::CompactWeilBasisFunction;
use crate::weil_sine_direct_q::CompactSineWeilFunction;

const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;

/// Minimal contract needed by the finite compact mixed Weil decomposition.
pub trait CompactWeilEvaluand {
    fn bump(&self) -> CompactArchimedeanBump;

    fn value(&self, rho: f64) -> Result<f64, WeilBoundaryError>;

    fn boundary_moments(
        &self,
        quadrature_order: usize,
    ) -> Result<WeilBoundaryMoments, WeilBoundaryError>;
}

impl CompactWeilEvaluand for CompactWeilBasisFunction {
    #[inline]
    fn bump(&self) -> CompactArchimedeanBump {
        (*self).bump()
    }

    #[inline]
    fn value(&self, rho: f64) -> Result<f64, WeilBoundaryError> {
        (*self).value(rho)
    }

    #[inline]
    fn boundary_moments(
        &self,
        quadrature_order: usize,
    ) -> Result<WeilBoundaryMoments, WeilBoundaryError> {
        (*self).boundary_moments(quadrature_order)
    }
}

impl CompactWeilEvaluand for CompactSineWeilFunction {
    #[inline]
    fn bump(&self) -> CompactArchimedeanBump {
        (*self).bump()
    }

    #[inline]
    fn value(&self, rho: f64) -> Result<f64, WeilBoundaryError> {
        (*self).q_value(rho)
    }

    #[inline]
    fn boundary_moments(
        &self,
        quadrature_order: usize,
    ) -> Result<WeilBoundaryMoments, WeilBoundaryError> {
        (*self).boundary_moments(quadrature_order)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CompactWeilPairingConfig {
    correlation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
}

impl CompactWeilPairingConfig {
    #[inline]
    pub const fn new(
        correlation_order: usize,
        archimedean_order: usize,
        boundary_order: usize,
    ) -> Self {
        Self {
            correlation_order,
            archimedean_order,
            boundary_order,
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
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteCompactWeilPairingAudit {
    config: CompactWeilPairingConfig,
    max_prime_power_argument: u64,
    left_boundary_residual: f64,
    right_boundary_residual: f64,
    pole_term: f64,
    archimedean_term: f64,
    prime_total: f64,
    value: f64,
}

impl FiniteCompactWeilPairingAudit {
    #[inline]
    pub const fn config(self) -> CompactWeilPairingConfig {
        self.config
    }

    #[inline]
    pub const fn max_prime_power_argument(self) -> u64 {
        self.max_prime_power_argument
    }

    #[inline]
    pub const fn left_boundary_residual(self) -> f64 {
        self.left_boundary_residual
    }

    #[inline]
    pub const fn right_boundary_residual(self) -> f64 {
        self.right_boundary_residual
    }

    #[inline]
    pub const fn pole_term(self) -> f64 {
        self.pole_term
    }

    #[inline]
    pub const fn archimedean_term(self) -> f64 {
        self.archimedean_term
    }

    #[inline]
    pub const fn prime_total(self) -> f64 {
        self.prime_total
    }

    #[inline]
    pub const fn value(self) -> f64 {
        self.value
    }
}

#[derive(Debug)]
pub enum CompactWeilPairingError {
    Quadrature(QuadratureError),
    Boundary(WeilBoundaryError),
    SupportMismatch,
    SupportRatioTooLarge { floor: u128 },
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for CompactWeilPairingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Quadrature(error) => {
                write!(f, "compact Weil pairing quadrature failed: {error:?}")
            }
            Self::Boundary(error) => write!(f, "compact Weil evaluand failed: {error}"),
            Self::SupportMismatch => write!(
                f,
                "mixed compact Weil evaluands must share one exact support"
            ),
            Self::SupportRatioTooLarge { floor } => write!(
                f,
                "compact support ratio requires a prime-power bound larger than u64: floor={floor}"
            ),
            Self::NonFiniteEvaluation { stage, value } => {
                write!(
                    f,
                    "non-finite compact Weil pairing value at {stage}: {value}"
                )
            }
        }
    }
}

impl std::error::Error for CompactWeilPairingError {}

impl From<QuadratureError> for CompactWeilPairingError {
    fn from(value: QuadratureError) -> Self {
        Self::Quadrature(value)
    }
}

impl From<WeilBoundaryError> for CompactWeilPairingError {
    fn from(value: WeilBoundaryError) -> Self {
        Self::Boundary(value)
    }
}

#[derive(Clone, Copy)]
struct ExactSupportRatio {
    numerator: u128,
    denominator: u128,
    floor: u64,
}

impl ExactSupportRatio {
    fn from_bump(bump: CompactArchimedeanBump) -> Result<Self, CompactWeilPairingError> {
        let lower = bump.lower();
        let upper = bump.upper();
        let numerator = u128::from(upper.numerator()) * u128::from(lower.denominator());
        let denominator = u128::from(upper.denominator()) * u128::from(lower.numerator());
        let floor = numerator / denominator;
        if floor > u128::from(u64::MAX) {
            return Err(CompactWeilPairingError::SupportRatioTooLarge { floor });
        }
        Ok(Self {
            numerator,
            denominator,
            floor: floor as u64,
        })
    }

    #[inline]
    fn integer_is_boundary(self, integer: u64) -> bool {
        u128::from(integer) * self.denominator == self.numerator
    }

    #[inline]
    fn as_f64(self) -> f64 {
        self.numerator as f64 / self.denominator as f64
    }
}

struct MixedLogCorrelation<'a, L: ?Sized, R: ?Sized> {
    left: &'a L,
    right: &'a R,
    quadrature: GaussLegendreUnit,
    log_lower: f64,
    log_upper: f64,
    log_span: f64,
}

impl<'a, L, R> MixedLogCorrelation<'a, L, R>
where
    L: CompactWeilEvaluand + ?Sized,
    R: CompactWeilEvaluand + ?Sized,
{
    fn new(
        left: &'a L,
        right: &'a R,
        quadrature_order: usize,
    ) -> Result<Self, CompactWeilPairingError> {
        if left.bump() != right.bump() {
            return Err(CompactWeilPairingError::SupportMismatch);
        }
        let support = left.bump().support();
        let log_lower = support.log_lower();
        let log_upper = support.log_upper();
        Ok(Self {
            left,
            right,
            quadrature: GaussLegendreUnit::new(quadrature_order)?,
            log_lower,
            log_upper,
            log_span: log_upper - log_lower,
        })
    }

    fn value(&self, shift: f64) -> Result<f64, CompactWeilPairingError> {
        if !shift.is_finite() {
            return Err(CompactWeilPairingError::NonFiniteEvaluation {
                stage: "mixed correlation shift",
                value: shift,
            });
        }
        if shift.abs() >= self.log_span {
            return Ok(0.0);
        }

        let lower = self.log_lower.max(self.log_lower - shift);
        let upper = self.log_upper.min(self.log_upper - shift);
        if upper <= lower {
            return Ok(0.0);
        }

        let span = upper - lower;
        let mut total = 0.0_f64;
        for (&node, &weight) in self
            .quadrature
            .nodes()
            .iter()
            .zip(self.quadrature.weights().iter())
        {
            let u = lower + span * node;
            let left = self.left.value(u.exp())?;
            let right = self.right.value((u + shift).exp())?;
            total += weight * left * right;
        }
        let value = span * total;
        checked_finite("mixed log correlation", value)?;
        Ok(value)
    }
}

pub fn audit_compact_weil_pairing<L, R>(
    left: &L,
    right: &R,
    config: CompactWeilPairingConfig,
) -> Result<FiniteCompactWeilPairingAudit, CompactWeilPairingError>
where
    L: CompactWeilEvaluand + ?Sized,
    R: CompactWeilEvaluand + ?Sized,
{
    if left.bump() != right.bump() {
        return Err(CompactWeilPairingError::SupportMismatch);
    }
    let ratio = ExactSupportRatio::from_bump(left.bump())?;
    let correlation = MixedLogCorrelation::new(left, right, config.correlation_order())?;
    let left_moments = left.boundary_moments(config.boundary_order())?;
    let right_moments = right.boundary_moments(config.boundary_order())?;
    let left_boundary_residual = left_moments
        .plus_half
        .abs()
        .max(left_moments.minus_half.abs());
    let right_boundary_residual = right_moments
        .plus_half
        .abs()
        .max(right_moments.minus_half.abs());

    let pole_term = left_moments.minus_half * right_moments.plus_half
        + left_moments.plus_half * right_moments.minus_half;
    checked_finite("mixed critical pole term", pole_term)?;

    let theta_zero = correlation.value(0.0)?;
    let archimedean_term =
        mixed_archimedean_term(&correlation, theta_zero, ratio, config.archimedean_order())?;

    let mut prime_total = 0.0_f64;
    for integer in 2..=ratio.floor {
        let Some((prime, _exponent)) = prime_power_decomposition(integer) else {
            continue;
        };
        let contribution = if ratio.integer_is_boundary(integer) {
            0.0
        } else {
            let shift = (integer as f64).ln();
            let theta_positive = correlation.value(shift)?;
            let theta_negative = correlation.value(-shift)?;
            (prime as f64).ln() / (integer as f64).sqrt() * (theta_positive + theta_negative)
        };
        checked_finite("mixed prime-power contribution", contribution)?;
        prime_total += contribution;
    }
    checked_finite("mixed prime-power total", prime_total)?;

    let value = pole_term - archimedean_term - prime_total;
    checked_finite("mixed Weil pairing", value)?;
    Ok(FiniteCompactWeilPairingAudit {
        config,
        max_prime_power_argument: ratio.floor,
        left_boundary_residual,
        right_boundary_residual,
        pole_term,
        archimedean_term,
        prime_total,
        value,
    })
}

fn mixed_archimedean_term<L, R>(
    correlation: &MixedLogCorrelation<'_, L, R>,
    theta_zero: f64,
    ratio: ExactSupportRatio,
    quadrature_order: usize,
) -> Result<f64, CompactWeilPairingError>
where
    L: CompactWeilEvaluand + ?Sized,
    R: CompactWeilEvaluand + ?Sized,
{
    let quadrature = GaussLegendreUnit::new(quadrature_order)?;
    let theta_sym_zero = 2.0 * theta_zero;
    let ratio_f64 = ratio.as_f64();
    let coefficient = EULER_MASCHERONI + (4.0 * PI * (ratio_f64 - 1.0) / (ratio_f64 + 1.0)).ln();
    checked_finite("mixed archimedean coefficient", coefficient)?;

    let mut weighted_sum = 0.0_f64;
    for (&node, &weight) in quadrature.nodes().iter().zip(quadrature.weights().iter()) {
        let t = correlation.log_span * node;
        let theta_symmetric = correlation.value(t)? + correlation.value(-t)?;
        let exp_half = (0.5 * t).exp();
        let numerator =
            (0.5 * t).exp_m1() * theta_sym_zero + exp_half * (theta_symmetric - theta_sym_zero);
        let denominator = 2.0 * t.sinh();
        let integrand = numerator / denominator;
        checked_finite("mixed archimedean integrand", integrand)?;
        weighted_sum += weight * integrand;
    }

    let value = 0.5 * theta_sym_zero * coefficient + correlation.log_span * weighted_sum;
    checked_finite("mixed archimedean term", value)?;
    Ok(value)
}

fn prime_power_decomposition(mut value: u64) -> Option<(u64, u32)> {
    if value < 2 {
        return None;
    }
    let original = value;
    let mut divisor = 2_u64;
    while divisor <= value / divisor {
        if value.is_multiple_of(divisor) {
            let prime = divisor;
            let mut exponent = 0_u32;
            while value.is_multiple_of(prime) {
                value /= prime;
                exponent += 1;
            }
            return (value == 1).then_some((prime, exponent));
        }
        divisor = if divisor == 2 { 3 } else { divisor + 2 };
    }
    Some((original, 1))
}

fn checked_finite(stage: &'static str, value: f64) -> Result<(), CompactWeilPairingError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(CompactWeilPairingError::NonFiniteEvaluation { stage, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::PositiveRational;
    use crate::weil_quadratic_matrix::audit_finite_weil_quadratic_matrix;

    fn bump() -> CompactArchimedeanBump {
        CompactArchimedeanBump::new(
            PositiveRational::new(1, 2).unwrap(),
            PositiveRational::new(7, 2).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn generic_pairing_regresses_existing_legendre_matrix() {
        let bump = bump();
        let config = CompactWeilPairingConfig::new(20, 20, 28);
        let matrix = audit_finite_weil_quadratic_matrix(bump, 2, 20, 20, 28).unwrap();
        let degree_zero = CompactWeilBasisFunction::new(bump, 0);
        let degree_one = CompactWeilBasisFunction::new(bump, 1);

        let diagonal = audit_compact_weil_pairing(&degree_zero, &degree_zero, config).unwrap();
        let expected_diagonal = matrix.entry(0, 0).unwrap();
        assert!(
            (diagonal.value() - expected_diagonal).abs()
                <= 5.0e-13 * expected_diagonal.abs().max(1.0)
        );

        let forward = audit_compact_weil_pairing(&degree_zero, &degree_one, config).unwrap();
        let reverse = audit_compact_weil_pairing(&degree_one, &degree_zero, config).unwrap();
        let symmetric = 0.5 * (forward.value() + reverse.value());
        let expected_mixed = matrix.entry(0, 1).unwrap();
        assert!((symmetric - expected_mixed).abs() <= 5.0e-13 * expected_mixed.abs().max(1.0));
    }
}
