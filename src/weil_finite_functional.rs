//! Finite Riemann--Weil functional audit for one compact convolution square.
//!
//! Let `h1` be the compact test function validated by
//! `semilocal_compact_weil`, written in logarithmic coordinate as
//! `phi(t) = h1(exp(t))`. For real `h1` its multiplicative convolution square is
//!
//! `F(exp(t)) = theta(t) = integral phi(u+t) phi(u) du`.
//!
//! This module evaluates the source decomposition
//!
//! `psi(F) = Fhat(i/2) + Fhat(-i/2) - W_R(F) - sum_p W_p(F)`
//!
//! using the compact-support formulas of Connes--Consani. Compact support makes
//! the prime-power sum finite. The result is a numerical audit of one
//! manufactured convolution square; its sign is not promoted to a general Weil
//! positivity statement.

use std::f64::consts::PI;
use std::fmt;

use crate::quadrature::{GaussLegendreUnit, QuadratureError};
use crate::semilocal_compact_weil::{CompactWeilBoundaryAudit, CompactWeilTestFunction};
use crate::weil_boundary::WeilBoundaryError;

const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;

/// One non-archimedean prime-power term in the finite Weil functional.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrimePowerWeilTerm {
    integer: u64,
    prime: u64,
    exponent: u32,
    on_support_boundary: bool,
    theta_positive: f64,
    theta_negative: f64,
    contribution: f64,
}

impl PrimePowerWeilTerm {
    #[inline]
    pub const fn integer(self) -> u64 {
        self.integer
    }

    #[inline]
    pub const fn prime(self) -> u64 {
        self.prime
    }

    #[inline]
    pub const fn exponent(self) -> u32 {
        self.exponent
    }

    #[inline]
    pub const fn on_support_boundary(self) -> bool {
        self.on_support_boundary
    }

    #[inline]
    pub const fn theta_positive(self) -> f64 {
        self.theta_positive
    }

    #[inline]
    pub const fn theta_negative(self) -> f64 {
        self.theta_negative
    }

    #[inline]
    pub const fn contribution(self) -> f64 {
        self.contribution
    }

    #[inline]
    pub fn symmetry_residual(self) -> f64 {
        self.theta_positive - self.theta_negative
    }
}

/// Auditable source decomposition for one finite compact convolution square.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilFunctionalAudit {
    log_support_radius: f64,
    max_prime_power_argument: u64,
    autocorrelation_zero: f64,
    boundary: CompactWeilBoundaryAudit,
    pole_term: f64,
    archimedean_term: f64,
    prime_terms: Vec<PrimePowerWeilTerm>,
    prime_total: f64,
    functional_value: f64,
}

impl FiniteWeilFunctionalAudit {
    #[inline]
    pub fn log_support_radius(&self) -> f64 {
        self.log_support_radius
    }

    #[inline]
    pub fn max_prime_power_argument(&self) -> u64 {
        self.max_prime_power_argument
    }

    #[inline]
    pub fn autocorrelation_zero(&self) -> f64 {
        self.autocorrelation_zero
    }

    #[inline]
    pub fn boundary(&self) -> CompactWeilBoundaryAudit {
        self.boundary
    }

    #[inline]
    pub fn pole_term(&self) -> f64 {
        self.pole_term
    }

    /// Source archimedean distribution `W_R(F)`.
    #[inline]
    pub fn archimedean_term(&self) -> f64 {
        self.archimedean_term
    }

    #[inline]
    pub fn prime_terms(&self) -> &[PrimePowerWeilTerm] {
        &self.prime_terms
    }

    /// Sum of all finite non-archimedean terms represented by the support.
    #[inline]
    pub fn prime_total(&self) -> f64 {
        self.prime_total
    }

    /// `psi(F) = pole - W_R(F) - sum_p W_p(F)` in the source convention.
    #[inline]
    pub fn functional_value(&self) -> f64 {
        self.functional_value
    }

    pub fn max_prime_symmetry_residual(&self) -> f64 {
        self.prime_terms
            .iter()
            .map(|term| term.symmetry_residual().abs())
            .fold(0.0_f64, f64::max)
    }
}

#[derive(Debug)]
pub enum FiniteWeilFunctionalError {
    Quadrature(QuadratureError),
    Boundary(WeilBoundaryError),
    SupportRatioTooLarge { floor: u128 },
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilFunctionalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Quadrature(error) => write!(f, "Weil quadrature construction failed: {error:?}"),
            Self::Boundary(error) => {
                write!(f, "compact Weil test-function evaluation failed: {error}")
            }
            Self::SupportRatioTooLarge { floor } => write!(
                f,
                "compact support ratio requires a prime-power bound larger than u64: floor={floor}"
            ),
            Self::NonFiniteEvaluation { stage, value } => {
                write!(f, "non-finite finite-Weil value at {stage}: {value}")
            }
        }
    }
}

impl std::error::Error for FiniteWeilFunctionalError {}

impl From<QuadratureError> for FiniteWeilFunctionalError {
    fn from(value: QuadratureError) -> Self {
        Self::Quadrature(value)
    }
}

impl From<WeilBoundaryError> for FiniteWeilFunctionalError {
    fn from(value: WeilBoundaryError) -> Self {
        Self::Boundary(value)
    }
}

#[derive(Clone, Copy, Debug)]
struct ExactSupportRatio {
    numerator: u128,
    denominator: u128,
    floor: u64,
}

impl ExactSupportRatio {
    fn from_test_function(
        test_function: CompactWeilTestFunction,
    ) -> Result<Self, FiniteWeilFunctionalError> {
        let generator = test_function.generator();
        let lower = generator.lower();
        let upper = generator.upper();
        let numerator = u128::from(upper.numerator()) * u128::from(lower.denominator());
        let denominator = u128::from(upper.denominator()) * u128::from(lower.numerator());
        let floor = numerator / denominator;
        if floor > u128::from(u64::MAX) {
            return Err(FiniteWeilFunctionalError::SupportRatioTooLarge { floor });
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
}

struct CompactLogAutocorrelation {
    test_function: CompactWeilTestFunction,
    quadrature: GaussLegendreUnit,
    log_lower: f64,
    log_upper: f64,
    log_span: f64,
}

impl CompactLogAutocorrelation {
    fn new(
        test_function: CompactWeilTestFunction,
        quadrature_order: usize,
    ) -> Result<Self, FiniteWeilFunctionalError> {
        let support = test_function.generator().support();
        let log_lower = support.log_lower();
        let log_upper = support.log_upper();
        let log_span = log_upper - log_lower;
        let quadrature = GaussLegendreUnit::new(quadrature_order)?;
        Ok(Self {
            test_function,
            quadrature,
            log_lower,
            log_upper,
            log_span,
        })
    }

    fn value(&self, shift: f64) -> Result<f64, FiniteWeilFunctionalError> {
        if !shift.is_finite() {
            return Err(FiniteWeilFunctionalError::NonFiniteEvaluation {
                stage: "autocorrelation shift",
                value: shift,
            });
        }
        if shift.abs() >= self.log_span {
            return Ok(0.0);
        }

        // phi(u) and phi(u+shift) must both lie in [log_lower, log_upper].
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
            let left = self.test_function.q_value((u + shift).exp())?;
            let right = self.test_function.q_value(u.exp())?;
            total += weight * left * right;
        }
        let value = span * total;
        checked_finite("log-autocorrelation", value)?;
        Ok(value)
    }
}

/// Evaluate one compact convolution square with the finite Riemann--Weil source
/// decomposition.
///
/// `autocorrelation_order` controls the inner log-coordinate convolution,
/// `archimedean_order` controls the compact integral for `W_R`, and
/// `boundary_order` controls the two critical Mellin moments of `h1`.
pub fn audit_finite_weil_functional(
    test_function: CompactWeilTestFunction,
    autocorrelation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
) -> Result<FiniteWeilFunctionalAudit, FiniteWeilFunctionalError> {
    let ratio = ExactSupportRatio::from_test_function(test_function)?;
    let autocorrelation = CompactLogAutocorrelation::new(test_function, autocorrelation_order)?;
    let log_support_radius = autocorrelation.log_span;
    let theta_zero = autocorrelation.value(0.0)?;
    let boundary = test_function.audit_boundary_moments(boundary_order)?;
    let moments = boundary.moments();

    // For real h1 and F = h1^* * h1, multiplicative Fourier convolution gives
    // Fhat(+i/2) = Fhat(-i/2) = M_+(h1) M_-(h1).
    let pole_term = 2.0 * moments.plus_half * moments.minus_half;
    checked_finite("critical pole term", pole_term)?;

    let archimedean_term =
        archimedean_weil_term(&autocorrelation, theta_zero, archimedean_order, ratio)?;

    let mut prime_terms = Vec::new();
    let mut prime_total = 0.0_f64;
    for integer in 2..=ratio.floor {
        let Some((prime, exponent)) = prime_power_decomposition(integer) else {
            continue;
        };

        let on_support_boundary = ratio.integer_is_boundary(integer);
        let (theta_positive, theta_negative) = if on_support_boundary {
            // Exact compact support says the autocorrelation vanishes at the
            // endpoint even if binary64 logs differ by a last-bit rounding.
            (0.0, 0.0)
        } else {
            let shift = (integer as f64).ln();
            (
                autocorrelation.value(shift)?,
                autocorrelation.value(-shift)?,
            )
        };
        let contribution =
            (prime as f64).ln() / (integer as f64).sqrt() * (theta_positive + theta_negative);
        checked_finite("prime-power contribution", contribution)?;
        prime_total += contribution;
        prime_terms.push(PrimePowerWeilTerm {
            integer,
            prime,
            exponent,
            on_support_boundary,
            theta_positive,
            theta_negative,
            contribution,
        });
    }
    checked_finite("prime-power total", prime_total)?;

    let functional_value = pole_term - archimedean_term - prime_total;
    checked_finite("finite Weil functional", functional_value)?;

    Ok(FiniteWeilFunctionalAudit {
        log_support_radius,
        max_prime_power_argument: ratio.floor,
        autocorrelation_zero: theta_zero,
        boundary,
        pole_term,
        archimedean_term,
        prime_terms,
        prime_total,
        functional_value,
    })
}

/// Source formula (2.32) for the real-place distribution specialized to a
/// compact log-autocorrelation.
fn archimedean_weil_term(
    autocorrelation: &CompactLogAutocorrelation,
    theta_zero: f64,
    quadrature_order: usize,
    ratio: ExactSupportRatio,
) -> Result<f64, FiniteWeilFunctionalError> {
    let quadrature = GaussLegendreUnit::new(quadrature_order)?;
    let log_span = autocorrelation.log_span;
    let theta_sym_zero = 2.0 * theta_zero;

    let ratio_f64 = ratio.numerator as f64 / ratio.denominator as f64;
    let coefficient = EULER_MASCHERONI + (4.0 * PI * (ratio_f64 - 1.0) / (ratio_f64 + 1.0)).ln();
    checked_finite("archimedean compact coefficient", coefficient)?;

    let mut weighted_sum = 0.0_f64;
    for (&node, &weight) in quadrature.nodes().iter().zip(quadrature.weights().iter()) {
        let t = log_span * node;
        let theta_symmetric = autocorrelation.value(t)? + autocorrelation.value(-t)?;
        let exp_half = (0.5 * t).exp();

        // Algebraically this is exp(t/2)*theta_sym(t)-theta_sym(0).
        // Splitting off expm1 reduces cancellation near t=0.
        let numerator =
            (0.5 * t).exp_m1() * theta_sym_zero + exp_half * (theta_symmetric - theta_sym_zero);
        let denominator = 2.0 * t.sinh();
        let integrand = numerator / denominator;
        checked_finite("archimedean compact integrand", integrand)?;
        weighted_sum += weight * integrand;
    }

    let integral = log_span * weighted_sum;
    let value = 0.5 * theta_sym_zero * coefficient + integral;
    checked_finite("archimedean Weil term", value)?;
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

fn checked_finite(stage: &'static str, value: f64) -> Result<(), FiniteWeilFunctionalError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilFunctionalError::NonFiniteEvaluation { stage, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::{CompactArchimedeanBump, PositiveRational};

    fn fixture() -> CompactWeilTestFunction {
        let bump = CompactArchimedeanBump::new(
            PositiveRational::new(1, 2).unwrap(),
            PositiveRational::new(7, 2).unwrap(),
        )
        .unwrap();
        CompactWeilTestFunction::new(bump)
    }

    #[test]
    fn prime_power_classifier_rejects_mixed_composites() {
        assert_eq!(prime_power_decomposition(2), Some((2, 1)));
        assert_eq!(prime_power_decomposition(4), Some((2, 2)));
        assert_eq!(prime_power_decomposition(8), Some((2, 3)));
        assert_eq!(prime_power_decomposition(9), Some((3, 2)));
        assert_eq!(prime_power_decomposition(5), Some((5, 1)));
        assert_eq!(prime_power_decomposition(6), None);
        assert_eq!(prime_power_decomposition(12), None);
    }

    #[test]
    fn compact_ratio_drives_exact_prime_power_window() {
        let audit = audit_finite_weil_functional(fixture(), 96, 96, 128).unwrap();
        assert_eq!(audit.max_prime_power_argument(), 7);
        assert_eq!(
            audit
                .prime_terms()
                .iter()
                .map(|term| term.integer())
                .collect::<Vec<_>>(),
            vec![2, 3, 4, 5, 7]
        );
        let endpoint = audit.prime_terms().last().unwrap();
        assert!(endpoint.on_support_boundary());
        assert_eq!(endpoint.contribution(), 0.0);
    }

    #[test]
    fn real_autocorrelation_is_numerically_symmetric_at_prime_powers() {
        let audit = audit_finite_weil_functional(fixture(), 128, 96, 128).unwrap();
        assert!(audit.max_prime_symmetry_residual() <= 2.0e-14);
        assert!(audit.autocorrelation_zero() > 0.0);
    }

    #[test]
    fn manufactured_weil_decomposition_is_stable_and_finite() {
        let audit = audit_finite_weil_functional(fixture(), 128, 128, 128).unwrap();

        assert!(audit.boundary().satisfies(5.0e-12));
        assert!(audit.pole_term().abs() <= 1.0e-20);
        assert!((audit.archimedean_term() - (-0.014_576_55)).abs() <= 2.0e-7);
        assert!((audit.prime_total() - (-0.002_092_49)).abs() <= 2.0e-7);
        assert!((audit.functional_value() - 0.016_669_04).abs() <= 3.0e-7);
        assert_eq!(
            audit.functional_value(),
            audit.pole_term() - audit.archimedean_term() - audit.prime_total()
        );
    }

    #[test]
    fn functional_value_converges_under_quadrature_refinement() {
        let coarse = audit_finite_weil_functional(fixture(), 64, 64, 64).unwrap();
        let medium = audit_finite_weil_functional(fixture(), 96, 96, 96).unwrap();
        let fine = audit_finite_weil_functional(fixture(), 128, 128, 128).unwrap();

        let coarse_error = (coarse.functional_value() - fine.functional_value()).abs();
        let medium_error = (medium.functional_value() - fine.functional_value()).abs();
        assert!(medium_error < coarse_error);
        assert!(medium_error <= 5.0e-7);
    }
}
