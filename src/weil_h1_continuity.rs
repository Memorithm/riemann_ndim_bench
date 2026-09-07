//! Deterministic finite-support H1 continuity envelope for the mixed Weil pairing.
//!
//! This module turns upper bounds on `L2(d*rho)` norms and logarithmic
//! derivatives into an explicit upper bound for the change of the finite
//! compact Riemann--Weil pairing used by `weil_compact_pairing`.
//!
//! In logarithmic coordinate `u = log rho`, write
//!
//! `C(f,g)(t) = integral f(u) g(u+t) du`.
//!
//! For target/approximant pairs `(f, f_N)` and `(g, g_N)`, Cauchy--Schwarz gives
//!
//! `|Delta C(t)| <= A`,
//!
//! with `A = e_f N_g + N_f e_g`, where `e_f=||f-f_N||_2` and
//! `N_f=max(||f||_2,||f_N||_2)`, and similarly on the right.  If logarithmic
//! derivatives are controlled, the derivative of the correlation difference
//! is bounded by `D`, the smaller of the valid left- and right-derivative
//! Cauchy--Schwarz estimates.
//!
//! These two quantities yield explicit bounds for the pole, real-place, and
//! prime-power pieces of the source-locked finite Weil decomposition.  The
//! formulas below are mathematical upper envelopes provided that the caller's
//! norm inputs really are upper bounds.  Their binary64 evaluation is not a
//! formal proof assistant certificate, and this module does not prove density,
//! complete-space Weil positivity, Conjecture 4.1, or RH.

use std::f64::consts::PI;
use std::fmt;

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;

const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;

/// Declared upper bounds for one target function and one approximant.
///
/// Derivatives are with respect to `u = log rho`, i.e. `d/du = rho d/drho`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct H1ApproximationUpperBounds {
    target_l2: f64,
    approximant_l2: f64,
    error_l2: f64,
    target_log_derivative_l2: f64,
    approximant_log_derivative_l2: f64,
    derivative_error_l2: f64,
}

impl H1ApproximationUpperBounds {
    pub fn new(
        target_l2: f64,
        approximant_l2: f64,
        error_l2: f64,
        target_log_derivative_l2: f64,
        approximant_log_derivative_l2: f64,
        derivative_error_l2: f64,
    ) -> Result<Self, FiniteWeilH1ContinuityError> {
        for (field, value) in [
            ("target_l2", target_l2),
            ("approximant_l2", approximant_l2),
            ("error_l2", error_l2),
            ("target_log_derivative_l2", target_log_derivative_l2),
            (
                "approximant_log_derivative_l2",
                approximant_log_derivative_l2,
            ),
            ("derivative_error_l2", derivative_error_l2),
        ] {
            checked_nonnegative(field, value)?;
        }
        Ok(Self {
            target_l2,
            approximant_l2,
            error_l2,
            target_log_derivative_l2,
            approximant_log_derivative_l2,
            derivative_error_l2,
        })
    }

    #[inline]
    pub const fn target_l2(self) -> f64 {
        self.target_l2
    }

    #[inline]
    pub const fn approximant_l2(self) -> f64 {
        self.approximant_l2
    }

    #[inline]
    pub const fn error_l2(self) -> f64 {
        self.error_l2
    }

    #[inline]
    pub const fn target_log_derivative_l2(self) -> f64 {
        self.target_log_derivative_l2
    }

    #[inline]
    pub const fn approximant_log_derivative_l2(self) -> f64 {
        self.approximant_log_derivative_l2
    }

    #[inline]
    pub const fn derivative_error_l2(self) -> f64 {
        self.derivative_error_l2
    }

    #[inline]
    pub fn max_l2(self) -> f64 {
        self.target_l2.max(self.approximant_l2)
    }

    #[inline]
    pub fn max_log_derivative_l2(self) -> f64 {
        self.target_log_derivative_l2
            .max(self.approximant_log_derivative_l2)
    }
}

/// Support-dependent coefficients in the H1 continuity envelope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WeilH1ContinuityConstants {
    support_ratio: f64,
    log_span: f64,
    pole_value_coefficient: f64,
    prime_value_coefficient: f64,
    archimedean_value_coefficient: f64,
    archimedean_derivative_coefficient: f64,
    total_value_coefficient: f64,
}

impl WeilH1ContinuityConstants {
    #[inline]
    pub const fn support_ratio(self) -> f64 {
        self.support_ratio
    }

    #[inline]
    pub const fn log_span(self) -> f64 {
        self.log_span
    }

    #[inline]
    pub const fn pole_value_coefficient(self) -> f64 {
        self.pole_value_coefficient
    }

    #[inline]
    pub const fn prime_value_coefficient(self) -> f64 {
        self.prime_value_coefficient
    }

    #[inline]
    pub const fn archimedean_value_coefficient(self) -> f64 {
        self.archimedean_value_coefficient
    }

    #[inline]
    pub const fn archimedean_derivative_coefficient(self) -> f64 {
        self.archimedean_derivative_coefficient
    }

    #[inline]
    pub const fn total_value_coefficient(self) -> f64 {
        self.total_value_coefficient
    }
}

/// Auditable component-wise upper envelope for one mixed-pairing perturbation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteWeilH1ContinuityBound {
    left: H1ApproximationUpperBounds,
    right: H1ApproximationUpperBounds,
    constants: WeilH1ContinuityConstants,
    correlation_value_error_bound: f64,
    correlation_derivative_error_bound: f64,
    pole_error_bound: f64,
    archimedean_error_bound: f64,
    prime_error_bound: f64,
    total_pairing_error_bound: f64,
}

impl FiniteWeilH1ContinuityBound {
    #[inline]
    pub const fn left(self) -> H1ApproximationUpperBounds {
        self.left
    }

    #[inline]
    pub const fn right(self) -> H1ApproximationUpperBounds {
        self.right
    }

    #[inline]
    pub const fn constants(self) -> WeilH1ContinuityConstants {
        self.constants
    }

    #[inline]
    pub const fn correlation_value_error_bound(self) -> f64 {
        self.correlation_value_error_bound
    }

    #[inline]
    pub const fn correlation_derivative_error_bound(self) -> f64 {
        self.correlation_derivative_error_bound
    }

    #[inline]
    pub const fn pole_error_bound(self) -> f64 {
        self.pole_error_bound
    }

    #[inline]
    pub const fn archimedean_error_bound(self) -> f64 {
        self.archimedean_error_bound
    }

    #[inline]
    pub const fn prime_error_bound(self) -> f64 {
        self.prime_error_bound
    }

    #[inline]
    pub const fn total_pairing_error_bound(self) -> f64 {
        self.total_pairing_error_bound
    }
}

#[derive(Debug)]
pub enum FiniteWeilH1ContinuityError {
    InvalidUpperBound { field: &'static str, value: f64 },
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilH1ContinuityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidUpperBound { field, value } => write!(
                f,
                "finite Weil H1 continuity input must be finite and non-negative: {field}={value}"
            ),
            Self::NonFiniteEvaluation { stage, value } => write!(
                f,
                "non-finite finite-Weil H1 continuity value at {stage}: {value}"
            ),
        }
    }
}

impl std::error::Error for FiniteWeilH1ContinuityError {}

/// Compute the deterministic finite-support H1 continuity envelope.
///
/// The caller is responsible for supplying genuine upper bounds in `left` and
/// `right`.  Under that contract the returned component bounds follow from
/// Cauchy--Schwarz, translation invariance of `L2(R)`, the one-dimensional weak
/// derivative identity for correlations, and `sinh(t) >= t` for `t >= 0`.
pub fn bound_finite_weil_h1_pairing_perturbation(
    bump: CompactArchimedeanBump,
    left: H1ApproximationUpperBounds,
    right: H1ApproximationUpperBounds,
) -> Result<FiniteWeilH1ContinuityBound, FiniteWeilH1ContinuityError> {
    let constants = continuity_constants(bump)?;

    let correlation_value_error_bound =
        left.error_l2 * right.max_l2() + left.max_l2() * right.error_l2;

    // Differentiate the translated right factor, or integrate by parts and
    // differentiate the left factor. Both are valid H1 bounds, so their minimum
    // is still a valid upper bound.
    let derivative_on_right =
        left.error_l2 * right.max_log_derivative_l2() + left.max_l2() * right.derivative_error_l2;
    let derivative_on_left =
        left.derivative_error_l2 * right.max_l2() + left.max_log_derivative_l2() * right.error_l2;
    let correlation_derivative_error_bound = derivative_on_right.min(derivative_on_left);

    let pole_error_bound = constants.pole_value_coefficient * correlation_value_error_bound;
    let prime_error_bound = constants.prime_value_coefficient * correlation_value_error_bound;
    let archimedean_error_bound = constants.archimedean_value_coefficient
        * correlation_value_error_bound
        + constants.archimedean_derivative_coefficient * correlation_derivative_error_bound;
    let total_pairing_error_bound = pole_error_bound + archimedean_error_bound + prime_error_bound;

    for (stage, value) in [
        (
            "correlation value error bound",
            correlation_value_error_bound,
        ),
        (
            "correlation derivative error bound",
            correlation_derivative_error_bound,
        ),
        ("pole error bound", pole_error_bound),
        ("archimedean error bound", archimedean_error_bound),
        ("prime error bound", prime_error_bound),
        ("total pairing error bound", total_pairing_error_bound),
    ] {
        checked_finite(stage, value)?;
    }

    Ok(FiniteWeilH1ContinuityBound {
        left,
        right,
        constants,
        correlation_value_error_bound,
        correlation_derivative_error_bound,
        pole_error_bound,
        archimedean_error_bound,
        prime_error_bound,
        total_pairing_error_bound,
    })
}

fn continuity_constants(
    bump: CompactArchimedeanBump,
) -> Result<WeilH1ContinuityConstants, FiniteWeilH1ContinuityError> {
    let lower = bump.lower().as_f64();
    let upper = bump.upper().as_f64();
    let support_ratio = upper / lower;
    let log_span = support_ratio.ln();
    let sqrt_ratio = support_ratio.sqrt();

    // Critical Mellin moments have L2 weights exp(+-u/2). Their two weight
    // norms multiply to sqrt(R) - 1/sqrt(R), and the pole term has two products.
    let pole_value_coefficient = 2.0 * (sqrt_ratio - sqrt_ratio.recip());

    // The source sum contains only prime powers n < R. Bounding ln(p) <= ln(R)
    // and sum_{n=2}^m n^(-1/2) <= 2 sqrt(m) <= 2 sqrt(R), with the two
    // correlation shifts, gives the safe all-integer envelope below.
    let prime_value_coefficient = 4.0 * sqrt_ratio * log_span;

    // For the real-place term, Delta theta_sym(0) <= 2A and
    // |Delta theta_sym(t)-Delta theta_sym(0)| <= 2Dt. Using sinh(t) >= t,
    // exp(t/2)-1 <= (t/2) exp(L/2), and exp(t/2) <= exp(L/2) on [0,L]
    // yields these explicit coefficients.
    let source_coefficient =
        EULER_MASCHERONI + (4.0 * PI * (support_ratio - 1.0) / (support_ratio + 1.0)).ln();
    let archimedean_value_coefficient = source_coefficient.abs() + 0.5 * log_span * sqrt_ratio;
    let archimedean_derivative_coefficient = log_span * sqrt_ratio;
    let total_value_coefficient =
        pole_value_coefficient + prime_value_coefficient + archimedean_value_coefficient;

    for (stage, value) in [
        ("support ratio", support_ratio),
        ("log support span", log_span),
        ("pole value coefficient", pole_value_coefficient),
        ("prime value coefficient", prime_value_coefficient),
        (
            "archimedean value coefficient",
            archimedean_value_coefficient,
        ),
        (
            "archimedean derivative coefficient",
            archimedean_derivative_coefficient,
        ),
        ("total value coefficient", total_value_coefficient),
    ] {
        checked_finite(stage, value)?;
    }

    Ok(WeilH1ContinuityConstants {
        support_ratio,
        log_span,
        pole_value_coefficient,
        prime_value_coefficient,
        archimedean_value_coefficient,
        archimedean_derivative_coefficient,
        total_value_coefficient,
    })
}

fn checked_nonnegative(field: &'static str, value: f64) -> Result<(), FiniteWeilH1ContinuityError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(FiniteWeilH1ContinuityError::InvalidUpperBound { field, value })
    }
}

fn checked_finite(stage: &'static str, value: f64) -> Result<(), FiniteWeilH1ContinuityError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(FiniteWeilH1ContinuityError::NonFiniteEvaluation { stage, value })
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

    fn bounds(
        target_l2: f64,
        approximant_l2: f64,
        error_l2: f64,
        target_d1: f64,
        approximant_d1: f64,
        error_d1: f64,
    ) -> H1ApproximationUpperBounds {
        H1ApproximationUpperBounds::new(
            target_l2,
            approximant_l2,
            error_l2,
            target_d1,
            approximant_d1,
            error_d1,
        )
        .unwrap()
    }

    #[test]
    fn ratio_seven_constants_match_the_closed_envelope() {
        let left = bounds(1.0, 1.0, 0.0, 1.0, 1.0, 0.0);
        let audit = bound_finite_weil_h1_pairing_perturbation(bump(), left, left).unwrap();
        let constants = audit.constants();
        assert!((constants.support_ratio() - 7.0).abs() <= 1.0e-14);
        assert!((constants.log_span() - 7.0_f64.ln()).abs() <= 1.0e-14);
        assert!((constants.pole_value_coefficient() - 4.535_573_676_110_727).abs() <= 1.0e-12);
        assert!((constants.prime_value_coefficient() - 20.593_577_312_307_954).abs() <= 1.0e-12);
        assert!(
            (constants.archimedean_value_coefficient() - 5.394_755_003_457_536).abs() <= 1.0e-12
        );
        assert!(
            (constants.archimedean_derivative_coefficient() - 5.148_394_328_076_988).abs()
                <= 1.0e-12
        );
    }

    #[test]
    fn exact_zero_input_errors_force_zero_pairing_envelope() {
        let left = bounds(2.0, 2.0, 0.0, 3.0, 3.0, 0.0);
        let right = bounds(4.0, 4.0, 0.0, 5.0, 5.0, 0.0);
        let audit = bound_finite_weil_h1_pairing_perturbation(bump(), left, right).unwrap();
        assert_eq!(audit.correlation_value_error_bound(), 0.0);
        assert_eq!(audit.correlation_derivative_error_bound(), 0.0);
        assert_eq!(audit.pole_error_bound(), 0.0);
        assert_eq!(audit.archimedean_error_bound(), 0.0);
        assert_eq!(audit.prime_error_bound(), 0.0);
        assert_eq!(audit.total_pairing_error_bound(), 0.0);
    }

    #[test]
    fn mixed_envelope_is_symmetric_under_left_right_exchange() {
        let left = bounds(2.0, 2.1, 0.1, 3.0, 3.2, 0.25);
        let right = bounds(1.3, 1.4, 0.08, 2.2, 2.3, 0.12);
        let forward = bound_finite_weil_h1_pairing_perturbation(bump(), left, right).unwrap();
        let reverse = bound_finite_weil_h1_pairing_perturbation(bump(), right, left).unwrap();
        assert!(
            (forward.correlation_value_error_bound() - reverse.correlation_value_error_bound())
                .abs()
                <= 1.0e-14
        );
        assert!(
            (forward.correlation_derivative_error_bound()
                - reverse.correlation_derivative_error_bound())
            .abs()
                <= 1.0e-14
        );
        assert!(
            (forward.total_pairing_error_bound() - reverse.total_pairing_error_bound()).abs()
                <= 1.0e-12
        );
    }

    #[test]
    fn positive_declared_errors_produce_finite_component_bounds() {
        let left = bounds(2.0, 1.9, 0.2, 4.0, 3.8, 0.3);
        let right = bounds(1.5, 1.6, 0.1, 2.5, 2.7, 0.2);
        let audit = bound_finite_weil_h1_pairing_perturbation(bump(), left, right).unwrap();
        assert!(audit.correlation_value_error_bound().is_finite());
        assert!(audit.correlation_derivative_error_bound().is_finite());
        assert!(audit.pole_error_bound().is_finite());
        assert!(audit.archimedean_error_bound().is_finite());
        assert!(audit.prime_error_bound().is_finite());
        assert!(audit.total_pairing_error_bound().is_finite());
        assert!(audit.total_pairing_error_bound() >= audit.pole_error_bound());
        assert!(audit.total_pairing_error_bound() >= audit.archimedean_error_bound());
        assert!(audit.total_pairing_error_bound() >= audit.prime_error_bound());
    }

    #[test]
    fn non_finite_or_negative_declared_bounds_fail_closed() {
        assert!(matches!(
            H1ApproximationUpperBounds::new(1.0, 1.0, -1.0, 1.0, 1.0, 1.0),
            Err(FiniteWeilH1ContinuityError::InvalidUpperBound {
                field: "error_l2",
                ..
            })
        ));
        assert!(matches!(
            H1ApproximationUpperBounds::new(1.0, f64::NAN, 1.0, 1.0, 1.0, 1.0),
            Err(FiniteWeilH1ContinuityError::InvalidUpperBound {
                field: "approximant_l2",
                ..
            })
        ));
    }
}
