//! Closed-form H1 envelopes for compact sine Taylor approximants.
//!
//! This module supplies concrete upper-bound inputs to `weil_h1_continuity`
//! without estimating the needed norms by quadrature.  In the affine support
//! coordinate `t=(rho-a)/(b-a)`, it compares
//!
//! `p(t) = sin(mode*pi*t)`
//!
//! with the degree-`n` Taylor polynomial of `p` about `t=1/2`.  Taylor's
//! theorem gives explicit uniform bounds for the remainder and its first three
//! `t` derivatives.  Combining those with closed majorants for the first three
//! derivatives of the standard compact bump yields upper envelopes for
//!
//! `h = Q(bump*p)`, `D h`, `h-h_n`, and `D(h-h_n)`,
//!
//! where `D=rho d/drho` and `Q=-D^2+1/4`.
//!
//! The mathematical inequalities are analytic.  Their constants are evaluated
//! here in binary64 and are not outward-rounded interval certificates or proof
//! assistant objects.  The module therefore does not claim machine-certified
//! real arithmetic, density/completeness, full Weil positivity, Conjecture 4.1,
//! or RH.

use std::f64::consts::PI;
use std::fmt;

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_h1_continuity::{
    FiniteWeilH1ContinuityBound, FiniteWeilH1ContinuityError, H1ApproximationUpperBounds,
    bound_finite_weil_h1_pairing_perturbation,
};

/// Uniform derivative majorants for `exp(-1/(t(1-t)))` on `0<t<1`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StandardBumpDerivativeEnvelope {
    value: f64,
    first: f64,
    second: f64,
    third: f64,
}

impl StandardBumpDerivativeEnvelope {
    #[inline]
    pub const fn value(self) -> f64 {
        self.value
    }

    #[inline]
    pub const fn first(self) -> f64 {
        self.first
    }

    #[inline]
    pub const fn second(self) -> f64 {
        self.second
    }

    #[inline]
    pub const fn third(self) -> f64 {
        self.third
    }

    #[inline]
    fn as_array(self) -> [f64; 4] {
        [self.value, self.first, self.second, self.third]
    }
}

/// Uniform Taylor-remainder bounds for derivatives of orders 0 through 3.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SineTaylorRemainderEnvelope {
    derivative_bounds: [f64; 4],
}

impl SineTaylorRemainderEnvelope {
    #[inline]
    pub const fn value(self) -> f64 {
        self.derivative_bounds[0]
    }

    #[inline]
    pub const fn first(self) -> f64 {
        self.derivative_bounds[1]
    }

    #[inline]
    pub const fn second(self) -> f64 {
        self.derivative_bounds[2]
    }

    #[inline]
    pub const fn third(self) -> f64 {
        self.derivative_bounds[3]
    }

    #[inline]
    pub const fn derivative_bounds(self) -> [f64; 4] {
        self.derivative_bounds
    }
}

/// Closed upper envelopes for one sine target and one centered Taylor approximant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSineTaylorH1Envelope {
    mode: usize,
    degree: usize,
    omega: f64,
    log_support_span: f64,
    max_log_scale: f64,
    bump_derivatives: StandardBumpDerivativeEnvelope,
    remainder: SineTaylorRemainderEnvelope,
    target_q_sup: f64,
    target_log_derivative_q_sup: f64,
    error_q_sup: f64,
    error_log_derivative_q_sup: f64,
    h1_bounds: H1ApproximationUpperBounds,
}

impl FiniteSineTaylorH1Envelope {
    #[inline]
    pub const fn mode(self) -> usize {
        self.mode
    }

    #[inline]
    pub const fn degree(self) -> usize {
        self.degree
    }

    #[inline]
    pub const fn omega(self) -> f64 {
        self.omega
    }

    #[inline]
    pub const fn log_support_span(self) -> f64 {
        self.log_support_span
    }

    #[inline]
    pub const fn max_log_scale(self) -> f64 {
        self.max_log_scale
    }

    #[inline]
    pub const fn bump_derivatives(self) -> StandardBumpDerivativeEnvelope {
        self.bump_derivatives
    }

    #[inline]
    pub const fn remainder(self) -> SineTaylorRemainderEnvelope {
        self.remainder
    }

    #[inline]
    pub const fn target_q_sup(self) -> f64 {
        self.target_q_sup
    }

    #[inline]
    pub const fn target_log_derivative_q_sup(self) -> f64 {
        self.target_log_derivative_q_sup
    }

    #[inline]
    pub const fn error_q_sup(self) -> f64 {
        self.error_q_sup
    }

    #[inline]
    pub const fn error_log_derivative_q_sup(self) -> f64 {
        self.error_log_derivative_q_sup
    }

    #[inline]
    pub const fn h1_bounds(self) -> H1ApproximationUpperBounds {
        self.h1_bounds
    }
}

/// One direct bridge from two concrete Taylor envelopes to the H1 Weil bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSineTaylorPairingEnvelope {
    left: FiniteSineTaylorH1Envelope,
    right: FiniteSineTaylorH1Envelope,
    pairing: FiniteWeilH1ContinuityBound,
}

impl FiniteSineTaylorPairingEnvelope {
    #[inline]
    pub const fn left(self) -> FiniteSineTaylorH1Envelope {
        self.left
    }

    #[inline]
    pub const fn right(self) -> FiniteSineTaylorH1Envelope {
        self.right
    }

    #[inline]
    pub const fn pairing(self) -> FiniteWeilH1ContinuityBound {
        self.pairing
    }
}

#[derive(Debug)]
pub enum FiniteSineTaylorH1Error {
    ZeroMode,
    DegreeTooSmall { degree: usize },
    NonFiniteEvaluation { stage: &'static str, value: f64 },
    H1(FiniteWeilH1ContinuityError),
}

impl fmt::Display for FiniteSineTaylorH1Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroMode => write!(f, "sine Taylor H1 mode must be positive"),
            Self::DegreeTooSmall { degree } => write!(
                f,
                "sine Taylor H1 degree must be at least two to control three derivatives: {degree}"
            ),
            Self::NonFiniteEvaluation { stage, value } => {
                write!(f, "non-finite sine Taylor H1 envelope at {stage}: {value}")
            }
            Self::H1(error) => write!(f, "finite Weil H1 envelope failed: {error}"),
        }
    }
}

impl std::error::Error for FiniteSineTaylorH1Error {}

impl From<FiniteWeilH1ContinuityError> for FiniteSineTaylorH1Error {
    fn from(value: FiniteWeilH1ContinuityError) -> Self {
        Self::H1(value)
    }
}

/// Derive closed H1 upper envelopes for one compact sine and its Taylor polynomial.
pub fn derive_sine_taylor_h1_envelope(
    bump: CompactArchimedeanBump,
    mode: usize,
    degree: usize,
) -> Result<FiniteSineTaylorH1Envelope, FiniteSineTaylorH1Error> {
    if mode == 0 {
        return Err(FiniteSineTaylorH1Error::ZeroMode);
    }
    if degree < 2 {
        return Err(FiniteSineTaylorH1Error::DegreeTooSmall { degree });
    }

    let omega = mode as f64 * PI;
    checked_nonnegative("angular frequency", omega)?;

    let support = bump.support();
    let width = support.upper() - support.lower();
    let log_support_span = support.log_upper() - support.log_lower();
    let max_log_scale = support.upper() / width;
    let sqrt_log_span = log_support_span.sqrt();
    for (stage, value) in [
        ("support width", width),
        ("log support span", log_support_span),
        ("maximum logarithmic scale", max_log_scale),
        ("square root log support span", sqrt_log_span),
    ] {
        checked_positive(stage, value)?;
    }

    let bump_derivatives = standard_bump_derivative_envelope()?;
    let target_profile = [1.0, omega, omega * omega, omega * omega * omega];
    for (stage, value) in [
        ("target first derivative bound", target_profile[1]),
        ("target second derivative bound", target_profile[2]),
        ("target third derivative bound", target_profile[3]),
    ] {
        checked_nonnegative(stage, value)?;
    }
    let target_product = product_derivative_envelope(bump_derivatives.as_array(), target_profile)?;
    let target_q_sup = q_sup_envelope(target_product, max_log_scale)?;
    let target_log_derivative_q_sup = log_derivative_q_sup_envelope(target_product, max_log_scale)?;

    let remainder = sine_taylor_remainder_envelope(omega, degree)?;
    let error_product =
        product_derivative_envelope(bump_derivatives.as_array(), remainder.derivative_bounds())?;
    let error_q_sup = q_sup_envelope(error_product, max_log_scale)?;
    let error_log_derivative_q_sup = log_derivative_q_sup_envelope(error_product, max_log_scale)?;

    let target_l2 = sqrt_log_span * target_q_sup;
    let error_l2 = sqrt_log_span * error_q_sup;
    let target_log_derivative_l2 = sqrt_log_span * target_log_derivative_q_sup;
    let derivative_error_l2 = sqrt_log_span * error_log_derivative_q_sup;
    let approximant_l2 = target_l2 + error_l2;
    let approximant_log_derivative_l2 = target_log_derivative_l2 + derivative_error_l2;

    for (stage, value) in [
        ("target L2 bound", target_l2),
        ("approximant L2 bound", approximant_l2),
        ("error L2 bound", error_l2),
        ("target log derivative L2 bound", target_log_derivative_l2),
        (
            "approximant log derivative L2 bound",
            approximant_log_derivative_l2,
        ),
        ("log derivative error L2 bound", derivative_error_l2),
    ] {
        checked_nonnegative(stage, value)?;
    }

    let h1_bounds = H1ApproximationUpperBounds::new(
        target_l2,
        approximant_l2,
        error_l2,
        target_log_derivative_l2,
        approximant_log_derivative_l2,
        derivative_error_l2,
    )?;

    Ok(FiniteSineTaylorH1Envelope {
        mode,
        degree,
        omega,
        log_support_span,
        max_log_scale,
        bump_derivatives,
        remainder,
        target_q_sup,
        target_log_derivative_q_sup,
        error_q_sup,
        error_log_derivative_q_sup,
        h1_bounds,
    })
}

/// Feed two concrete sine Taylor H1 envelopes into the finite Weil continuity bound.
pub fn bound_sine_taylor_pairing_perturbation(
    bump: CompactArchimedeanBump,
    left_mode: usize,
    left_degree: usize,
    right_mode: usize,
    right_degree: usize,
) -> Result<FiniteSineTaylorPairingEnvelope, FiniteSineTaylorH1Error> {
    let left = derive_sine_taylor_h1_envelope(bump, left_mode, left_degree)?;
    let right = derive_sine_taylor_h1_envelope(bump, right_mode, right_degree)?;
    let pairing =
        bound_finite_weil_h1_pairing_perturbation(bump, left.h1_bounds(), right.h1_bounds())?;
    Ok(FiniteSineTaylorPairingEnvelope {
        left,
        right,
        pairing,
    })
}

/// Majorants for the standard bump and its first three `t` derivatives.
///
/// Put `x=1/(t(1-t)) >= 4`.  Every derivative is `exp(-x)` times a polynomial
/// in `x` and `1-2t`.  Replacing `|1-2t|` by one and maximizing each monomial
/// `x^k exp(-x)` on `x>=4` gives the closed bounds below.
pub fn standard_bump_derivative_envelope()
-> Result<StandardBumpDerivativeEnvelope, FiniteSineTaylorH1Error> {
    let m0 = x_power_exp_sup(0)?;
    let m2 = x_power_exp_sup(2)?;
    let m3 = x_power_exp_sup(3)?;
    let m4 = x_power_exp_sup(4)?;
    let m5 = x_power_exp_sup(5)?;
    let m6 = x_power_exp_sup(6)?;

    let value = m0;
    let first = m2;
    let second = 2.0 * m2 + 2.0 * m3 + m4;
    let third = 12.0 * m3 + 12.0 * m4 + 6.0 * m5 + m6;
    for (stage, bound) in [
        ("bump value derivative envelope", value),
        ("bump first derivative envelope", first),
        ("bump second derivative envelope", second),
        ("bump third derivative envelope", third),
    ] {
        checked_nonnegative(stage, bound)?;
    }
    Ok(StandardBumpDerivativeEnvelope {
        value,
        first,
        second,
        third,
    })
}

fn sine_taylor_remainder_envelope(
    omega: f64,
    degree: usize,
) -> Result<SineTaylorRemainderEnvelope, FiniteSineTaylorH1Error> {
    let mut derivative_bounds = [0.0_f64; 4];
    for (derivative_order, derivative_bound) in derivative_bounds.iter_mut().enumerate() {
        *derivative_bound =
            sine_taylor_remainder_derivative_bound(omega, degree, derivative_order)?;
    }
    Ok(SineTaylorRemainderEnvelope { derivative_bounds })
}

/// Taylor remainder for the `r`th derivative on `|t-1/2|<=1/2`:
///
/// `omega^(n+1) (1/2)^(n+1-r) / (n+1-r)!`.
fn sine_taylor_remainder_derivative_bound(
    omega: f64,
    degree: usize,
    derivative_order: usize,
) -> Result<f64, FiniteSineTaylorH1Error> {
    debug_assert!(derivative_order <= 3);
    debug_assert!(degree >= 2);
    let remaining_order = degree + 1 - derivative_order;
    let mut bound = match derivative_order {
        0 => 1.0,
        1 => omega,
        2 => omega * omega,
        3 => omega * omega * omega,
        _ => unreachable!(),
    };
    checked_nonnegative("Taylor derivative frequency factor", bound)?;
    let half_omega = 0.5 * omega;
    for denominator in 1..=remaining_order {
        bound *= half_omega / denominator as f64;
        checked_nonnegative("Taylor derivative remainder", bound)?;
    }
    Ok(bound)
}

/// Product-rule majorants for derivatives zero through three.
fn product_derivative_envelope(
    left: [f64; 4],
    right: [f64; 4],
) -> Result<[f64; 4], FiniteSineTaylorH1Error> {
    let result = [
        left[0] * right[0],
        left[1] * right[0] + left[0] * right[1],
        left[2] * right[0] + 2.0 * left[1] * right[1] + left[0] * right[2],
        left[3] * right[0]
            + 3.0 * left[2] * right[1]
            + 3.0 * left[1] * right[2]
            + left[0] * right[3],
    ];
    for value in result {
        checked_nonnegative("product derivative envelope", value)?;
    }
    Ok(result)
}

/// `Q g = -A g' - A^2 g'' + g/4`, with `A=rho/(b-a)`.
fn q_sup_envelope(
    derivatives: [f64; 4],
    max_log_scale: f64,
) -> Result<f64, FiniteSineTaylorH1Error> {
    let bound = 0.25 * derivatives[0]
        + max_log_scale * derivatives[1]
        + max_log_scale * max_log_scale * derivatives[2];
    checked_nonnegative("Q supremum envelope", bound)?;
    Ok(bound)
}

/// `D(Qg) = -(3/4)A g' - 3 A^2 g'' - A^3 g'''`.
fn log_derivative_q_sup_envelope(
    derivatives: [f64; 4],
    max_log_scale: f64,
) -> Result<f64, FiniteSineTaylorH1Error> {
    let scale2 = max_log_scale * max_log_scale;
    let bound = 0.75 * max_log_scale * derivatives[1]
        + 3.0 * scale2 * derivatives[2]
        + scale2 * max_log_scale * derivatives[3];
    checked_nonnegative("log derivative Q supremum envelope", bound)?;
    Ok(bound)
}

fn x_power_exp_sup(power: usize) -> Result<f64, FiniteSineTaylorH1Error> {
    let maximizer = if power <= 4 { 4.0 } else { power as f64 };
    let bound = maximizer.powi(power as i32) * (-maximizer).exp();
    checked_nonnegative("x^k exp(-x) envelope", bound)?;
    Ok(bound)
}

fn checked_positive(stage: &'static str, value: f64) -> Result<(), FiniteSineTaylorH1Error> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(FiniteSineTaylorH1Error::NonFiniteEvaluation { stage, value })
    }
}

fn checked_nonnegative(stage: &'static str, value: f64) -> Result<(), FiniteSineTaylorH1Error> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(FiniteSineTaylorH1Error::NonFiniteEvaluation { stage, value })
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
    fn bump_derivative_majorants_match_closed_form_construction() {
        let bounds = standard_bump_derivative_envelope().unwrap();
        let exp_minus_four = (-4.0_f64).exp();
        assert!((bounds.value() - exp_minus_four).abs() <= 1.0e-15);
        assert!((bounds.first() - 16.0 * exp_minus_four).abs() <= 1.0e-14);
        assert!((bounds.second() - 416.0 * exp_minus_four).abs() <= 1.0e-12);
        assert!(bounds.third().is_finite());
        assert!(bounds.third() > bounds.second());
    }

    #[test]
    fn centered_taylor_remainders_contract_with_degree_for_first_mode() {
        let low = derive_sine_taylor_h1_envelope(bump(), 1, 6).unwrap();
        let high = derive_sine_taylor_h1_envelope(bump(), 1, 10).unwrap();
        for (high_bound, low_bound) in high
            .remainder()
            .derivative_bounds()
            .into_iter()
            .zip(low.remainder().derivative_bounds())
        {
            assert!(high_bound < low_bound);
        }
        assert!(high.h1_bounds().error_l2() < low.h1_bounds().error_l2());
        assert!(high.h1_bounds().derivative_error_l2() < low.h1_bounds().derivative_error_l2());
    }

    #[test]
    fn target_bounds_do_not_depend_on_taylor_degree() {
        let low = derive_sine_taylor_h1_envelope(bump(), 2, 6).unwrap();
        let high = derive_sine_taylor_h1_envelope(bump(), 2, 12).unwrap();
        assert_eq!(low.h1_bounds().target_l2(), high.h1_bounds().target_l2());
        assert_eq!(
            low.h1_bounds().target_log_derivative_l2(),
            high.h1_bounds().target_log_derivative_l2()
        );
    }

    #[test]
    fn concrete_pairing_envelope_tightens_with_taylor_degree() {
        let low = bound_sine_taylor_pairing_perturbation(bump(), 1, 6, 2, 8).unwrap();
        let high = bound_sine_taylor_pairing_perturbation(bump(), 1, 10, 2, 12).unwrap();
        assert!(
            high.pairing().total_pairing_error_bound() < low.pairing().total_pairing_error_bound()
        );
        assert!(high.pairing().total_pairing_error_bound().is_finite());
    }

    #[test]
    fn invalid_mode_and_degree_fail_closed() {
        assert!(matches!(
            derive_sine_taylor_h1_envelope(bump(), 0, 6),
            Err(FiniteSineTaylorH1Error::ZeroMode)
        ));
        assert!(matches!(
            derive_sine_taylor_h1_envelope(bump(), 1, 1),
            Err(FiniteSineTaylorH1Error::DegreeTooSmall { degree: 1 })
        ));
    }
}
