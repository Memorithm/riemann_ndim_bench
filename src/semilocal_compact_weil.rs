//! Numerical Weil-boundary audit for the exact-support compact archimedean
//! fixture.
//!
//! Connes--Consani use
//!
//! `Q = -(rho d/drho)^2 + 1/4`
//!
//! to impose the two critical Mellin boundary conditions while preserving
//! compact support. This module applies that operator to the manufactured smooth
//! compact bump from `semilocal_compact_archimedean` using an analytic second
//! logarithmic derivative, then audits the two critical moments with the
//! existing source-locked quadrature layer.
//!
//! This is a test-function validation layer. Small numerical boundary moments
//! do not establish Weil positivity, the semilocal trace formula, Conjecture
//! 4.1, or RH.

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_boundary::{
    WeilBoundaryError, WeilBoundaryMoments, critical_boundary_moments,
    q_from_log_second_derivative, q_on_support,
};

const LOG_UNDERFLOW_GUARD: f64 = -700.0;

/// Source-boundary wrapper around one compact archimedean generator `g`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactWeilTestFunction {
    generator: CompactArchimedeanBump,
}

impl CompactWeilTestFunction {
    #[inline]
    pub const fn new(generator: CompactArchimedeanBump) -> Self {
        Self { generator }
    }

    #[inline]
    pub const fn generator(self) -> CompactArchimedeanBump {
        self.generator
    }

    /// Numerically evaluate the compact generator `g(rho)` using the same
    /// profile as `CompactArchimedeanBump` but for an arbitrary positive `rho`.
    pub fn generator_value(self, rho: f64) -> Result<f64, WeilBoundaryError> {
        if !rho.is_finite() || rho <= 0.0 {
            return Err(WeilBoundaryError::InvalidRho { rho });
        }
        if !self.generator.support().contains(rho) {
            return Ok(0.0);
        }
        Ok(self.generator_value_and_second_log_derivative_inside(rho).0)
    }

    /// Evaluate `Qg(rho)` with
    ///
    /// `Q = -(rho d/drho)^2 + 1/4`.
    ///
    /// The second logarithmic derivative is analytic; finite differences are
    /// not used in this audit.
    pub fn q_value(self, rho: f64) -> Result<f64, WeilBoundaryError> {
        q_on_support(self.generator.support(), rho, |rho| {
            self.generator_value_and_second_log_derivative_inside(rho)
        })
    }

    /// Audit the two critical Mellin moments of `Qg` using the existing
    /// log-coordinate Gauss--Legendre quadrature.
    pub fn audit_boundary_moments(
        self,
        quadrature_order: usize,
    ) -> Result<CompactWeilBoundaryAudit, WeilBoundaryError> {
        let moments =
            critical_boundary_moments(self.generator.support(), quadrature_order, |rho| {
                self.q_value_inside(rho)
            })?;
        Ok(CompactWeilBoundaryAudit {
            quadrature_order,
            moments,
        })
    }

    fn q_value_inside(self, rho: f64) -> f64 {
        let (value, second_log_derivative) =
            self.generator_value_and_second_log_derivative_inside(rho);
        q_from_log_second_derivative(value, second_log_derivative)
    }

    /// Return `(g, d^2 g / d(log rho)^2)` for an interior support point.
    ///
    /// Write `t=(rho-a)/(b-a)`, `A=rho/(b-a)` and
    /// `g=exp(h)`, `h=-1/(t(1-t))`. Since
    ///
    /// `rho d/drho = A d/dt`,
    ///
    /// one has
    ///
    /// `d_log^2 g = A g_t + A^2 g_tt`.
    fn generator_value_and_second_log_derivative_inside(self, rho: f64) -> (f64, f64) {
        let support = self.generator.support();
        let lower = support.lower();
        let upper = support.upper();
        let width = upper - lower;
        let t = (rho - lower) / width;

        if !(0.0 < t && t < 1.0) {
            return (0.0, 0.0);
        }

        let d = t * (1.0 - t);
        let exponent = -1.0 / d;
        if exponent < LOG_UNDERFLOW_GUARD {
            // The C-infinity bump and all derivatives vanish super-polynomially
            // at the boundary. Returning numerical zero here prevents a
            // `0 * infinity` indeterminate form after exponential underflow.
            return (0.0, 0.0);
        }

        let value = exponent.exp();
        let first_h = (1.0 - 2.0 * t) / d.powi(2);
        let second_h = -2.0 / d.powi(2) - 2.0 * (1.0 - 2.0 * t).powi(2) / d.powi(3);
        let first_t = value * first_h;
        let second_t = value * (second_h + first_h * first_h);
        let scale = rho / width;
        let second_log_derivative = scale * first_t + scale * scale * second_t;

        debug_assert!(value.is_finite());
        debug_assert!(second_log_derivative.is_finite());
        (value, second_log_derivative)
    }
}

/// Numerical boundary-condition audit for `Qg`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactWeilBoundaryAudit {
    quadrature_order: usize,
    moments: WeilBoundaryMoments,
}

impl CompactWeilBoundaryAudit {
    #[inline]
    pub const fn quadrature_order(self) -> usize {
        self.quadrature_order
    }

    #[inline]
    pub const fn moments(self) -> WeilBoundaryMoments {
        self.moments
    }

    #[inline]
    pub fn max_abs_moment(self) -> f64 {
        self.moments
            .plus_half
            .abs()
            .max(self.moments.minus_half.abs())
    }

    #[inline]
    pub fn satisfies(self, absolute_tolerance: f64) -> bool {
        absolute_tolerance.is_finite()
            && absolute_tolerance >= 0.0
            && self.max_abs_moment() <= absolute_tolerance
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::PositiveRational;

    fn fixture() -> (CompactArchimedeanBump, CompactWeilTestFunction) {
        let bump = CompactArchimedeanBump::new(
            PositiveRational::new(1, 2).unwrap(),
            PositiveRational::new(7, 2).unwrap(),
        )
        .unwrap();
        (bump, CompactWeilTestFunction::new(bump))
    }

    #[test]
    fn arbitrary_rho_generator_matches_discrete_compact_profile() {
        let (bump, test_function) = fixture();
        let x = PositiveRational::new(1, 2).unwrap();

        for m in 2..=6 {
            let rho = m as f64 * x.as_f64();
            let direct = test_function.generator_value(rho).unwrap();
            let discrete = bump.value_at_scaled_m(m, x);
            assert!((direct - discrete).abs() <= 2.0e-15 * direct.abs().max(1.0));
        }
    }

    #[test]
    fn q_profile_preserves_compact_support() {
        let (bump, test_function) = fixture();
        let support = bump.support();

        assert_eq!(test_function.q_value(support.lower()).unwrap(), 0.0);
        assert_eq!(test_function.q_value(support.upper()).unwrap(), 0.0);
        assert_eq!(test_function.q_value(0.5 * support.lower()).unwrap(), 0.0);
        assert_eq!(test_function.q_value(2.0 * support.upper()).unwrap(), 0.0);
        assert!(
            test_function
                .q_value(0.5 * (support.lower() + support.upper()))
                .unwrap()
                .is_finite()
        );
    }

    #[test]
    fn q_image_annihilates_the_two_critical_mellin_moments_numerically() {
        let (_, test_function) = fixture();
        let audit = test_function.audit_boundary_moments(128).unwrap();

        assert_eq!(audit.quadrature_order(), 128);
        assert!(
            audit.satisfies(5.0e-12),
            "plus={:.3e} minus={:.3e}",
            audit.moments().plus_half,
            audit.moments().minus_half
        );
    }

    #[test]
    fn boundary_moment_residual_converges_under_quadrature_refinement() {
        let (_, test_function) = fixture();
        let coarse = test_function.audit_boundary_moments(32).unwrap();
        let medium = test_function.audit_boundary_moments(64).unwrap();
        let fine = test_function.audit_boundary_moments(128).unwrap();

        assert!(medium.max_abs_moment() < coarse.max_abs_moment());
        assert!(fine.max_abs_moment() < medium.max_abs_moment());
        assert!(fine.satisfies(5.0e-12));
    }

    #[test]
    fn invalid_rho_is_rejected_by_the_source_boundary_contract() {
        let (_, test_function) = fixture();
        assert!(matches!(
            test_function.q_value(0.0),
            Err(WeilBoundaryError::InvalidRho { rho }) if rho == 0.0
        ));
        assert!(matches!(
            test_function.q_value(f64::NAN),
            Err(WeilBoundaryError::InvalidRho { rho }) if rho.is_nan()
        ));
    }
}
