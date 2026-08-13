//! Core mathematical primitives for the Riemann N-dimensional research bench.
//!
//! This crate deliberately separates exact identities from experimental
//! geometric interpretations. At this stage we do not evaluate zeta(s) and we
//! do not claim to test or prove the Riemann hypothesis.

use std::f64::consts::PI;

/// Real part of the critical line.
pub const CRITICAL_SIGMA: f64 = 0.5;

/// A point `s = sigma + i t` in the complex plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralPoint {
    pub sigma: f64,
    pub t: f64,
}

impl SpectralPoint {
    pub const fn new(sigma: f64, t: f64) -> Self {
        Self { sigma, t }
    }

    /// The map appearing in the completed zeta functional equation:
    /// `s -> 1 - s`.
    pub fn functional_reflection(self) -> Self {
        Self::new(1.0 - self.sigma, -self.t)
    }

    /// Complex conjugation: `s -> conjugate(s)`.
    pub fn conjugate(self) -> Self {
        Self::new(self.sigma, -self.t)
    }

    /// Geometric reflection across the critical line `Re(s) = 1/2`:
    /// `s -> 1 - conjugate(s)`.
    pub fn critical_line_reflection(self) -> Self {
        Self::new(1.0 - self.sigma, self.t)
    }

    /// Signed Euclidean displacement from the critical line in the real
    /// direction.
    pub fn critical_displacement(self) -> f64 {
        self.sigma - CRITICAL_SIGMA
    }
}

/// Geometry induced by the exact factor `pi^(-s/2)` occurring in the
/// completed zeta function.
///
/// The formulas below are exact consequences of choosing the modulus
/// `R_pi(sigma) = pi^(-sigma/2)` as a radial coordinate. Interpreting that
/// coordinate as a physical/geometric radius is experimental and is not a
/// theorem about zeta(s).
pub struct PiRadialGeometry;

impl PiRadialGeometry {
    /// Raw radial scale `pi^(-sigma/2)`.
    pub fn raw_radius(sigma: f64) -> f64 {
        PI.powf(-0.5 * sigma)
    }

    /// Raw radius at the critical line: `pi^(-1/4)`.
    pub fn critical_radius() -> f64 {
        PI.powf(-0.25)
    }

    /// Radius normalized so that the critical line has radius 1:
    ///
    /// `rho(sigma) = R_pi(sigma) / R_pi(1/2)
    ///             = pi^((1 - 2 sigma)/4)`.
    pub fn normalized_radius(sigma: f64) -> f64 {
        PI.powf((1.0 - 2.0 * sigma) / 4.0)
    }

    /// Logarithmic radial coordinate.
    ///
    /// This is antisymmetric under `sigma -> 1 - sigma`.
    pub fn log_normalized_radius(sigma: f64) -> f64 {
        ((1.0 - 2.0 * sigma) / 4.0) * PI.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1.0e-12;

    fn approx_eq(a: f64, b: f64) {
        assert!((a - b).abs() <= EPS, "left={a:.16e}, right={b:.16e}");
    }

    #[test]
    fn functional_reflection_is_an_involution() {
        let s = SpectralPoint::new(0.37, 14.0);
        assert_eq!(s.functional_reflection().functional_reflection(), s);
    }

    #[test]
    fn critical_line_reflection_is_an_involution() {
        let s = SpectralPoint::new(0.37, 14.0);
        assert_eq!(s.critical_line_reflection().critical_line_reflection(), s);
    }

    #[test]
    fn functional_reflection_equals_conjugation_after_line_reflection() {
        let s = SpectralPoint::new(0.37, 14.0);
        assert_eq!(
            s.critical_line_reflection().conjugate(),
            s.functional_reflection()
        );
    }

    #[test]
    fn critical_line_is_fixed_by_geometric_reflection() {
        let s = SpectralPoint::new(CRITICAL_SIGMA, 42.0);
        assert_eq!(s.critical_line_reflection(), s);
    }

    #[test]
    fn normalized_radius_is_one_on_critical_line() {
        approx_eq(PiRadialGeometry::normalized_radius(CRITICAL_SIGMA), 1.0);
    }

    #[test]
    fn normalized_radii_are_reciprocal_under_critical_reflection() {
        for sigma in [0.1, 0.3, 0.49, 0.5, 0.51, 0.7, 0.9] {
            let rho = PiRadialGeometry::normalized_radius(sigma);
            let rho_reflected = PiRadialGeometry::normalized_radius(1.0 - sigma);
            approx_eq(rho * rho_reflected, 1.0);
        }
    }

    #[test]
    fn logarithmic_radius_changes_sign_under_reflection() {
        for sigma in [0.1, 0.3, 0.49, 0.51, 0.7, 0.9] {
            let q = PiRadialGeometry::log_normalized_radius(sigma);
            let q_reflected = PiRadialGeometry::log_normalized_radius(1.0 - sigma);
            approx_eq(q_reflected, -q);
        }
    }

    #[test]
    fn raw_critical_radius_matches_pi_to_minus_one_quarter() {
        approx_eq(
            PiRadialGeometry::critical_radius(),
            PI.powf(-0.25),
        );
    }
}
