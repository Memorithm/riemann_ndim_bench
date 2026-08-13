use std::f64::consts::PI;

pub const CRITICAL_SIGMA: f64 = 0.5;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralPoint {
    pub sigma: f64,
    pub t: f64,
}

impl SpectralPoint {
    pub const fn new(sigma: f64, t: f64) -> Self {
        Self { sigma, t }
    }

    pub fn functional_reflection(self) -> Self {
        Self::new(1.0 - self.sigma, -self.t)
    }

    pub fn conjugate(self) -> Self {
        Self::new(self.sigma, -self.t)
    }

    pub fn critical_line_reflection(self) -> Self {
        Self::new(1.0 - self.sigma, self.t)
    }

    pub fn critical_displacement(self) -> f64 {
        self.sigma - CRITICAL_SIGMA
    }
}

pub struct PiRadialGeometry;

impl PiRadialGeometry {
    pub fn raw_radius(sigma: f64) -> f64 {
        PI.powf(-0.5 * sigma)
    }

    pub fn critical_radius() -> f64 {
        PI.powf(-0.25)
    }

    pub fn normalized_radius(sigma: f64) -> f64 {
        PI.powf((1.0 - 2.0 * sigma) / 4.0)
    }

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
            let reflected = PiRadialGeometry::normalized_radius(1.0 - sigma);
            approx_eq(rho * reflected, 1.0);
        }
    }

    #[test]
    fn logarithmic_radius_changes_sign_under_reflection() {
        for sigma in [0.1, 0.3, 0.49, 0.51, 0.7, 0.9] {
            let q = PiRadialGeometry::log_normalized_radius(sigma);
            let reflected = PiRadialGeometry::log_normalized_radius(1.0 - sigma);
            approx_eq(reflected, -q);
        }
    }

    #[test]
    fn raw_critical_radius_matches_pi_to_minus_one_quarter() {
        approx_eq(PiRadialGeometry::critical_radius(), PI.powf(-0.25));
    }
}
