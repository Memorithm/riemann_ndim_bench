use num_complex::Complex64;
use std::f64::consts::PI;

const BERN_COEFF: [f64; 10] = [
    1.0 / 12.0,
    -1.0 / 720.0,
    1.0 / 30_240.0,
    -1.0 / 1_209_600.0,
    1.0 / 47_900_160.0,
    -691.0 / 1_307_674_368_000.0,
    1.0 / 74_724_249_600.0,
    -3617.0 / 10_670_622_842_880_000.0,
    43_867.0 / 5_109_094_217_170_944_000.0,
    -174_611.0 / 802_857_662_698_291_200_000.0,
];

const LANCZOS_G: f64 = 7.0;
const SQRT_2_PI: f64 = 2.506_628_274_631_000_7;
const LANCZOS_COEFF: [f64; 9] = [
    0.999_999_999_999_809_9,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_1,
    -176.615_029_162_140_6,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NumericError {
    PoleAtOne,
    InvalidConfiguration,
}

#[derive(Debug, Clone, Copy)]
pub struct EulerMaclaurinConfig {
    pub n: usize,
    pub bernoulli_terms: usize,
}

impl EulerMaclaurinConfig {
    pub const fn new(n: usize, bernoulli_terms: usize) -> Self {
        Self { n, bernoulli_terms }
    }

    pub fn for_point(s: Complex64) -> Self {
        let n = ((1.5 * s.im.abs()).ceil() as usize + 16).max(32);
        Self::new(n, 8)
    }

    fn validate(self) -> Result<Self, NumericError> {
        if self.n < 2 || self.bernoulli_terms > BERN_COEFF.len() {
            return Err(NumericError::InvalidConfiguration);
        }
        Ok(self)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ZetaEstimate {
    pub value: Complex64,
    pub cross_resolution_delta: f64,
    pub fine_n: usize,
    pub bernoulli_terms: usize,
}

fn real_base_complex_power(base: f64, exponent: Complex64) -> Complex64 {
    (exponent * base.ln()).exp()
}

fn pochhammer(s: Complex64, order: usize) -> Complex64 {
    let mut product = Complex64::new(1.0, 0.0);
    for j in 0..order {
        product *= s + Complex64::new(j as f64, 0.0);
    }
    product
}

fn correction_term(s: Complex64, n: usize, k: usize) -> Complex64 {
    let exponent = Complex64::new(1.0 - 2.0 * k as f64, 0.0) - s;
    let n_power = real_base_complex_power(n as f64, exponent);
    pochhammer(s, 2 * k - 1) * n_power * BERN_COEFF[k - 1]
}

/// Euler-Maclaurin evaluation of the analytically continued Riemann zeta
/// function. The finite correction follows the standard Bernoulli expansion.
pub fn zeta_euler_maclaurin(
    s: Complex64,
    config: EulerMaclaurinConfig,
) -> Result<Complex64, NumericError> {
    let config = config.validate()?;
    if s == Complex64::new(1.0, 0.0) {
        return Err(NumericError::PoleAtOne);
    }

    let mut sum = Complex64::new(0.0, 0.0);
    for n in 1..config.n {
        sum += real_base_complex_power(n as f64, -s);
    }

    let n = config.n as f64;
    let tail_integral =
        real_base_complex_power(n, Complex64::new(1.0, 0.0) - s) / (s - Complex64::new(1.0, 0.0));
    let endpoint = real_base_complex_power(n, -s) * 0.5;

    sum += tail_integral + endpoint;

    for k in 1..=config.bernoulli_terms {
        sum += correction_term(s, config.n, k);
    }

    Ok(sum)
}

pub fn zeta_checked(s: Complex64) -> Result<ZetaEstimate, NumericError> {
    let coarse_config = EulerMaclaurinConfig::for_point(s);
    let fine_config = EulerMaclaurinConfig::new(coarse_config.n * 2, BERN_COEFF.len());

    let coarse_value = zeta_euler_maclaurin(s, coarse_config)?;
    let value = zeta_euler_maclaurin(s, fine_config)?;
    let cross_resolution_delta = (value - coarse_value).norm();

    Ok(ZetaEstimate {
        value,
        cross_resolution_delta,
        fine_n: fine_config.n,
        bernoulli_terms: fine_config.bernoulli_terms,
    })
}

/// Complex Gamma using the Lanczos approximation and the reflection formula.
pub fn gamma_lanczos(z: Complex64) -> Complex64 {
    if z.re < 0.5 {
        let pi_z = z * PI;
        let denominator = pi_z.sin() * gamma_lanczos(Complex64::new(1.0, 0.0) - z);
        return Complex64::new(PI, 0.0) / denominator;
    }

    let shifted = z - Complex64::new(1.0, 0.0);
    let mut x = Complex64::new(LANCZOS_COEFF[0], 0.0);

    for (index, coefficient) in LANCZOS_COEFF.iter().enumerate().skip(1) {
        x += Complex64::new(*coefficient, 0.0) / (shifted + Complex64::new(index as f64, 0.0));
    }

    let t = shifted + Complex64::new(LANCZOS_G + 0.5, 0.0);
    let exponent = (shifted + Complex64::new(0.5, 0.0)) * t.ln() - t;
    x * exponent.exp() * SQRT_2_PI
}

pub fn xi(s: Complex64) -> Result<Complex64, NumericError> {
    let zeta = zeta_checked(s)?.value;
    let gamma = gamma_lanczos(s * 0.5);
    let pi_factor = (-s * (0.5 * PI.ln())).exp();
    Ok(s * (s - Complex64::new(1.0, 0.0)) * 0.5 * gamma * pi_factor * zeta)
}

pub fn xi_symmetry_residual(s: Complex64) -> Result<f64, NumericError> {
    let left = xi(s)?;
    let right = xi(Complex64::new(1.0, 0.0) - s)?;
    let scale = left.norm().max(right.norm()).max(1.0);
    Ok((left - right).norm() / scale)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zeta_two_matches_pi_squared_over_six() {
        let s = Complex64::new(2.0, 0.0);
        let value = zeta_checked(s).unwrap().value;
        let expected = PI * PI / 6.0;
        assert!((value.re - expected).abs() < 2.0e-13);
        assert!(value.im.abs() < 2.0e-13);
    }

    #[test]
    fn zeta_zero_matches_minus_one_half() {
        let s = Complex64::new(0.0, 0.0);
        let value = zeta_checked(s).unwrap().value;
        assert!((value.re + 0.5).abs() < 2.0e-13);
        assert!(value.im.abs() < 2.0e-13);
    }

    #[test]
    fn gamma_one_is_one() {
        let value = gamma_lanczos(Complex64::new(1.0, 0.0));
        assert!((value.re - 1.0).abs() < 2.0e-14);
        assert!(value.im.abs() < 2.0e-14);
    }

    #[test]
    fn gamma_half_is_sqrt_pi() {
        let value = gamma_lanczos(Complex64::new(0.5, 0.0));
        assert!((value.re - PI.sqrt()).abs() < 2.0e-14);
        assert!(value.im.abs() < 2.0e-14);
    }

    #[test]
    fn zeta_respects_complex_conjugation() {
        let s = Complex64::new(0.37, 9.25);
        let z = zeta_checked(s).unwrap().value;
        let z_conj = zeta_checked(s.conj()).unwrap().value;
        assert!((z_conj - z.conj()).norm() < 2.0e-12);
    }

    #[test]
    fn xi_respects_functional_symmetry() {
        let s = Complex64::new(0.37, 9.25);
        let residual = xi_symmetry_residual(s).unwrap();
        assert!(residual < 2.0e-10, "residual={residual:.3e}");
    }

    #[test]
    fn checked_zeta_converges_across_resolution() {
        let s = Complex64::new(0.5, 14.0);
        let estimate = zeta_checked(s).unwrap();
        assert!(
            estimate.cross_resolution_delta < 1.0e-10,
            "delta={:.3e}",
            estimate.cross_resolution_delta
        );
    }
}
