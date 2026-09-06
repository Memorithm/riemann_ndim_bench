//! Finite-dimensional Riemann--Weil quadratic-form matrix on a compact,
//! boundary-admissible test-function family.
//!
//! The basis is built as
//!
//! `g_j(rho) = bump(rho) P_j(2t-1)`,
//!
//! where `t=(rho-a)/(b-a)` and `P_j` is the Legendre polynomial of degree `j`.
//! Each basis function is then
//!
//! `h_j = Q g_j`, `Q = -(rho d/drho)^2 + 1/4`.
//!
//! Thus every basis vector inherits the same compact support and is audited
//! against the two critical Mellin boundary conditions. Mixed entries evaluate
//! the source pairing `psi(h_i^* * h_j)` using the same compact Riemann--Weil
//! decomposition as `weil_finite_functional`.
//!
//! A finite matrix, even when positive semidefinite, is only a restriction of
//! the Weil quadratic form to the declared finite basis. No finite-dimensional
//! eigenvalue sign is promoted to Weil positivity or RH.

use std::f64::consts::PI;
use std::fmt;

use faer::{Mat, Side, linalg::solvers::SelfAdjointEigen};

use crate::quadrature::{GaussLegendreUnit, QuadratureError};
use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::semilocal_compact_weil::CompactWeilTestFunction;
use crate::weil_boundary::{WeilBoundaryError, WeilBoundaryMoments, critical_boundary_moments};

const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;
const LOG_UNDERFLOW_GUARD: f64 = -700.0;

/// One compact boundary-admissible basis vector `h_j = Q g_j`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactWeilBasisFunction {
    bump: CompactArchimedeanBump,
    degree: usize,
}

impl CompactWeilBasisFunction {
    #[inline]
    pub const fn new(bump: CompactArchimedeanBump, degree: usize) -> Self {
        Self { bump, degree }
    }

    #[inline]
    pub const fn bump(self) -> CompactArchimedeanBump {
        self.bump
    }

    #[inline]
    pub const fn degree(self) -> usize {
        self.degree
    }

    /// Evaluate `h_j(rho)=Q g_j(rho)`, returning exact numerical zero outside
    /// the declared compact support.
    pub fn value(self, rho: f64) -> Result<f64, WeilBoundaryError> {
        if !rho.is_finite() || rho <= 0.0 {
            return Err(WeilBoundaryError::InvalidRho { rho });
        }
        if !self.bump.support().contains(rho) {
            return Ok(0.0);
        }
        Ok(self.value_inside(rho))
    }

    /// Numerically audit the two critical Mellin moments of this basis vector.
    pub fn boundary_moments(
        self,
        quadrature_order: usize,
    ) -> Result<WeilBoundaryMoments, WeilBoundaryError> {
        critical_boundary_moments(self.bump.support(), quadrature_order, |rho| {
            self.value_inside(rho)
        })
    }

    fn value_inside(self, rho: f64) -> f64 {
        let support = self.bump.support();
        let lower = support.lower();
        let upper = support.upper();
        let width = upper - lower;
        let t = (rho - lower) / width;
        if !(0.0 < t && t < 1.0) {
            return 0.0;
        }

        let d = t * (1.0 - t);
        let exponent = -1.0 / d;
        if exponent < LOG_UNDERFLOW_GUARD {
            return 0.0;
        }

        let bump = exponent.exp();
        let first_h = (1.0 - 2.0 * t) / d.powi(2);
        let second_h = -2.0 / d.powi(2) - 2.0 * (1.0 - 2.0 * t).powi(2) / d.powi(3);
        let bump_t = bump * first_h;
        let bump_tt = bump * (second_h + first_h * first_h);

        let y = 2.0 * t - 1.0;
        let (polynomial, derivative_y, second_derivative_y) = legendre_with_derivatives(self.degree, y);
        let polynomial_t = 2.0 * derivative_y;
        let polynomial_tt = 4.0 * second_derivative_y;

        let generator = bump * polynomial;
        let generator_t = bump_t * polynomial + bump * polynomial_t;
        let generator_tt = bump_tt * polynomial
            + 2.0 * bump_t * polynomial_t
            + bump * polynomial_tt;

        // In the affine coordinate t, rho d/drho = A d/dt with A=rho/width,
        // hence (rho d/drho)^2 g = A g_t + A^2 g_tt.
        let scale = rho / width;
        let second_log_derivative = scale * generator_t + scale * scale * generator_tt;
        let value = -second_log_derivative + 0.25 * generator;

        debug_assert!(value.is_finite());
        value
    }
}

/// Source decomposition of one mixed basis pairing `psi(h_i^* * h_j)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteWeilPairingAudit {
    left_degree: usize,
    right_degree: usize,
    pole_term: f64,
    archimedean_term: f64,
    prime_total: f64,
    value: f64,
}

impl FiniteWeilPairingAudit {
    #[inline]
    pub const fn left_degree(self) -> usize {
        self.left_degree
    }

    #[inline]
    pub const fn right_degree(self) -> usize {
        self.right_degree
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

/// Finite Hermitian (real-symmetric here) Weil quadratic-form matrix audit.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilQuadraticMatrixAudit {
    dimension: usize,
    entries: Vec<f64>,
    eigenvalues: Vec<f64>,
    boundary_residuals: Vec<f64>,
    max_raw_pairing_asymmetry: f64,
}

impl FiniteWeilQuadraticMatrixAudit {
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn entry(&self, row: usize, col: usize) -> Option<f64> {
        (row < self.dimension && col < self.dimension)
            .then_some(self.entries[row * self.dimension + col])
    }

    #[inline]
    pub fn eigenvalues(&self) -> &[f64] {
        &self.eigenvalues
    }

    #[inline]
    pub fn minimum_eigenvalue(&self) -> f64 {
        self.eigenvalues[0]
    }

    #[inline]
    pub fn boundary_residuals(&self) -> &[f64] {
        &self.boundary_residuals
    }

    #[inline]
    pub fn max_boundary_residual(&self) -> f64 {
        self.boundary_residuals
            .iter()
            .copied()
            .fold(0.0_f64, f64::max)
    }

    #[inline]
    pub fn max_raw_pairing_asymmetry(&self) -> f64 {
        self.max_raw_pairing_asymmetry
    }

    /// Smallest eigenvalue of the leading `size x size` principal submatrix.
    pub fn principal_minimum_eigenvalue(
        &self,
        size: usize,
    ) -> Result<f64, FiniteWeilMatrixError> {
        if size == 0 || size > self.dimension {
            return Err(FiniteWeilMatrixError::InvalidPrincipalDimension {
                requested: size,
                available: self.dimension,
            });
        }
        let matrix = Mat::from_fn(size, size, |i, j| self.entries[i * self.dimension + j]);
        smallest_self_adjoint_eigenvalue(matrix)
    }
}

#[derive(Debug)]
pub enum FiniteWeilMatrixError {
    EmptyBasis,
    Quadrature(QuadratureError),
    Boundary(WeilBoundaryError),
    SupportRatioTooLarge { floor: u128 },
    DecompositionFailed,
    InvalidPrincipalDimension { requested: usize, available: usize },
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilMatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyBasis => write!(f, "finite Weil quadratic-form basis must be non-empty"),
            Self::Quadrature(error) => write!(f, "quadrature construction failed: {error:?}"),
            Self::Boundary(error) => write!(f, "boundary-admissible basis evaluation failed: {error}"),
            Self::SupportRatioTooLarge { floor } => write!(
                f,
                "compact support ratio requires a prime-power bound larger than u64: floor={floor}"
            ),
            Self::DecompositionFailed => write!(f, "self-adjoint eigendecomposition failed"),
            Self::InvalidPrincipalDimension { requested, available } => write!(
                f,
                "invalid principal matrix dimension {requested}; available dimension is {available}"
            ),
            Self::NonFiniteEvaluation { stage, value } => {
                write!(f, "non-finite Weil matrix value at {stage}: {value}")
            }
        }
    }
}

impl std::error::Error for FiniteWeilMatrixError {}

impl From<QuadratureError> for FiniteWeilMatrixError {
    fn from(value: QuadratureError) -> Self {
        Self::Quadrature(value)
    }
}

impl From<WeilBoundaryError> for FiniteWeilMatrixError {
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
    fn from_bump(bump: CompactArchimedeanBump) -> Result<Self, FiniteWeilMatrixError> {
        let lower = bump.lower();
        let upper = bump.upper();
        let numerator = u128::from(upper.numerator()) * u128::from(lower.denominator());
        let denominator = u128::from(upper.denominator()) * u128::from(lower.numerator());
        let floor = numerator / denominator;
        if floor > u128::from(u64::MAX) {
            return Err(FiniteWeilMatrixError::SupportRatioTooLarge { floor });
        }
        Ok(Self {
            numerator,
            denominator,
            floor: floor as u64,
        })
    }

    fn integer_is_boundary(self, integer: u64) -> bool {
        u128::from(integer) * self.denominator == self.numerator
    }

    fn as_f64(self) -> f64 {
        self.numerator as f64 / self.denominator as f64
    }
}

struct MixedLogCorrelation {
    starred: CompactWeilBasisFunction,
    unstarred: CompactWeilBasisFunction,
    quadrature: GaussLegendreUnit,
    log_lower: f64,
    log_upper: f64,
    log_span: f64,
}

impl MixedLogCorrelation {
    fn new(
        starred: CompactWeilBasisFunction,
        unstarred: CompactWeilBasisFunction,
        quadrature_order: usize,
    ) -> Result<Self, FiniteWeilMatrixError> {
        let support = starred.bump().support();
        let log_lower = support.log_lower();
        let log_upper = support.log_upper();
        Ok(Self {
            starred,
            unstarred,
            quadrature: GaussLegendreUnit::new(quadrature_order)?,
            log_lower,
            log_upper,
            log_span: log_upper - log_lower,
        })
    }

    /// Log-coordinate representation of `(starred^* * unstarred)(exp(shift))`.
    fn value(&self, shift: f64) -> Result<f64, FiniteWeilMatrixError> {
        if !shift.is_finite() {
            return Err(FiniteWeilMatrixError::NonFiniteEvaluation {
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
            // For F = f^* * g one obtains integral f(exp(u)) g(exp(u+t)) du.
            let left = self.starred.value(u.exp())?;
            let right = self.unstarred.value((u + shift).exp())?;
            total += weight * left * right;
        }
        let value = span * total;
        checked_finite("mixed log correlation", value)?;
        Ok(value)
    }
}

/// Build the real-symmetric matrix of `psi(h_i^* * h_j)` for Legendre-weighted
/// compact generators of degrees `0..dimension`.
pub fn audit_finite_weil_quadratic_matrix(
    bump: CompactArchimedeanBump,
    dimension: usize,
    correlation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
) -> Result<FiniteWeilQuadraticMatrixAudit, FiniteWeilMatrixError> {
    if dimension == 0 {
        return Err(FiniteWeilMatrixError::EmptyBasis);
    }

    let ratio = ExactSupportRatio::from_bump(bump)?;
    let basis = (0..dimension)
        .map(|degree| CompactWeilBasisFunction::new(bump, degree))
        .collect::<Vec<_>>();

    let mut boundary_residuals = Vec::with_capacity(dimension);
    for &function in &basis {
        let moments = function.boundary_moments(boundary_order)?;
        boundary_residuals.push(moments.plus_half.abs().max(moments.minus_half.abs()));
    }

    let mut entries = vec![0.0_f64; dimension * dimension];
    let mut max_raw_pairing_asymmetry = 0.0_f64;
    for i in 0..dimension {
        for j in i..dimension {
            let forward = audit_pairing(
                basis[i],
                basis[j],
                ratio,
                correlation_order,
                archimedean_order,
                boundary_order,
            )?;
            let value = if i == j {
                forward.value()
            } else {
                let reverse = audit_pairing(
                    basis[j],
                    basis[i],
                    ratio,
                    correlation_order,
                    archimedean_order,
                    boundary_order,
                )?;
                max_raw_pairing_asymmetry = max_raw_pairing_asymmetry
                    .max((forward.value() - reverse.value()).abs());
                0.5 * (forward.value() + reverse.value())
            };
            entries[i * dimension + j] = value;
            entries[j * dimension + i] = value;
        }
    }

    let matrix = Mat::from_fn(dimension, dimension, |i, j| entries[i * dimension + j]);
    let decomposition = SelfAdjointEigen::new(matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilMatrixError::DecompositionFailed)?;
    let diagonal = decomposition.S().column_vector();
    let mut eigenvalues = (0..dimension)
        .map(|index| diagonal[index])
        .collect::<Vec<_>>();
    eigenvalues.sort_by(f64::total_cmp);

    Ok(FiniteWeilQuadraticMatrixAudit {
        dimension,
        entries,
        eigenvalues,
        boundary_residuals,
        max_raw_pairing_asymmetry,
    })
}

fn audit_pairing(
    left: CompactWeilBasisFunction,
    right: CompactWeilBasisFunction,
    ratio: ExactSupportRatio,
    correlation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
) -> Result<FiniteWeilPairingAudit, FiniteWeilMatrixError> {
    let correlation = MixedLogCorrelation::new(left, right, correlation_order)?;
    let left_moments = left.boundary_moments(boundary_order)?;
    let right_moments = right.boundary_moments(boundary_order)?;
    let pole_term = left_moments.minus_half * right_moments.plus_half
        + left_moments.plus_half * right_moments.minus_half;
    checked_finite("mixed critical pole term", pole_term)?;

    let theta_zero = correlation.value(0.0)?;
    let archimedean_term = mixed_archimedean_term(
        &correlation,
        theta_zero,
        ratio,
        archimedean_order,
    )?;

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
            (prime as f64).ln() / (integer as f64).sqrt()
                * (theta_positive + theta_negative)
        };
        checked_finite("mixed prime-power contribution", contribution)?;
        prime_total += contribution;
    }
    checked_finite("mixed prime-power total", prime_total)?;

    let value = pole_term - archimedean_term - prime_total;
    checked_finite("mixed Weil pairing", value)?;
    Ok(FiniteWeilPairingAudit {
        left_degree: left.degree(),
        right_degree: right.degree(),
        pole_term,
        archimedean_term,
        prime_total,
        value,
    })
}

fn mixed_archimedean_term(
    correlation: &MixedLogCorrelation,
    theta_zero: f64,
    ratio: ExactSupportRatio,
    quadrature_order: usize,
) -> Result<f64, FiniteWeilMatrixError> {
    let quadrature = GaussLegendreUnit::new(quadrature_order)?;
    let theta_sym_zero = 2.0 * theta_zero;
    let ratio_f64 = ratio.as_f64();
    let coefficient = EULER_MASCHERONI
        + (4.0 * PI * (ratio_f64 - 1.0) / (ratio_f64 + 1.0)).ln();
    checked_finite("mixed archimedean coefficient", coefficient)?;

    let mut weighted_sum = 0.0_f64;
    for (&node, &weight) in quadrature.nodes().iter().zip(quadrature.weights().iter()) {
        let t = correlation.log_span * node;
        let theta_symmetric = correlation.value(t)? + correlation.value(-t)?;
        let exp_half = (0.5 * t).exp();
        let numerator = (0.5 * t).exp_m1() * theta_sym_zero
            + exp_half * (theta_symmetric - theta_sym_zero);
        let denominator = 2.0 * t.sinh();
        let integrand = numerator / denominator;
        checked_finite("mixed archimedean integrand", integrand)?;
        weighted_sum += weight * integrand;
    }

    let value = 0.5 * theta_sym_zero * coefficient + correlation.log_span * weighted_sum;
    checked_finite("mixed archimedean term", value)?;
    Ok(value)
}

fn smallest_self_adjoint_eigenvalue(matrix: Mat<f64>) -> Result<f64, FiniteWeilMatrixError> {
    let dimension = matrix.nrows();
    let decomposition = SelfAdjointEigen::new(matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilMatrixError::DecompositionFailed)?;
    let diagonal = decomposition.S().column_vector();
    (0..dimension)
        .map(|index| diagonal[index])
        .min_by(f64::total_cmp)
        .ok_or(FiniteWeilMatrixError::EmptyBasis)
}

fn legendre_with_derivatives(degree: usize, x: f64) -> (f64, f64, f64) {
    if degree == 0 {
        return (1.0, 0.0, 0.0);
    }
    if degree == 1 {
        return (x, 1.0, 0.0);
    }

    let mut p_nm2 = 1.0;
    let mut p_nm1 = x;
    let mut d_nm2 = 0.0;
    let mut d_nm1 = 1.0;
    let mut dd_nm2 = 0.0;
    let mut dd_nm1 = 0.0;

    for n in 2..=degree {
        let n_f = n as f64;
        let a = (2 * n - 1) as f64;
        let b = (n - 1) as f64;
        let p = (a * x * p_nm1 - b * p_nm2) / n_f;
        let d = (a * (p_nm1 + x * d_nm1) - b * d_nm2) / n_f;
        let dd = (a * (2.0 * d_nm1 + x * dd_nm1) - b * dd_nm2) / n_f;
        p_nm2 = p_nm1;
        p_nm1 = p;
        d_nm2 = d_nm1;
        d_nm1 = d;
        dd_nm2 = dd_nm1;
        dd_nm1 = dd;
    }

    (p_nm1, d_nm1, dd_nm1)
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

fn checked_finite(stage: &'static str, value: f64) -> Result<(), FiniteWeilMatrixError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilMatrixError::NonFiniteEvaluation { stage, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::PositiveRational;
    use crate::weil_finite_functional::audit_finite_weil_functional;

    fn bump() -> CompactArchimedeanBump {
        CompactArchimedeanBump::new(
            PositiveRational::new(1, 2).unwrap(),
            PositiveRational::new(7, 2).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn degree_zero_basis_matches_existing_compact_q_profile() {
        let bump = bump();
        let basis = CompactWeilBasisFunction::new(bump, 0);
        let scalar = CompactWeilTestFunction::new(bump);
        for rho in [0.75_f64, 1.0, 1.5, 2.0, 3.0] {
            let left = basis.value(rho).unwrap();
            let right = scalar.q_value(rho).unwrap();
            assert!((left - right).abs() <= 2.0e-14 * left.abs().max(right.abs()).max(1.0));
        }
    }

    #[test]
    fn one_dimensional_matrix_reproduces_scalar_weil_audit() {
        let bump = bump();
        let matrix = audit_finite_weil_quadratic_matrix(bump, 1, 96, 96, 128).unwrap();
        let scalar = audit_finite_weil_functional(
            CompactWeilTestFunction::new(bump),
            96,
            96,
            128,
        )
        .unwrap();

        assert_eq!(matrix.dimension(), 1);
        assert!((matrix.entry(0, 0).unwrap() - scalar.functional_value()).abs() <= 2.0e-10);
        assert!((matrix.minimum_eigenvalue() - scalar.functional_value()).abs() <= 2.0e-10);
    }

    #[test]
    fn finite_pairing_matrix_is_numerically_symmetric_and_boundary_admissible() {
        let audit = audit_finite_weil_quadratic_matrix(bump(), 3, 96, 96, 128).unwrap();
        assert_eq!(audit.dimension(), 3);
        assert_eq!(audit.eigenvalues().len(), 3);
        assert!(audit.max_raw_pairing_asymmetry() <= 5.0e-12);
        assert!(audit.max_boundary_residual() <= 2.0e-10);
        assert!(audit.eigenvalues().iter().all(|value| value.is_finite()));
    }

    #[test]
    fn rayleigh_ritz_minimum_is_nonincreasing_with_basis_size() {
        let audit = audit_finite_weil_quadratic_matrix(bump(), 4, 96, 96, 128).unwrap();
        let lambda1 = audit.principal_minimum_eigenvalue(1).unwrap();
        let lambda2 = audit.principal_minimum_eigenvalue(2).unwrap();
        let lambda3 = audit.principal_minimum_eigenvalue(3).unwrap();
        let lambda4 = audit.principal_minimum_eigenvalue(4).unwrap();
        let tolerance = 2.0e-10;
        assert!(lambda2 <= lambda1 + tolerance);
        assert!(lambda3 <= lambda2 + tolerance);
        assert!(lambda4 <= lambda3 + tolerance);
    }

    #[test]
    fn minimum_eigenvalue_stabilizes_under_quadrature_refinement() {
        let coarse = audit_finite_weil_quadratic_matrix(bump(), 2, 48, 48, 64).unwrap();
        let medium = audit_finite_weil_quadratic_matrix(bump(), 2, 72, 72, 96).unwrap();
        let fine = audit_finite_weil_quadratic_matrix(bump(), 2, 96, 96, 128).unwrap();

        let coarse_error = (coarse.minimum_eigenvalue() - fine.minimum_eigenvalue()).abs();
        let medium_error = (medium.minimum_eigenvalue() - fine.minimum_eigenvalue()).abs();
        assert!(medium_error < coarse_error);
        assert!(medium_error <= 2.0e-6);
    }
}
