//! Localized Bernstein probes for the finite Riemann--Weil quadratic form.
//!
//! For a fixed polynomial degree `D`, the Bernstein enrichment
//!
//! `B_{D,k}(t) = binom(D,k) t^k (1-t)^(D-k)`
//!
//! defines the compact generator `g_{D,k}(rho) = bump(rho) B_{D,k}(t)` and
//! the boundary-admissible direction `h_{D,k} = Q g_{D,k}`.  Because each
//! Bernstein polynomial has degree `D`, it is expanded in the already-audited
//! shifted Legendre parent basis and its one-dimensional Weil/Gram quantities
//! are obtained from the same parent matrices `A` and `G`.
//!
//! A complete Bernstein family of degree `D` spans the same polynomial space
//! as Legendre degrees `0..D`; this module therefore probes different localized
//! directions, not a new complete function space.  Finite positive Rayleigh
//! quotients are numerical evidence only and do not establish complete-space
//! Weil positivity, density/completeness, Conjecture 4.1, or RH.

use std::fmt;

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_boundary::WeilBoundaryError;
use crate::weil_generalized_spectrum::{
    FiniteWeilGeneralizedSpectrumError, audit_finite_weil_generalized_spectrum,
};
use crate::weil_quadratic_matrix::CompactWeilBasisFunction;

const RECONSTRUCTION_INTERVALS: usize = 128;

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilBernsteinProbe {
    index: usize,
    legendre_coefficients: Vec<f64>,
    coefficient_l1_norm: f64,
    raw_quadratic_value: f64,
    gram_norm_squared: f64,
    generalized_rayleigh_quotient: f64,
    boundary_plus_half: f64,
    boundary_minus_half: f64,
    max_boundary_residual: f64,
    max_reconstruction_residual: f64,
}

impl FiniteWeilBernsteinProbe {
    #[inline]
    pub fn index(&self) -> usize {
        self.index
    }

    #[inline]
    pub fn legendre_coefficients(&self) -> &[f64] {
        &self.legendre_coefficients
    }

    #[inline]
    pub fn coefficient_l1_norm(&self) -> f64 {
        self.coefficient_l1_norm
    }

    #[inline]
    pub fn raw_quadratic_value(&self) -> f64 {
        self.raw_quadratic_value
    }

    #[inline]
    pub fn gram_norm_squared(&self) -> f64 {
        self.gram_norm_squared
    }

    #[inline]
    pub fn generalized_rayleigh_quotient(&self) -> f64 {
        self.generalized_rayleigh_quotient
    }

    #[inline]
    pub fn boundary_plus_half(&self) -> f64 {
        self.boundary_plus_half
    }

    #[inline]
    pub fn boundary_minus_half(&self) -> f64 {
        self.boundary_minus_half
    }

    #[inline]
    pub fn max_boundary_residual(&self) -> f64 {
        self.max_boundary_residual
    }

    #[inline]
    pub fn max_reconstruction_residual(&self) -> f64 {
        self.max_reconstruction_residual
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilBernsteinProbeAudit {
    degree: usize,
    parent_dimension: usize,
    parent_gram_condition_number: f64,
    parent_max_raw_pairing_asymmetry: f64,
    parent_max_whitened_asymmetry: f64,
    probes: Vec<FiniteWeilBernsteinProbe>,
}

impl FiniteWeilBernsteinProbeAudit {
    #[inline]
    pub fn degree(&self) -> usize {
        self.degree
    }

    #[inline]
    pub fn parent_dimension(&self) -> usize {
        self.parent_dimension
    }

    #[inline]
    pub fn parent_gram_condition_number(&self) -> f64 {
        self.parent_gram_condition_number
    }

    #[inline]
    pub fn parent_max_raw_pairing_asymmetry(&self) -> f64 {
        self.parent_max_raw_pairing_asymmetry
    }

    #[inline]
    pub fn parent_max_whitened_asymmetry(&self) -> f64 {
        self.parent_max_whitened_asymmetry
    }

    #[inline]
    pub fn probes(&self) -> &[FiniteWeilBernsteinProbe] {
        &self.probes
    }
}

#[derive(Debug)]
pub enum FiniteWeilBernsteinProbeError {
    EmptyIndexSet,
    IndicesNotStrictlyIncreasing { previous: usize, next: usize },
    IndexOutOfRange { degree: usize, index: usize },
    ParentDimensionOverflow { degree: usize },
    Parent(FiniteWeilGeneralizedSpectrumError),
    Boundary(WeilBoundaryError),
    NonPositiveGramNorm { index: usize, value: f64 },
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilBernsteinProbeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyIndexSet => write!(f, "Bernstein probe index set must be non-empty"),
            Self::IndicesNotStrictlyIncreasing { previous, next } => write!(
                f,
                "Bernstein probe indices must be strictly increasing: previous={previous}, next={next}"
            ),
            Self::IndexOutOfRange { degree, index } => write!(
                f,
                "Bernstein index {index} is outside the degree-{degree} family"
            ),
            Self::ParentDimensionOverflow { degree } => write!(
                f,
                "Bernstein degree {degree} cannot be converted to a parent dimension"
            ),
            Self::Parent(error) => write!(f, "parent generalized-spectrum audit failed: {error}"),
            Self::Boundary(error) => write!(f, "Bernstein boundary audit failed: {error}"),
            Self::NonPositiveGramNorm { index, value } => write!(
                f,
                "Bernstein probe {index} has non-positive numerical Gram norm squared: {value}"
            ),
            Self::NonFiniteEvaluation { stage, value } => {
                write!(f, "non-finite Bernstein probe value at {stage}: {value}")
            }
        }
    }
}

impl std::error::Error for FiniteWeilBernsteinProbeError {}

impl From<FiniteWeilGeneralizedSpectrumError> for FiniteWeilBernsteinProbeError {
    fn from(value: FiniteWeilGeneralizedSpectrumError) -> Self {
        Self::Parent(value)
    }
}

impl From<WeilBoundaryError> for FiniteWeilBernsteinProbeError {
    fn from(value: WeilBoundaryError) -> Self {
        Self::Boundary(value)
    }
}

/// Return the coefficients `c_j` in
/// `B_{D,k}(t) = sum_{j=0}^D c_j P_j(2t-1)`.
pub fn bernstein_legendre_coefficients(
    degree: usize,
    index: usize,
) -> Result<Vec<f64>, FiniteWeilBernsteinProbeError> {
    if index > degree {
        return Err(FiniteWeilBernsteinProbeError::IndexOutOfRange { degree, index });
    }

    let mut coefficients = Vec::with_capacity(degree + 1);
    for legendre_degree in 0..=degree {
        let mut integral = 0.0_f64;
        for monomial_degree in 0..=legendre_degree {
            let sign = if (legendre_degree + monomial_degree).is_multiple_of(2) {
                1.0
            } else {
                -1.0
            };
            let shifted_legendre_coefficient = sign
                * binomial_as_f64(legendre_degree, monomial_degree)
                * binomial_as_f64(legendre_degree + monomial_degree, monomial_degree);
            let moment = bernstein_monomial_moment(degree, index, monomial_degree);
            integral += shifted_legendre_coefficient * moment;
        }
        let coefficient = (2 * legendre_degree + 1) as f64 * integral;
        checked_finite("Bernstein Legendre coefficient", coefficient)?;
        coefficients.push(coefficient);
    }
    Ok(coefficients)
}

/// Evaluate localized Bernstein directions against one shared parent finite Weil
/// pairing/Gram computation of dimension `degree + 1`.
pub fn audit_finite_weil_bernstein_probes(
    bump: CompactArchimedeanBump,
    degree: usize,
    indices: &[usize],
    correlation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
    gram_order: usize,
) -> Result<FiniteWeilBernsteinProbeAudit, FiniteWeilBernsteinProbeError> {
    validate_indices(degree, indices)?;
    let parent_dimension = degree
        .checked_add(1)
        .ok_or(FiniteWeilBernsteinProbeError::ParentDimensionOverflow { degree })?;
    let parent = audit_finite_weil_generalized_spectrum(
        bump,
        parent_dimension,
        correlation_order,
        archimedean_order,
        boundary_order,
        gram_order,
    )?;

    let mut parent_moments = Vec::with_capacity(parent_dimension);
    for legendre_degree in 0..parent_dimension {
        parent_moments.push(
            CompactWeilBasisFunction::new(bump, legendre_degree)
                .boundary_moments(boundary_order)?,
        );
    }

    let mut probes = Vec::with_capacity(indices.len());
    for &index in indices {
        let coefficients = bernstein_legendre_coefficients(degree, index)?;
        let raw_quadratic_value =
            quadratic_form_from_parent(&coefficients, parent_dimension, |row, col| {
                parent
                    .pairing()
                    .entry(row, col)
                    .expect("Bernstein coefficient index is inside parent pairing matrix")
            });
        checked_finite("Bernstein raw quadratic value", raw_quadratic_value)?;

        let gram_norm_squared =
            quadratic_form_from_parent(&coefficients, parent_dimension, |row, col| {
                parent
                    .gram_entry(row, col)
                    .expect("Bernstein coefficient index is inside parent Gram matrix")
            });
        checked_finite("Bernstein Gram norm squared", gram_norm_squared)?;
        if gram_norm_squared <= 0.0 {
            return Err(FiniteWeilBernsteinProbeError::NonPositiveGramNorm {
                index,
                value: gram_norm_squared,
            });
        }

        let generalized_rayleigh_quotient = raw_quadratic_value / gram_norm_squared;
        checked_finite(
            "Bernstein generalized Rayleigh quotient",
            generalized_rayleigh_quotient,
        )?;

        let mut boundary_plus_half = 0.0_f64;
        let mut boundary_minus_half = 0.0_f64;
        for (&coefficient, moments) in coefficients.iter().zip(parent_moments.iter()) {
            boundary_plus_half += coefficient * moments.plus_half;
            boundary_minus_half += coefficient * moments.minus_half;
        }
        checked_finite("Bernstein +1/2 boundary moment", boundary_plus_half)?;
        checked_finite("Bernstein -1/2 boundary moment", boundary_minus_half)?;
        let max_boundary_residual = boundary_plus_half.abs().max(boundary_minus_half.abs());
        let max_reconstruction_residual =
            bernstein_reconstruction_residual(degree, index, &coefficients);
        checked_finite(
            "Bernstein reconstruction residual",
            max_reconstruction_residual,
        )?;
        let coefficient_l1_norm = coefficients.iter().map(|value| value.abs()).sum::<f64>();
        checked_finite("Bernstein coefficient L1 norm", coefficient_l1_norm)?;

        probes.push(FiniteWeilBernsteinProbe {
            index,
            legendre_coefficients: coefficients,
            coefficient_l1_norm,
            raw_quadratic_value,
            gram_norm_squared,
            generalized_rayleigh_quotient,
            boundary_plus_half,
            boundary_minus_half,
            max_boundary_residual,
            max_reconstruction_residual,
        });
    }

    Ok(FiniteWeilBernsteinProbeAudit {
        degree,
        parent_dimension,
        parent_gram_condition_number: parent.gram_condition_number(),
        parent_max_raw_pairing_asymmetry: parent.pairing().max_raw_pairing_asymmetry(),
        parent_max_whitened_asymmetry: parent.max_whitened_asymmetry(),
        probes,
    })
}

fn validate_indices(degree: usize, indices: &[usize]) -> Result<(), FiniteWeilBernsteinProbeError> {
    if indices.is_empty() {
        return Err(FiniteWeilBernsteinProbeError::EmptyIndexSet);
    }
    for &index in indices {
        if index > degree {
            return Err(FiniteWeilBernsteinProbeError::IndexOutOfRange { degree, index });
        }
    }
    for pair in indices.windows(2) {
        if pair[0] >= pair[1] {
            return Err(
                FiniteWeilBernsteinProbeError::IndicesNotStrictlyIncreasing {
                    previous: pair[0],
                    next: pair[1],
                },
            );
        }
    }
    Ok(())
}

fn quadratic_form_from_parent(
    coefficients: &[f64],
    dimension: usize,
    mut entry: impl FnMut(usize, usize) -> f64,
) -> f64 {
    let mut total = 0.0_f64;
    for row in 0..dimension {
        for col in 0..dimension {
            total += coefficients[row] * entry(row, col) * coefficients[col];
        }
    }
    total
}

fn bernstein_monomial_moment(degree: usize, index: usize, monomial_degree: usize) -> f64 {
    let mut moment = 1.0 / (degree + 1) as f64;
    for step in 1..=monomial_degree {
        moment *= (index + step) as f64 / (degree + step + 1) as f64;
    }
    moment
}

fn binomial_as_f64(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    let mut value = 1.0_f64;
    for step in 1..=k {
        value *= (n - k + step) as f64 / step as f64;
    }
    value
}

fn bernstein_reconstruction_residual(degree: usize, index: usize, coefficients: &[f64]) -> f64 {
    let mut maximum = 0.0_f64;
    for sample in 0..=RECONSTRUCTION_INTERVALS {
        let t = sample as f64 / RECONSTRUCTION_INTERVALS as f64;
        let direct = bernstein_value(degree, index, t);
        let expanded = coefficients
            .iter()
            .enumerate()
            .map(|(legendre_degree, &coefficient)| {
                coefficient * shifted_legendre_value(legendre_degree, t)
            })
            .sum::<f64>();
        maximum = maximum.max((direct - expanded).abs());
    }
    maximum
}

fn bernstein_value(degree: usize, index: usize, t: f64) -> f64 {
    binomial_as_f64(degree, index)
        * integer_power(t, index)
        * integer_power(1.0 - t, degree - index)
}

fn integer_power(base: f64, exponent: usize) -> f64 {
    let mut value = 1.0_f64;
    for _ in 0..exponent {
        value *= base;
    }
    value
}

fn shifted_legendre_value(degree: usize, t: f64) -> f64 {
    let x = 2.0 * t - 1.0;
    if degree == 0 {
        return 1.0;
    }
    if degree == 1 {
        return x;
    }

    let mut p_nm2 = 1.0_f64;
    let mut p_nm1 = x;
    for n in 2..=degree {
        let p = ((2 * n - 1) as f64 * x * p_nm1 - (n - 1) as f64 * p_nm2) / n as f64;
        p_nm2 = p_nm1;
        p_nm1 = p;
    }
    p_nm1
}

fn checked_finite(stage: &'static str, value: f64) -> Result<(), FiniteWeilBernsteinProbeError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilBernsteinProbeError::NonFiniteEvaluation { stage, value })
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
    fn degree_two_coefficients_match_closed_forms() {
        let left = bernstein_legendre_coefficients(2, 0).unwrap();
        let middle = bernstein_legendre_coefficients(2, 1).unwrap();
        let right = bernstein_legendre_coefficients(2, 2).unwrap();

        let tolerance = 2.0e-14;
        for (&actual, expected) in left.iter().zip([1.0 / 3.0, -0.5, 1.0 / 6.0]) {
            assert!((actual - expected).abs() <= tolerance);
        }
        for (&actual, expected) in middle.iter().zip([1.0 / 3.0, 0.0, -1.0 / 3.0]) {
            assert!((actual - expected).abs() <= tolerance);
        }
        for (&actual, expected) in right.iter().zip([1.0 / 3.0, 0.5, 1.0 / 6.0]) {
            assert!((actual - expected).abs() <= tolerance);
        }
    }

    #[test]
    fn localized_coefficients_reconstruct_bernstein_polynomials() {
        for index in [0_usize, 2, 4, 6] {
            let coefficients = bernstein_legendre_coefficients(6, index).unwrap();
            assert!(bernstein_reconstruction_residual(6, index, &coefficients) <= 2.0e-12);
        }
    }

    #[test]
    fn finite_bernstein_probes_are_audited_without_sign_assumption() {
        let audit =
            audit_finite_weil_bernstein_probes(bump(), 2, &[0, 1, 2], 32, 32, 48, 48).unwrap();
        assert_eq!(audit.degree(), 2);
        assert_eq!(audit.parent_dimension(), 3);
        assert_eq!(audit.probes().len(), 3);
        assert!(audit.parent_gram_condition_number().is_finite());
        for probe in audit.probes() {
            assert!(probe.raw_quadratic_value().is_finite());
            assert!(probe.gram_norm_squared() > 0.0);
            assert!(probe.generalized_rayleigh_quotient().is_finite());
            assert!(probe.max_boundary_residual().is_finite());
            assert!(probe.max_reconstruction_residual() <= 2.0e-12);
        }
    }
}
