//! Gram-normalized finite Riemann--Weil spectrum.
//!
//! The raw pairing matrix from `weil_quadratic_matrix` is a matrix of a
//! quadratic form in a non-orthonormal basis. Its inertia is meaningful on the
//! declared finite span, but ordinary eigenvalue magnitudes depend on basis
//! scaling. This module adds an explicit positive Gram matrix
//!
//! `G_ij = integral h_i(rho) h_j(rho) d^*rho`,
//!
//! with `d^*rho = d rho / rho`, and computes the equivalent self-adjoint
//! normalized matrix
//!
//! `B = G^(-1/2) A G^(-1/2)`.
//!
//! The eigenvalues of `B` are the generalized eigenvalues of `A v = lambda G v`.
//! The multiplicative `L^2(d^*rho)` norm used here is an explicit numerical
//! normalization choice. It is not identified with the source semilocal
//! Hilbert space `L^2(X_S)` and does not upgrade finite-basis evidence to Weil
//! positivity or RH.

use std::fmt;

use faer::{Mat, Side, linalg::solvers::SelfAdjointEigen};

use crate::quadrature::{GaussLegendreUnit, QuadratureError};
use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_boundary::WeilBoundaryError;
use crate::weil_quadratic_matrix::{
    CompactWeilBasisFunction, FiniteWeilMatrixError, FiniteWeilQuadraticMatrixAudit,
    audit_finite_weil_quadratic_matrix,
};

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilGeneralizedSpectrumAudit {
    pairing: FiniteWeilQuadraticMatrixAudit,
    gram_entries: Vec<f64>,
    gram_eigenvalues: Vec<f64>,
    gram_condition_number: f64,
    generalized_eigenvalues: Vec<f64>,
    max_whitened_asymmetry: f64,
}

impl FiniteWeilGeneralizedSpectrumAudit {
    #[inline]
    pub fn dimension(&self) -> usize {
        self.pairing.dimension()
    }

    #[inline]
    pub fn pairing(&self) -> &FiniteWeilQuadraticMatrixAudit {
        &self.pairing
    }

    pub fn gram_entry(&self, row: usize, col: usize) -> Option<f64> {
        let dimension = self.dimension();
        (row < dimension && col < dimension).then_some(self.gram_entries[row * dimension + col])
    }

    #[inline]
    pub fn gram_eigenvalues(&self) -> &[f64] {
        &self.gram_eigenvalues
    }

    #[inline]
    pub fn gram_condition_number(&self) -> f64 {
        self.gram_condition_number
    }

    #[inline]
    pub fn generalized_eigenvalues(&self) -> &[f64] {
        &self.generalized_eigenvalues
    }

    #[inline]
    pub fn minimum_generalized_eigenvalue(&self) -> f64 {
        self.generalized_eigenvalues[0]
    }

    #[inline]
    pub fn max_whitened_asymmetry(&self) -> f64 {
        self.max_whitened_asymmetry
    }
}

#[derive(Debug)]
pub enum FiniteWeilGeneralizedSpectrumError {
    Pairing(FiniteWeilMatrixError),
    Quadrature(QuadratureError),
    Boundary(WeilBoundaryError),
    GramDecompositionFailed,
    GramNotPositiveDefinite { minimum_eigenvalue: f64 },
    NormalizedDecompositionFailed,
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilGeneralizedSpectrumError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Pairing(error) => write!(f, "finite Weil pairing matrix failed: {error}"),
            Self::Quadrature(error) => write!(f, "Gram quadrature construction failed: {error:?}"),
            Self::Boundary(error) => write!(f, "Gram basis evaluation failed: {error}"),
            Self::GramDecompositionFailed => write!(f, "Gram eigendecomposition failed"),
            Self::GramNotPositiveDefinite { minimum_eigenvalue } => write!(
                f,
                "Gram matrix is not numerically positive definite: lambda_min={minimum_eigenvalue}"
            ),
            Self::NormalizedDecompositionFailed => {
                write!(f, "normalized Weil eigendecomposition failed")
            }
            Self::NonFiniteEvaluation { stage, value } => {
                write!(
                    f,
                    "non-finite generalized-spectrum value at {stage}: {value}"
                )
            }
        }
    }
}

impl std::error::Error for FiniteWeilGeneralizedSpectrumError {}

impl From<FiniteWeilMatrixError> for FiniteWeilGeneralizedSpectrumError {
    fn from(value: FiniteWeilMatrixError) -> Self {
        Self::Pairing(value)
    }
}

impl From<QuadratureError> for FiniteWeilGeneralizedSpectrumError {
    fn from(value: QuadratureError) -> Self {
        Self::Quadrature(value)
    }
}

impl From<WeilBoundaryError> for FiniteWeilGeneralizedSpectrumError {
    fn from(value: WeilBoundaryError) -> Self {
        Self::Boundary(value)
    }
}

/// Build the raw finite Weil pairing matrix and normalize it with the
/// multiplicative `L^2(d rho / rho)` Gram matrix of the same basis.
pub fn audit_finite_weil_generalized_spectrum(
    bump: CompactArchimedeanBump,
    dimension: usize,
    correlation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
    gram_order: usize,
) -> Result<FiniteWeilGeneralizedSpectrumAudit, FiniteWeilGeneralizedSpectrumError> {
    let pairing = audit_finite_weil_quadratic_matrix(
        bump,
        dimension,
        correlation_order,
        archimedean_order,
        boundary_order,
    )?;

    let gram_entries = build_multiplicative_gram(bump, dimension, gram_order)?;
    let raw_entries = dense_entries_from_pairing(&pairing);
    let normalized = generalized_spectrum_from_dense_pair(&raw_entries, &gram_entries, dimension)?;

    Ok(FiniteWeilGeneralizedSpectrumAudit {
        pairing,
        gram_entries,
        gram_eigenvalues: normalized.gram_eigenvalues,
        gram_condition_number: normalized.gram_condition_number,
        generalized_eigenvalues: normalized.generalized_eigenvalues,
        max_whitened_asymmetry: normalized.max_whitened_asymmetry,
    })
}

fn build_multiplicative_gram(
    bump: CompactArchimedeanBump,
    dimension: usize,
    quadrature_order: usize,
) -> Result<Vec<f64>, FiniteWeilGeneralizedSpectrumError> {
    let quadrature = GaussLegendreUnit::new(quadrature_order)?;
    let support = bump.support();
    let log_lower = support.log_lower();
    let log_span = support.log_upper() - log_lower;
    checked_finite("Gram log-span", log_span)?;

    let basis = (0..dimension)
        .map(|degree| CompactWeilBasisFunction::new(bump, degree))
        .collect::<Vec<_>>();
    let mut sampled = vec![0.0_f64; dimension * quadrature.order()];

    for (node_index, &node) in quadrature.nodes().iter().enumerate() {
        let rho = (log_lower + log_span * node).exp();
        for (degree, &function) in basis.iter().enumerate() {
            let value = function.value(rho)?;
            checked_finite("Gram basis sample", value)?;
            sampled[degree * quadrature.order() + node_index] = value;
        }
    }

    let mut gram = vec![0.0_f64; dimension * dimension];
    for i in 0..dimension {
        for j in i..dimension {
            let mut sum = 0.0_f64;
            for (node_index, &weight) in quadrature.weights().iter().enumerate() {
                sum += weight
                    * sampled[i * quadrature.order() + node_index]
                    * sampled[j * quadrature.order() + node_index];
            }
            let value = log_span * sum;
            checked_finite("Gram entry", value)?;
            gram[i * dimension + j] = value;
            gram[j * dimension + i] = value;
        }
    }
    Ok(gram)
}

fn dense_entries_from_pairing(pairing: &FiniteWeilQuadraticMatrixAudit) -> Vec<f64> {
    let dimension = pairing.dimension();
    let mut entries = vec![0.0_f64; dimension * dimension];
    for i in 0..dimension {
        for j in 0..dimension {
            entries[i * dimension + j] = pairing
                .entry(i, j)
                .expect("indices are inside the finite pairing matrix");
        }
    }
    entries
}

struct NormalizedSpectrum {
    gram_eigenvalues: Vec<f64>,
    gram_condition_number: f64,
    generalized_eigenvalues: Vec<f64>,
    max_whitened_asymmetry: f64,
}

fn generalized_spectrum_from_dense_pair(
    a: &[f64],
    g: &[f64],
    dimension: usize,
) -> Result<NormalizedSpectrum, FiniteWeilGeneralizedSpectrumError> {
    let gram_matrix = Mat::from_fn(dimension, dimension, |i, j| g[i * dimension + j]);
    let gram_decomposition = SelfAdjointEigen::new(gram_matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilGeneralizedSpectrumError::GramDecompositionFailed)?;
    let gram_diagonal = gram_decomposition.S().column_vector();
    let gram_vectors = gram_decomposition.U();

    let raw_gram_eigenvalues = (0..dimension)
        .map(|index| gram_diagonal[index])
        .collect::<Vec<_>>();
    for &value in &raw_gram_eigenvalues {
        checked_finite("Gram eigenvalue", value)?;
    }
    let minimum_gram = raw_gram_eigenvalues
        .iter()
        .copied()
        .min_by(f64::total_cmp)
        .ok_or(FiniteWeilGeneralizedSpectrumError::GramDecompositionFailed)?;
    if minimum_gram <= 0.0 {
        return Err(
            FiniteWeilGeneralizedSpectrumError::GramNotPositiveDefinite {
                minimum_eigenvalue: minimum_gram,
            },
        );
    }
    let maximum_gram = raw_gram_eigenvalues
        .iter()
        .copied()
        .max_by(f64::total_cmp)
        .ok_or(FiniteWeilGeneralizedSpectrumError::GramDecompositionFailed)?;
    let gram_condition_number = maximum_gram / minimum_gram;
    checked_finite("Gram condition number", gram_condition_number)?;

    let mut inverse_sqrt = vec![0.0_f64; dimension * dimension];
    for i in 0..dimension {
        for j in 0..dimension {
            let mut sum = 0.0_f64;
            for k in 0..dimension {
                sum += gram_vectors[(i, k)] * gram_vectors[(j, k)] / raw_gram_eigenvalues[k].sqrt();
            }
            checked_finite("Gram inverse square root", sum)?;
            inverse_sqrt[i * dimension + j] = sum;
        }
    }

    let left_product = multiply_dense(&inverse_sqrt, a, dimension);
    let whitened_raw = multiply_dense(&left_product, &inverse_sqrt, dimension);
    let mut max_whitened_asymmetry = 0.0_f64;
    let mut whitened = whitened_raw.clone();
    for i in 0..dimension {
        for j in i..dimension {
            let forward = whitened_raw[i * dimension + j];
            let reverse = whitened_raw[j * dimension + i];
            max_whitened_asymmetry = max_whitened_asymmetry.max((forward - reverse).abs());
            let symmetric = 0.5 * (forward + reverse);
            whitened[i * dimension + j] = symmetric;
            whitened[j * dimension + i] = symmetric;
        }
    }
    checked_finite("whitened asymmetry", max_whitened_asymmetry)?;

    let whitened_matrix = Mat::from_fn(dimension, dimension, |i, j| whitened[i * dimension + j]);
    let decomposition = SelfAdjointEigen::new(whitened_matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilGeneralizedSpectrumError::NormalizedDecompositionFailed)?;
    let diagonal = decomposition.S().column_vector();
    let mut generalized_eigenvalues = (0..dimension)
        .map(|index| diagonal[index])
        .collect::<Vec<_>>();
    for &value in &generalized_eigenvalues {
        checked_finite("generalized eigenvalue", value)?;
    }
    generalized_eigenvalues.sort_by(f64::total_cmp);

    let mut gram_eigenvalues = raw_gram_eigenvalues;
    gram_eigenvalues.sort_by(f64::total_cmp);

    Ok(NormalizedSpectrum {
        gram_eigenvalues,
        gram_condition_number,
        generalized_eigenvalues,
        max_whitened_asymmetry,
    })
}

fn multiply_dense(left: &[f64], right: &[f64], dimension: usize) -> Vec<f64> {
    let mut output = vec![0.0_f64; dimension * dimension];
    for i in 0..dimension {
        for k in 0..dimension {
            let left_value = left[i * dimension + k];
            for j in 0..dimension {
                output[i * dimension + j] += left_value * right[k * dimension + j];
            }
        }
    }
    output
}

fn checked_finite(
    stage: &'static str,
    value: f64,
) -> Result<(), FiniteWeilGeneralizedSpectrumError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilGeneralizedSpectrumError::NonFiniteEvaluation { stage, value })
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
    fn one_dimensional_generalized_value_is_pairing_over_gram_norm() {
        let audit = audit_finite_weil_generalized_spectrum(bump(), 1, 64, 64, 96, 96).unwrap();
        let expected = audit.pairing().entry(0, 0).unwrap() / audit.gram_entry(0, 0).unwrap();
        assert!((audit.minimum_generalized_eigenvalue() - expected).abs() <= 5.0e-12);
        assert!(audit.gram_eigenvalues()[0] > 0.0);
    }

    #[test]
    fn generalized_spectrum_is_invariant_under_diagonal_basis_rescaling() {
        let dimension = 3;
        let a = vec![4.0, 1.0, -0.5, 1.0, 3.0, 0.25, -0.5, 0.25, 2.0];
        let g = vec![2.0, 0.2, 0.1, 0.2, 1.5, -0.05, 0.1, -0.05, 1.0];
        let original = generalized_spectrum_from_dense_pair(&a, &g, dimension).unwrap();

        let scales = [0.25_f64, 3.0, 7.0];
        let mut scaled_a = vec![0.0_f64; dimension * dimension];
        let mut scaled_g = vec![0.0_f64; dimension * dimension];
        for i in 0..dimension {
            for j in 0..dimension {
                let scale = scales[i] * scales[j];
                scaled_a[i * dimension + j] = scale * a[i * dimension + j];
                scaled_g[i * dimension + j] = scale * g[i * dimension + j];
            }
        }
        let scaled = generalized_spectrum_from_dense_pair(&scaled_a, &scaled_g, dimension).unwrap();

        for (&left, &right) in original
            .generalized_eigenvalues
            .iter()
            .zip(scaled.generalized_eigenvalues.iter())
        {
            assert!((left - right).abs() <= 2.0e-12 * left.abs().max(right.abs()).max(1.0));
        }
    }

    #[test]
    fn compact_four_dimensional_generalized_spectrum_matches_independent_regression() {
        let audit = audit_finite_weil_generalized_spectrum(bump(), 4, 96, 96, 128, 128).unwrap();
        assert_eq!(audit.dimension(), 4);
        assert!(audit.pairing().max_boundary_residual() <= 2.0e-10);
        assert!(audit.pairing().max_raw_pairing_asymmetry() <= 5.0e-12);
        assert!(audit.max_whitened_asymmetry() <= 5.0e-12);
        assert!(audit.gram_eigenvalues().iter().all(|value| *value > 0.0));
        assert!(audit.gram_condition_number().is_finite());

        let expected_minimum = 3.707_307_755_581_390_4e-3;
        assert!(
            (audit.minimum_generalized_eigenvalue() - expected_minimum).abs() <= 5.0e-7,
            "generalized lambda_min={:.15e}",
            audit.minimum_generalized_eigenvalue()
        );
    }
}
