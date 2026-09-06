//! Generic coefficient-defined finite Riemann--Weil subspace audit.
//!
//! Given an already computed parent Legendre pairing/Gram pair `(A,G)` and a
//! matrix `C` whose columns are coefficient vectors in that parent basis, this
//! module evaluates
//!
//! `A_C = C^T A C`, `G_C = C^T G C`
//!
//! and solves the corresponding Gram-normalized generalized spectrum.  This is
//! finite linear algebra on the declared parent space.  It does not identify a
//! complete-space operator or strengthen finite positivity into Weil positivity
//! or RH.

use std::fmt;

use faer::{Mat, Side, linalg::solvers::SelfAdjointEigen};

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_boundary::WeilBoundaryError;
use crate::weil_generalized_spectrum::{
    FiniteWeilGeneralizedSpectrumAudit, FiniteWeilGeneralizedSpectrumError,
};
use crate::weil_quadratic_matrix::CompactWeilBasisFunction;

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilCoefficientSubspaceAudit {
    parent_dimension: usize,
    dimension: usize,
    raw_eigenvalues: Vec<f64>,
    gram_eigenvalues: Vec<f64>,
    generalized_eigenvalues: Vec<f64>,
    gram_condition_number: f64,
    max_boundary_residual: f64,
    max_whitened_asymmetry: f64,
    parent_max_raw_pairing_asymmetry: f64,
    leading_legendre_generalized_minimum: f64,
    leading_legendre_gram_condition_number: f64,
}

impl FiniteWeilCoefficientSubspaceAudit {
    #[inline]
    pub const fn parent_dimension(&self) -> usize {
        self.parent_dimension
    }

    #[inline]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    #[inline]
    pub fn raw_eigenvalues(&self) -> &[f64] {
        &self.raw_eigenvalues
    }

    #[inline]
    pub fn gram_eigenvalues(&self) -> &[f64] {
        &self.gram_eigenvalues
    }

    #[inline]
    pub fn generalized_eigenvalues(&self) -> &[f64] {
        &self.generalized_eigenvalues
    }

    #[inline]
    pub fn minimum_raw_eigenvalue(&self) -> f64 {
        self.raw_eigenvalues[0]
    }

    #[inline]
    pub fn minimum_generalized_eigenvalue(&self) -> f64 {
        self.generalized_eigenvalues[0]
    }

    #[inline]
    pub const fn gram_condition_number(&self) -> f64 {
        self.gram_condition_number
    }

    #[inline]
    pub const fn max_boundary_residual(&self) -> f64 {
        self.max_boundary_residual
    }

    #[inline]
    pub const fn max_whitened_asymmetry(&self) -> f64 {
        self.max_whitened_asymmetry
    }

    #[inline]
    pub const fn parent_max_raw_pairing_asymmetry(&self) -> f64 {
        self.parent_max_raw_pairing_asymmetry
    }

    #[inline]
    pub const fn leading_legendre_generalized_minimum(&self) -> f64 {
        self.leading_legendre_generalized_minimum
    }

    #[inline]
    pub const fn leading_legendre_gram_condition_number(&self) -> f64 {
        self.leading_legendre_gram_condition_number
    }
}

#[derive(Debug)]
pub enum FiniteWeilCoefficientSubspaceError {
    EmptyCoefficientSet,
    EmptyCoefficientVector { column: usize },
    CoefficientDimensionMismatch {
        column: usize,
        expected: usize,
        actual: usize,
    },
    Boundary(WeilBoundaryError),
    Generalized(FiniteWeilGeneralizedSpectrumError),
    RawDecompositionFailed,
    GramDecompositionFailed,
    GramNotPositiveDefinite { minimum_eigenvalue: f64 },
    NormalizedDecompositionFailed,
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilCoefficientSubspaceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyCoefficientSet => {
                write!(f, "finite Weil coefficient subspace must contain at least one vector")
            }
            Self::EmptyCoefficientVector { column } => {
                write!(f, "finite Weil coefficient vector {column} is empty")
            }
            Self::CoefficientDimensionMismatch {
                column,
                expected,
                actual,
            } => write!(
                f,
                "finite Weil coefficient vector {column} has length {actual}; expected {expected}"
            ),
            Self::Boundary(error) => write!(f, "coefficient-subspace boundary audit failed: {error}"),
            Self::Generalized(error) => {
                write!(f, "coefficient-subspace leading control failed: {error}")
            }
            Self::RawDecompositionFailed => {
                write!(f, "coefficient-subspace raw eigendecomposition failed")
            }
            Self::GramDecompositionFailed => {
                write!(f, "coefficient-subspace Gram eigendecomposition failed")
            }
            Self::GramNotPositiveDefinite { minimum_eigenvalue } => write!(
                f,
                "coefficient-subspace Gram matrix is not numerically positive definite: lambda_min={minimum_eigenvalue}"
            ),
            Self::NormalizedDecompositionFailed => {
                write!(f, "coefficient-subspace normalized eigendecomposition failed")
            }
            Self::NonFiniteEvaluation { stage, value } => {
                write!(f, "non-finite coefficient-subspace value at {stage}: {value}")
            }
        }
    }
}

impl std::error::Error for FiniteWeilCoefficientSubspaceError {}

impl From<WeilBoundaryError> for FiniteWeilCoefficientSubspaceError {
    fn from(value: WeilBoundaryError) -> Self {
        Self::Boundary(value)
    }
}

impl From<FiniteWeilGeneralizedSpectrumError> for FiniteWeilCoefficientSubspaceError {
    fn from(value: FiniteWeilGeneralizedSpectrumError) -> Self {
        Self::Generalized(value)
    }
}

/// Restrict one already-computed parent finite Weil problem to coefficient
/// vectors expressed in its Legendre basis.
pub fn audit_finite_weil_coefficient_subspace(
    bump: CompactArchimedeanBump,
    parent: &FiniteWeilGeneralizedSpectrumAudit,
    coefficients: &[Vec<f64>],
    boundary_order: usize,
) -> Result<FiniteWeilCoefficientSubspaceAudit, FiniteWeilCoefficientSubspaceError> {
    if coefficients.is_empty() {
        return Err(FiniteWeilCoefficientSubspaceError::EmptyCoefficientSet);
    }
    let parent_dimension = parent.dimension();
    for (column, vector) in coefficients.iter().enumerate() {
        if vector.is_empty() {
            return Err(FiniteWeilCoefficientSubspaceError::EmptyCoefficientVector { column });
        }
        if vector.len() != parent_dimension {
            return Err(
                FiniteWeilCoefficientSubspaceError::CoefficientDimensionMismatch {
                    column,
                    expected: parent_dimension,
                    actual: vector.len(),
                },
            );
        }
        for &value in vector {
            checked_finite("coefficient", value)?;
        }
    }

    let dimension = coefficients.len();
    let mut transformed_a = vec![0.0_f64; dimension * dimension];
    let mut transformed_g = vec![0.0_f64; dimension * dimension];
    for row in 0..dimension {
        for col in row..dimension {
            let a_value = bilinear_parent_form(&coefficients[row], &coefficients[col], |i, j| {
                parent
                    .pairing()
                    .entry(i, j)
                    .expect("coefficient index is inside parent pairing matrix")
            });
            let g_value = bilinear_parent_form(&coefficients[row], &coefficients[col], |i, j| {
                parent
                    .gram_entry(i, j)
                    .expect("coefficient index is inside parent Gram matrix")
            });
            checked_finite("transformed coefficient pairing entry", a_value)?;
            checked_finite("transformed coefficient Gram entry", g_value)?;
            transformed_a[row * dimension + col] = a_value;
            transformed_a[col * dimension + row] = a_value;
            transformed_g[row * dimension + col] = g_value;
            transformed_g[col * dimension + row] = g_value;
        }
    }

    let mut parent_moments = Vec::with_capacity(parent_dimension);
    for degree in 0..parent_dimension {
        parent_moments
            .push(CompactWeilBasisFunction::new(bump, degree).boundary_moments(boundary_order)?);
    }
    let mut max_boundary_residual = 0.0_f64;
    for coefficient_vector in coefficients {
        let plus_half = coefficient_vector
            .iter()
            .zip(parent_moments.iter())
            .map(|(&coefficient, moments)| coefficient * moments.plus_half)
            .sum::<f64>();
        let minus_half = coefficient_vector
            .iter()
            .zip(parent_moments.iter())
            .map(|(&coefficient, moments)| coefficient * moments.minus_half)
            .sum::<f64>();
        checked_finite("coefficient-subspace +1/2 boundary moment", plus_half)?;
        checked_finite("coefficient-subspace -1/2 boundary moment", minus_half)?;
        max_boundary_residual = max_boundary_residual.max(plus_half.abs().max(minus_half.abs()));
    }

    let solved = solve_transformed_pair(&transformed_a, &transformed_g, dimension)?;
    let leading_control = parent.principal_spectrum(dimension)?;

    Ok(FiniteWeilCoefficientSubspaceAudit {
        parent_dimension,
        dimension,
        raw_eigenvalues: solved.raw_eigenvalues,
        gram_eigenvalues: solved.gram_eigenvalues,
        generalized_eigenvalues: solved.generalized_eigenvalues,
        gram_condition_number: solved.gram_condition_number,
        max_boundary_residual,
        max_whitened_asymmetry: solved.max_whitened_asymmetry,
        parent_max_raw_pairing_asymmetry: parent.pairing().max_raw_pairing_asymmetry(),
        leading_legendre_generalized_minimum: leading_control.generalized_minimum_eigenvalue(),
        leading_legendre_gram_condition_number: leading_control.gram_condition_number(),
    })
}

fn bilinear_parent_form(
    left: &[f64],
    right: &[f64],
    mut entry: impl FnMut(usize, usize) -> f64,
) -> f64 {
    let mut total = 0.0_f64;
    for (i, &left_value) in left.iter().enumerate() {
        for (j, &right_value) in right.iter().enumerate() {
            total += left_value * entry(i, j) * right_value;
        }
    }
    total
}

struct TransformedSpectrum {
    raw_eigenvalues: Vec<f64>,
    gram_eigenvalues: Vec<f64>,
    generalized_eigenvalues: Vec<f64>,
    gram_condition_number: f64,
    max_whitened_asymmetry: f64,
}

fn solve_transformed_pair(
    a: &[f64],
    g: &[f64],
    dimension: usize,
) -> Result<TransformedSpectrum, FiniteWeilCoefficientSubspaceError> {
    let raw_matrix = Mat::from_fn(dimension, dimension, |i, j| a[i * dimension + j]);
    let raw_decomposition = SelfAdjointEigen::new(raw_matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilCoefficientSubspaceError::RawDecompositionFailed)?;
    let raw_diagonal = raw_decomposition.S().column_vector();
    let mut raw_eigenvalues = (0..dimension)
        .map(|index| raw_diagonal[index])
        .collect::<Vec<_>>();
    for &value in &raw_eigenvalues {
        checked_finite("coefficient-subspace raw eigenvalue", value)?;
    }
    raw_eigenvalues.sort_by(f64::total_cmp);

    let gram_matrix = Mat::from_fn(dimension, dimension, |i, j| g[i * dimension + j]);
    let gram_decomposition = SelfAdjointEigen::new(gram_matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilCoefficientSubspaceError::GramDecompositionFailed)?;
    let gram_diagonal = gram_decomposition.S().column_vector();
    let gram_vectors = gram_decomposition.U();
    let raw_gram_eigenvalues = (0..dimension)
        .map(|index| gram_diagonal[index])
        .collect::<Vec<_>>();
    for &value in &raw_gram_eigenvalues {
        checked_finite("coefficient-subspace Gram eigenvalue", value)?;
    }
    let minimum_gram = raw_gram_eigenvalues
        .iter()
        .copied()
        .min_by(f64::total_cmp)
        .ok_or(FiniteWeilCoefficientSubspaceError::GramDecompositionFailed)?;
    if minimum_gram <= 0.0 {
        return Err(FiniteWeilCoefficientSubspaceError::GramNotPositiveDefinite {
            minimum_eigenvalue: minimum_gram,
        });
    }
    let maximum_gram = raw_gram_eigenvalues
        .iter()
        .copied()
        .max_by(f64::total_cmp)
        .ok_or(FiniteWeilCoefficientSubspaceError::GramDecompositionFailed)?;
    let gram_condition_number = maximum_gram / minimum_gram;
    checked_finite("coefficient-subspace Gram condition number", gram_condition_number)?;

    let mut inverse_sqrt = vec![0.0_f64; dimension * dimension];
    for i in 0..dimension {
        for j in 0..dimension {
            let mut sum = 0.0_f64;
            for k in 0..dimension {
                sum += gram_vectors[(i, k)] * gram_vectors[(j, k)] / raw_gram_eigenvalues[k].sqrt();
            }
            checked_finite("coefficient-subspace Gram inverse square root", sum)?;
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
    checked_finite("coefficient-subspace whitened asymmetry", max_whitened_asymmetry)?;

    let whitened_matrix = Mat::from_fn(dimension, dimension, |i, j| whitened[i * dimension + j]);
    let normalized_decomposition = SelfAdjointEigen::new(whitened_matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilCoefficientSubspaceError::NormalizedDecompositionFailed)?;
    let normalized_diagonal = normalized_decomposition.S().column_vector();
    let mut generalized_eigenvalues = (0..dimension)
        .map(|index| normalized_diagonal[index])
        .collect::<Vec<_>>();
    for &value in &generalized_eigenvalues {
        checked_finite("coefficient-subspace generalized eigenvalue", value)?;
    }
    generalized_eigenvalues.sort_by(f64::total_cmp);

    let mut gram_eigenvalues = raw_gram_eigenvalues;
    gram_eigenvalues.sort_by(f64::total_cmp);

    Ok(TransformedSpectrum {
        raw_eigenvalues,
        gram_eigenvalues,
        generalized_eigenvalues,
        gram_condition_number,
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
) -> Result<(), FiniteWeilCoefficientSubspaceError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilCoefficientSubspaceError::NonFiniteEvaluation { stage, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::PositiveRational;
    use crate::weil_generalized_spectrum::audit_finite_weil_generalized_spectrum;

    fn bump() -> CompactArchimedeanBump {
        CompactArchimedeanBump::new(
            PositiveRational::new(1, 2).unwrap(),
            PositiveRational::new(7, 2).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn identity_coefficients_reproduce_parent_generalized_spectrum() {
        let bump = bump();
        let parent = audit_finite_weil_generalized_spectrum(bump, 3, 28, 28, 40, 40).unwrap();
        let coefficients = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let audit =
            audit_finite_weil_coefficient_subspace(bump, &parent, &coefficients, 40).unwrap();

        assert_eq!(audit.dimension(), 3);
        assert_eq!(audit.parent_dimension(), 3);
        for (&left, &right) in audit
            .generalized_eigenvalues()
            .iter()
            .zip(parent.generalized_eigenvalues().iter())
        {
            assert!((left - right).abs() <= 2.0e-10 * left.abs().max(right.abs()).max(1.0));
        }
    }

    #[test]
    fn coefficient_dimension_mismatch_is_rejected() {
        let bump = bump();
        let parent = audit_finite_weil_generalized_spectrum(bump, 2, 20, 20, 28, 28).unwrap();
        let error = audit_finite_weil_coefficient_subspace(
            bump,
            &parent,
            &[vec![1.0, 0.0], vec![0.0]],
            28,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            FiniteWeilCoefficientSubspaceError::CoefficientDimensionMismatch { .. }
        ));
    }
}
