//! Generalized finite Riemann--Weil spectra on selected Bernstein subspaces.
//!
//! This layer strengthens the one-direction Bernstein probes by selecting several
//! localized Bernstein enrichments simultaneously.  If `C` contains their
//! shifted-Legendre coefficient vectors as columns, the restricted matrices are
//!
//! `A_B = C^T A C`, `G_B = C^T G C`.
//!
//! The parent `A` and `G` are computed once by the existing validated finite Weil
//! machinery.  A selected set with fewer than `D+1` degree-`D` Bernstein
//! functions generally spans a different finite subspace from the leading
//! Legendre space of the same dimension.  This is still a finite numerical
//! comparison and does not establish complete-space Weil positivity, density,
//! uniform convergence, Conjecture 4.1, or RH.

use std::fmt;

use faer::{Mat, Side, linalg::solvers::SelfAdjointEigen};

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_bernstein_probe::{FiniteWeilBernsteinProbeError, bernstein_legendre_coefficients};
use crate::weil_boundary::WeilBoundaryError;
use crate::weil_generalized_spectrum::{
    FiniteWeilGeneralizedSpectrumError, audit_finite_weil_generalized_spectrum,
};
use crate::weil_quadratic_matrix::CompactWeilBasisFunction;

/// A strictly increasing selection of Bernstein indices at one fixed degree.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BernsteinIndexSubspace {
    degree: usize,
    indices: Vec<usize>,
}

impl BernsteinIndexSubspace {
    pub fn new(
        degree: usize,
        indices: Vec<usize>,
    ) -> Result<Self, FiniteWeilBernsteinSubspaceError> {
        if indices.is_empty() {
            return Err(FiniteWeilBernsteinSubspaceError::EmptyIndexSet);
        }
        for &index in &indices {
            if index > degree {
                return Err(FiniteWeilBernsteinSubspaceError::IndexOutOfRange { degree, index });
            }
        }
        for pair in indices.windows(2) {
            if pair[0] >= pair[1] {
                return Err(
                    FiniteWeilBernsteinSubspaceError::IndicesNotStrictlyIncreasing {
                        previous: pair[0],
                        next: pair[1],
                    },
                );
            }
        }
        Ok(Self { degree, indices })
    }

    #[inline]
    pub fn degree(&self) -> usize {
        self.degree
    }

    #[inline]
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }

    #[inline]
    pub fn dimension(&self) -> usize {
        self.indices.len()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilBernsteinSubspaceAudit {
    degree: usize,
    indices: Vec<usize>,
    parent_dimension: usize,
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

impl FiniteWeilBernsteinSubspaceAudit {
    #[inline]
    pub fn degree(&self) -> usize {
        self.degree
    }

    #[inline]
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }

    #[inline]
    pub fn dimension(&self) -> usize {
        self.indices.len()
    }

    #[inline]
    pub fn parent_dimension(&self) -> usize {
        self.parent_dimension
    }

    #[inline]
    pub fn raw_eigenvalues(&self) -> &[f64] {
        &self.raw_eigenvalues
    }

    #[inline]
    pub fn minimum_raw_eigenvalue(&self) -> f64 {
        self.raw_eigenvalues[0]
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
    pub fn minimum_generalized_eigenvalue(&self) -> f64 {
        self.generalized_eigenvalues[0]
    }

    #[inline]
    pub fn gram_condition_number(&self) -> f64 {
        self.gram_condition_number
    }

    #[inline]
    pub fn max_boundary_residual(&self) -> f64 {
        self.max_boundary_residual
    }

    #[inline]
    pub fn max_whitened_asymmetry(&self) -> f64 {
        self.max_whitened_asymmetry
    }

    #[inline]
    pub fn parent_max_raw_pairing_asymmetry(&self) -> f64 {
        self.parent_max_raw_pairing_asymmetry
    }

    /// Dimension-matched leading Legendre control extracted from the same
    /// already-computed parent `(A,G)`.
    #[inline]
    pub fn leading_legendre_generalized_minimum(&self) -> f64 {
        self.leading_legendre_generalized_minimum
    }

    #[inline]
    pub fn leading_legendre_gram_condition_number(&self) -> f64 {
        self.leading_legendre_gram_condition_number
    }
}

#[derive(Debug)]
pub enum FiniteWeilBernsteinSubspaceError {
    EmptyIndexSet,
    IndicesNotStrictlyIncreasing { previous: usize, next: usize },
    IndexOutOfRange { degree: usize, index: usize },
    ParentDimensionOverflow { degree: usize },
    Parent(FiniteWeilGeneralizedSpectrumError),
    Probe(FiniteWeilBernsteinProbeError),
    Boundary(WeilBoundaryError),
    RawDecompositionFailed,
    GramDecompositionFailed,
    GramNotPositiveDefinite { minimum_eigenvalue: f64 },
    NormalizedDecompositionFailed,
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilBernsteinSubspaceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyIndexSet => write!(f, "Bernstein subspace index set must be non-empty"),
            Self::IndicesNotStrictlyIncreasing { previous, next } => write!(
                f,
                "Bernstein subspace indices must be strictly increasing: previous={previous}, next={next}"
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
            Self::Probe(error) => write!(f, "Bernstein coefficient construction failed: {error}"),
            Self::Boundary(error) => write!(f, "Bernstein subspace boundary audit failed: {error}"),
            Self::RawDecompositionFailed => {
                write!(f, "Bernstein subspace raw eigendecomposition failed")
            }
            Self::GramDecompositionFailed => {
                write!(f, "Bernstein subspace Gram eigendecomposition failed")
            }
            Self::GramNotPositiveDefinite { minimum_eigenvalue } => write!(
                f,
                "Bernstein subspace Gram matrix is not numerically positive definite: lambda_min={minimum_eigenvalue}"
            ),
            Self::NormalizedDecompositionFailed => {
                write!(f, "Bernstein subspace normalized eigendecomposition failed")
            }
            Self::NonFiniteEvaluation { stage, value } => {
                write!(f, "non-finite Bernstein subspace value at {stage}: {value}")
            }
        }
    }
}

impl std::error::Error for FiniteWeilBernsteinSubspaceError {}

impl From<FiniteWeilGeneralizedSpectrumError> for FiniteWeilBernsteinSubspaceError {
    fn from(value: FiniteWeilGeneralizedSpectrumError) -> Self {
        Self::Parent(value)
    }
}

impl From<FiniteWeilBernsteinProbeError> for FiniteWeilBernsteinSubspaceError {
    fn from(value: FiniteWeilBernsteinProbeError) -> Self {
        Self::Probe(value)
    }
}

impl From<WeilBoundaryError> for FiniteWeilBernsteinSubspaceError {
    fn from(value: WeilBoundaryError) -> Self {
        Self::Boundary(value)
    }
}

/// Build and solve the generalized finite Weil problem on a selected localized
/// Bernstein subspace, with a leading-Legendre control of the same dimension.
pub fn audit_finite_weil_bernstein_subspace(
    bump: CompactArchimedeanBump,
    subspace: &BernsteinIndexSubspace,
    correlation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
    gram_order: usize,
) -> Result<FiniteWeilBernsteinSubspaceAudit, FiniteWeilBernsteinSubspaceError> {
    let parent_dimension = subspace.degree().checked_add(1).ok_or(
        FiniteWeilBernsteinSubspaceError::ParentDimensionOverflow {
            degree: subspace.degree(),
        },
    )?;
    let parent = audit_finite_weil_generalized_spectrum(
        bump,
        parent_dimension,
        correlation_order,
        archimedean_order,
        boundary_order,
        gram_order,
    )?;

    let coefficients = subspace
        .indices()
        .iter()
        .map(|&index| bernstein_legendre_coefficients(subspace.degree(), index))
        .collect::<Result<Vec<_>, _>>()?;
    let dimension = subspace.dimension();

    let mut transformed_a = vec![0.0_f64; dimension * dimension];
    let mut transformed_g = vec![0.0_f64; dimension * dimension];
    for row in 0..dimension {
        for col in row..dimension {
            let a_value = bilinear_parent_form(
                &coefficients[row],
                &coefficients[col],
                parent_dimension,
                |i, j| {
                    parent
                        .pairing()
                        .entry(i, j)
                        .expect("Bernstein coefficient index is inside parent pairing matrix")
                },
            );
            let g_value = bilinear_parent_form(
                &coefficients[row],
                &coefficients[col],
                parent_dimension,
                |i, j| {
                    parent
                        .gram_entry(i, j)
                        .expect("Bernstein coefficient index is inside parent Gram matrix")
                },
            );
            checked_finite("transformed Bernstein pairing entry", a_value)?;
            checked_finite("transformed Bernstein Gram entry", g_value)?;
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
    for coefficient_vector in &coefficients {
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
        checked_finite("Bernstein subspace +1/2 boundary moment", plus_half)?;
        checked_finite("Bernstein subspace -1/2 boundary moment", minus_half)?;
        max_boundary_residual = max_boundary_residual.max(plus_half.abs().max(minus_half.abs()));
    }

    let solved = solve_transformed_pair(&transformed_a, &transformed_g, dimension)?;
    let leading_control = parent.principal_spectrum(dimension)?;

    Ok(FiniteWeilBernsteinSubspaceAudit {
        degree: subspace.degree(),
        indices: subspace.indices().to_vec(),
        parent_dimension,
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
    parent_dimension: usize,
    mut entry: impl FnMut(usize, usize) -> f64,
) -> f64 {
    let mut total = 0.0_f64;
    for i in 0..parent_dimension {
        for j in 0..parent_dimension {
            total += left[i] * entry(i, j) * right[j];
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
) -> Result<TransformedSpectrum, FiniteWeilBernsteinSubspaceError> {
    let raw_matrix = Mat::from_fn(dimension, dimension, |i, j| a[i * dimension + j]);
    let raw_decomposition = SelfAdjointEigen::new(raw_matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilBernsteinSubspaceError::RawDecompositionFailed)?;
    let raw_diagonal = raw_decomposition.S().column_vector();
    let mut raw_eigenvalues = (0..dimension)
        .map(|index| raw_diagonal[index])
        .collect::<Vec<_>>();
    for &value in &raw_eigenvalues {
        checked_finite("Bernstein subspace raw eigenvalue", value)?;
    }
    raw_eigenvalues.sort_by(f64::total_cmp);

    let gram_matrix = Mat::from_fn(dimension, dimension, |i, j| g[i * dimension + j]);
    let gram_decomposition = SelfAdjointEigen::new(gram_matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilBernsteinSubspaceError::GramDecompositionFailed)?;
    let gram_diagonal = gram_decomposition.S().column_vector();
    let gram_vectors = gram_decomposition.U();
    let raw_gram_eigenvalues = (0..dimension)
        .map(|index| gram_diagonal[index])
        .collect::<Vec<_>>();
    for &value in &raw_gram_eigenvalues {
        checked_finite("Bernstein subspace Gram eigenvalue", value)?;
    }
    let minimum_gram = raw_gram_eigenvalues
        .iter()
        .copied()
        .min_by(f64::total_cmp)
        .ok_or(FiniteWeilBernsteinSubspaceError::GramDecompositionFailed)?;
    if minimum_gram <= 0.0 {
        return Err(FiniteWeilBernsteinSubspaceError::GramNotPositiveDefinite {
            minimum_eigenvalue: minimum_gram,
        });
    }
    let maximum_gram = raw_gram_eigenvalues
        .iter()
        .copied()
        .max_by(f64::total_cmp)
        .ok_or(FiniteWeilBernsteinSubspaceError::GramDecompositionFailed)?;
    let gram_condition_number = maximum_gram / minimum_gram;
    checked_finite(
        "Bernstein subspace Gram condition number",
        gram_condition_number,
    )?;

    let mut inverse_sqrt = vec![0.0_f64; dimension * dimension];
    for i in 0..dimension {
        for j in 0..dimension {
            let mut sum = 0.0_f64;
            for k in 0..dimension {
                sum += gram_vectors[(i, k)] * gram_vectors[(j, k)] / raw_gram_eigenvalues[k].sqrt();
            }
            checked_finite("Bernstein subspace Gram inverse square root", sum)?;
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
    checked_finite(
        "Bernstein subspace whitened asymmetry",
        max_whitened_asymmetry,
    )?;

    let whitened_matrix = Mat::from_fn(dimension, dimension, |i, j| whitened[i * dimension + j]);
    let normalized_decomposition = SelfAdjointEigen::new(whitened_matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilBernsteinSubspaceError::NormalizedDecompositionFailed)?;
    let normalized_diagonal = normalized_decomposition.S().column_vector();
    let mut generalized_eigenvalues = (0..dimension)
        .map(|index| normalized_diagonal[index])
        .collect::<Vec<_>>();
    for &value in &generalized_eigenvalues {
        checked_finite("Bernstein subspace generalized eigenvalue", value)?;
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

fn checked_finite(stage: &'static str, value: f64) -> Result<(), FiniteWeilBernsteinSubspaceError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilBernsteinSubspaceError::NonFiniteEvaluation { stage, value })
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
    fn index_subspace_rejects_invalid_selections() {
        assert!(matches!(
            BernsteinIndexSubspace::new(4, Vec::new()),
            Err(FiniteWeilBernsteinSubspaceError::EmptyIndexSet)
        ));
        assert!(matches!(
            BernsteinIndexSubspace::new(4, vec![0, 2, 2]),
            Err(FiniteWeilBernsteinSubspaceError::IndicesNotStrictlyIncreasing { .. })
        ));
        assert!(matches!(
            BernsteinIndexSubspace::new(4, vec![0, 5]),
            Err(FiniteWeilBernsteinSubspaceError::IndexOutOfRange { .. })
        ));
    }

    #[test]
    fn complete_degree_two_bernstein_basis_preserves_generalized_spectrum() {
        let selected = BernsteinIndexSubspace::new(2, vec![0, 1, 2]).unwrap();
        let audit =
            audit_finite_weil_bernstein_subspace(bump(), &selected, 32, 32, 48, 48).unwrap();
        let parent = audit_finite_weil_generalized_spectrum(bump(), 3, 32, 32, 48, 48).unwrap();

        assert_eq!(audit.dimension(), 3);
        assert_eq!(audit.generalized_eigenvalues().len(), 3);
        for (&left, &right) in audit
            .generalized_eigenvalues()
            .iter()
            .zip(parent.generalized_eigenvalues().iter())
        {
            assert!((left - right).abs() <= 2.0e-10 * left.abs().max(right.abs()).max(1.0));
        }
    }

    #[test]
    fn selected_localized_subspace_is_audited_without_sign_assumption() {
        let selected = BernsteinIndexSubspace::new(4, vec![0, 2, 4]).unwrap();
        let audit =
            audit_finite_weil_bernstein_subspace(bump(), &selected, 28, 28, 40, 40).unwrap();

        assert_eq!(audit.degree(), 4);
        assert_eq!(audit.indices(), &[0, 2, 4]);
        assert_eq!(audit.dimension(), 3);
        assert_eq!(audit.parent_dimension(), 5);
        assert_eq!(audit.raw_eigenvalues().len(), 3);
        assert_eq!(audit.gram_eigenvalues().len(), 3);
        assert_eq!(audit.generalized_eigenvalues().len(), 3);
        assert!(audit.gram_eigenvalues().iter().all(|value| *value > 0.0));
        assert!(audit.gram_condition_number().is_finite());
        assert!(audit.minimum_raw_eigenvalue().is_finite());
        assert!(audit.minimum_generalized_eigenvalue().is_finite());
        assert!(audit.max_boundary_residual().is_finite());
        assert!(audit.max_whitened_asymmetry().is_finite());
        assert!(audit.leading_legendre_generalized_minimum().is_finite());
        assert!(audit.leading_legendre_gram_condition_number().is_finite());
    }
}
