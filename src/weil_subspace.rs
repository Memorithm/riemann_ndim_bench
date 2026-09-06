//! Finite Riemann--Weil sensitivity audit across explicit Legendre subspaces.
//!
//! The existing finite Weil machinery builds one validated pairing matrix `A`
//! and multiplicative Gram matrix `G` on the leading Legendre family
//!
//! `g_j(rho) = bump(rho) P_j(2t-1)`, `h_j = Q g_j`.
//!
//! This module reuses a single sufficiently large `(A,G)` computation and then
//! restricts it to explicitly declared degree sets.  This changes the finite
//! test-function subspace itself, rather than merely rescaling a fixed basis.
//! It therefore probes sensitivity to finite subspace choice while holding the
//! support, Riemann--Weil decomposition and quadrature levels fixed.
//!
//! These are finite numerical restrictions.  Agreement of signs or spectra
//! across several selected subspaces is stronger finite evidence only; it does
//! not establish density/completeness, complete-space Weil positivity, a
//! uniform approximation theorem, Conjecture 4.1, or RH.

use std::fmt;

use faer::{Mat, Side, linalg::solvers::SelfAdjointEigen};

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_generalized_spectrum::{
    FiniteWeilGeneralizedSpectrumError, audit_finite_weil_generalized_spectrum,
};

/// One finite subspace selected by strictly increasing Legendre degrees.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LegendreDegreeSubspace {
    degrees: Vec<usize>,
}

impl LegendreDegreeSubspace {
    pub fn new(degrees: Vec<usize>) -> Result<Self, FiniteWeilSubspaceError> {
        if degrees.is_empty() {
            return Err(FiniteWeilSubspaceError::EmptyDegreeSet);
        }
        for pair in degrees.windows(2) {
            if pair[0] >= pair[1] {
                return Err(FiniteWeilSubspaceError::DegreesNotStrictlyIncreasing {
                    previous: pair[0],
                    next: pair[1],
                });
            }
        }
        Ok(Self { degrees })
    }

    #[inline]
    pub fn degrees(&self) -> &[usize] {
        &self.degrees
    }

    #[inline]
    pub fn dimension(&self) -> usize {
        self.degrees.len()
    }

    #[inline]
    pub fn maximum_degree(&self) -> usize {
        *self
            .degrees
            .last()
            .expect("validated Legendre degree subspace is non-empty")
    }
}

/// Spectrum and numerical-quality diagnostics for one selected subspace.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilSubspaceSpectrum {
    degrees: Vec<usize>,
    raw_eigenvalues: Vec<f64>,
    gram_eigenvalues: Vec<f64>,
    generalized_eigenvalues: Vec<f64>,
    gram_condition_number: f64,
    max_boundary_residual: f64,
    max_whitened_asymmetry: f64,
}

impl FiniteWeilSubspaceSpectrum {
    #[inline]
    pub fn degrees(&self) -> &[usize] {
        &self.degrees
    }

    #[inline]
    pub fn dimension(&self) -> usize {
        self.degrees.len()
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
    pub fn max_boundary_residual(&self) -> f64 {
        self.max_boundary_residual
    }

    #[inline]
    pub fn max_whitened_asymmetry(&self) -> f64 {
        self.max_whitened_asymmetry
    }
}

/// One shared parent computation plus all requested finite subspace spectra.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilSubspaceAudit {
    computed_legendre_dimension: usize,
    full_matrix_max_raw_pairing_asymmetry: f64,
    spectra: Vec<FiniteWeilSubspaceSpectrum>,
}

impl FiniteWeilSubspaceAudit {
    #[inline]
    pub fn computed_legendre_dimension(&self) -> usize {
        self.computed_legendre_dimension
    }

    /// Conservative diagnostic inherited from the single full pairing matrix.
    /// It is the maximum directional pairing asymmetry over every degree in the
    /// parent computation, not only over an individual selected subspace.
    #[inline]
    pub fn full_matrix_max_raw_pairing_asymmetry(&self) -> f64 {
        self.full_matrix_max_raw_pairing_asymmetry
    }

    #[inline]
    pub fn spectra(&self) -> &[FiniteWeilSubspaceSpectrum] {
        &self.spectra
    }
}

#[derive(Debug)]
pub enum FiniteWeilSubspaceError {
    EmptySubspaceSet,
    EmptyDegreeSet,
    DegreesNotStrictlyIncreasing { previous: usize, next: usize },
    ParentDimensionOverflow { maximum_degree: usize },
    Parent(FiniteWeilGeneralizedSpectrumError),
    RawDecompositionFailed,
    GramDecompositionFailed,
    GramNotPositiveDefinite { minimum_eigenvalue: f64 },
    NormalizedDecompositionFailed,
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilSubspaceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySubspaceSet => write!(f, "finite Weil subspace set must be non-empty"),
            Self::EmptyDegreeSet => write!(f, "Legendre degree subspace must be non-empty"),
            Self::DegreesNotStrictlyIncreasing { previous, next } => write!(
                f,
                "Legendre degrees must be strictly increasing: previous={previous}, next={next}"
            ),
            Self::ParentDimensionOverflow { maximum_degree } => write!(
                f,
                "maximum Legendre degree {maximum_degree} cannot be converted to a parent dimension"
            ),
            Self::Parent(error) => write!(f, "parent generalized-spectrum audit failed: {error}"),
            Self::RawDecompositionFailed => write!(f, "raw subspace eigendecomposition failed"),
            Self::GramDecompositionFailed => write!(f, "subspace Gram eigendecomposition failed"),
            Self::GramNotPositiveDefinite { minimum_eigenvalue } => write!(
                f,
                "subspace Gram matrix is not numerically positive definite: lambda_min={minimum_eigenvalue}"
            ),
            Self::NormalizedDecompositionFailed => {
                write!(f, "normalized subspace eigendecomposition failed")
            }
            Self::NonFiniteEvaluation { stage, value } => {
                write!(f, "non-finite finite-subspace value at {stage}: {value}")
            }
        }
    }
}

impl std::error::Error for FiniteWeilSubspaceError {}

impl From<FiniteWeilGeneralizedSpectrumError> for FiniteWeilSubspaceError {
    fn from(value: FiniteWeilGeneralizedSpectrumError) -> Self {
        Self::Parent(value)
    }
}

/// Audit explicitly selected Legendre-degree subspaces while reusing one parent
/// finite Weil pairing/Gram computation up to the largest requested degree.
pub fn audit_finite_weil_legendre_subspaces(
    bump: CompactArchimedeanBump,
    subspaces: &[LegendreDegreeSubspace],
    correlation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
    gram_order: usize,
) -> Result<FiniteWeilSubspaceAudit, FiniteWeilSubspaceError> {
    if subspaces.is_empty() {
        return Err(FiniteWeilSubspaceError::EmptySubspaceSet);
    }

    let maximum_degree = subspaces
        .iter()
        .map(LegendreDegreeSubspace::maximum_degree)
        .max()
        .expect("validated non-empty subspace set has a maximum degree");
    let parent_dimension = maximum_degree
        .checked_add(1)
        .ok_or(FiniteWeilSubspaceError::ParentDimensionOverflow { maximum_degree })?;

    let parent = audit_finite_weil_generalized_spectrum(
        bump,
        parent_dimension,
        correlation_order,
        archimedean_order,
        boundary_order,
        gram_order,
    )?;

    let mut spectra = Vec::with_capacity(subspaces.len());
    for subspace in subspaces {
        let dimension = subspace.dimension();
        let mut a = vec![0.0_f64; dimension * dimension];
        let mut g = vec![0.0_f64; dimension * dimension];

        for (row, &source_row) in subspace.degrees().iter().enumerate() {
            for (col, &source_col) in subspace.degrees().iter().enumerate() {
                a[row * dimension + col] = parent
                    .pairing()
                    .entry(source_row, source_col)
                    .expect("selected degree is inside the shared parent pairing matrix");
                g[row * dimension + col] = parent
                    .gram_entry(source_row, source_col)
                    .expect("selected degree is inside the shared parent Gram matrix");
            }
        }

        let max_boundary_residual = subspace
            .degrees()
            .iter()
            .map(|&degree| parent.pairing().boundary_residuals()[degree])
            .fold(0.0_f64, f64::max);
        checked_finite("selected boundary residual", max_boundary_residual)?;

        let solved = solve_dense_subspace(&a, &g, dimension)?;
        spectra.push(FiniteWeilSubspaceSpectrum {
            degrees: subspace.degrees().to_vec(),
            raw_eigenvalues: solved.raw_eigenvalues,
            gram_eigenvalues: solved.gram_eigenvalues,
            generalized_eigenvalues: solved.generalized_eigenvalues,
            gram_condition_number: solved.gram_condition_number,
            max_boundary_residual,
            max_whitened_asymmetry: solved.max_whitened_asymmetry,
        });
    }

    Ok(FiniteWeilSubspaceAudit {
        computed_legendre_dimension: parent_dimension,
        full_matrix_max_raw_pairing_asymmetry: parent.pairing().max_raw_pairing_asymmetry(),
        spectra,
    })
}

struct DenseSubspaceSpectrum {
    raw_eigenvalues: Vec<f64>,
    gram_eigenvalues: Vec<f64>,
    generalized_eigenvalues: Vec<f64>,
    gram_condition_number: f64,
    max_whitened_asymmetry: f64,
}

fn solve_dense_subspace(
    a: &[f64],
    g: &[f64],
    dimension: usize,
) -> Result<DenseSubspaceSpectrum, FiniteWeilSubspaceError> {
    let raw_matrix = Mat::from_fn(dimension, dimension, |i, j| a[i * dimension + j]);
    let raw_decomposition = SelfAdjointEigen::new(raw_matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilSubspaceError::RawDecompositionFailed)?;
    let raw_diagonal = raw_decomposition.S().column_vector();
    let mut raw_eigenvalues = (0..dimension)
        .map(|index| raw_diagonal[index])
        .collect::<Vec<_>>();
    for &value in &raw_eigenvalues {
        checked_finite("raw subspace eigenvalue", value)?;
    }
    raw_eigenvalues.sort_by(f64::total_cmp);

    let gram_matrix = Mat::from_fn(dimension, dimension, |i, j| g[i * dimension + j]);
    let gram_decomposition = SelfAdjointEigen::new(gram_matrix.as_ref(), Side::Lower)
        .map_err(|_| FiniteWeilSubspaceError::GramDecompositionFailed)?;
    let gram_diagonal = gram_decomposition.S().column_vector();
    let gram_vectors = gram_decomposition.U();
    let raw_gram_eigenvalues = (0..dimension)
        .map(|index| gram_diagonal[index])
        .collect::<Vec<_>>();
    for &value in &raw_gram_eigenvalues {
        checked_finite("subspace Gram eigenvalue", value)?;
    }

    let minimum_gram = raw_gram_eigenvalues
        .iter()
        .copied()
        .min_by(f64::total_cmp)
        .ok_or(FiniteWeilSubspaceError::GramDecompositionFailed)?;
    if minimum_gram <= 0.0 {
        return Err(FiniteWeilSubspaceError::GramNotPositiveDefinite {
            minimum_eigenvalue: minimum_gram,
        });
    }
    let maximum_gram = raw_gram_eigenvalues
        .iter()
        .copied()
        .max_by(f64::total_cmp)
        .ok_or(FiniteWeilSubspaceError::GramDecompositionFailed)?;
    let gram_condition_number = maximum_gram / minimum_gram;
    checked_finite("subspace Gram condition number", gram_condition_number)?;

    let mut inverse_sqrt = vec![0.0_f64; dimension * dimension];
    for i in 0..dimension {
        for j in 0..dimension {
            let mut sum = 0.0_f64;
            for k in 0..dimension {
                sum += gram_vectors[(i, k)] * gram_vectors[(j, k)] / raw_gram_eigenvalues[k].sqrt();
            }
            checked_finite("subspace Gram inverse square root", sum)?;
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
    checked_finite("subspace whitened asymmetry", max_whitened_asymmetry)?;

    let normalized_matrix = Mat::from_fn(dimension, dimension, |i, j| whitened[i * dimension + j]);
    let normalized_decomposition =
        SelfAdjointEigen::new(normalized_matrix.as_ref(), Side::Lower)
            .map_err(|_| FiniteWeilSubspaceError::NormalizedDecompositionFailed)?;
    let normalized_diagonal = normalized_decomposition.S().column_vector();
    let mut generalized_eigenvalues = (0..dimension)
        .map(|index| normalized_diagonal[index])
        .collect::<Vec<_>>();
    for &value in &generalized_eigenvalues {
        checked_finite("generalized subspace eigenvalue", value)?;
    }
    generalized_eigenvalues.sort_by(f64::total_cmp);

    let mut gram_eigenvalues = raw_gram_eigenvalues;
    gram_eigenvalues.sort_by(f64::total_cmp);

    Ok(DenseSubspaceSpectrum {
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

fn checked_finite(stage: &'static str, value: f64) -> Result<(), FiniteWeilSubspaceError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilSubspaceError::NonFiniteEvaluation { stage, value })
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
    fn degree_sets_must_be_non_empty_and_strictly_increasing() {
        assert!(matches!(
            LegendreDegreeSubspace::new(Vec::new()),
            Err(FiniteWeilSubspaceError::EmptyDegreeSet)
        ));
        assert!(matches!(
            LegendreDegreeSubspace::new(vec![0, 2, 2]),
            Err(FiniteWeilSubspaceError::DegreesNotStrictlyIncreasing { .. })
        ));
        assert!(matches!(
            LegendreDegreeSubspace::new(vec![1, 0]),
            Err(FiniteWeilSubspaceError::DegreesNotStrictlyIncreasing { .. })
        ));
    }

    #[test]
    fn leading_subspace_matches_existing_principal_generalized_spectrum() {
        let selected = LegendreDegreeSubspace::new(vec![0, 1, 2]).unwrap();
        let audit =
            audit_finite_weil_legendre_subspaces(bump(), &[selected], 40, 40, 56, 56).unwrap();
        let direct = audit_finite_weil_generalized_spectrum(bump(), 3, 40, 40, 56, 56).unwrap();
        let principal = direct.principal_spectrum(3).unwrap();
        let selected = &audit.spectra()[0];

        assert_eq!(audit.computed_legendre_dimension(), 3);
        assert_eq!(selected.degrees(), &[0, 1, 2]);
        assert!(
            (selected.minimum_raw_eigenvalue() - principal.raw_minimum_eigenvalue()).abs()
                <= 2.0e-12
        );
        assert!(
            (selected.minimum_generalized_eigenvalue()
                - principal.generalized_minimum_eigenvalue())
            .abs()
                <= 2.0e-12
        );
        assert!(
            (selected.gram_condition_number() - principal.gram_condition_number()).abs()
                <= 2.0e-10 * principal.gram_condition_number().max(1.0)
        );
    }

    #[test]
    fn non_leading_even_degree_subspace_is_audited_without_sign_assumption() {
        let leading = LegendreDegreeSubspace::new(vec![0, 1, 2]).unwrap();
        let even = LegendreDegreeSubspace::new(vec![0, 2, 4]).unwrap();
        let audit =
            audit_finite_weil_legendre_subspaces(bump(), &[leading, even], 32, 32, 48, 48).unwrap();

        assert_eq!(audit.computed_legendre_dimension(), 5);
        assert_eq!(audit.spectra().len(), 2);
        for spectrum in audit.spectra() {
            assert_eq!(spectrum.dimension(), 3);
            assert_eq!(spectrum.raw_eigenvalues().len(), 3);
            assert_eq!(spectrum.gram_eigenvalues().len(), 3);
            assert_eq!(spectrum.generalized_eigenvalues().len(), 3);
            assert!(spectrum.gram_eigenvalues().iter().all(|value| *value > 0.0));
            assert!(spectrum.gram_condition_number().is_finite());
            assert!(spectrum.minimum_raw_eigenvalue().is_finite());
            assert!(spectrum.minimum_generalized_eigenvalue().is_finite());
            assert!(spectrum.max_boundary_residual().is_finite());
            assert!(spectrum.max_whitened_asymmetry().is_finite());
        }
    }
}
