//! Conditional Schur-complement test for finite Weil pairing blocks.
//!
//! Decompose a real-symmetric pairing as M = [A B; B^T D]. If
//! A >= alpha I for alpha > 0, D >= delta I, and ||B||_2 <= beta,
//! then the Schur complement satisfies
//!
//! D - B^T A^{-1} B >= (delta - beta^2 / alpha) I.
//!
//! Thus alpha > 0 and delta - beta^2/alpha > 0 suffice for strict
//! positivity of M. The lemma also applies to bounded cross-couplings
//! with an infinite-dimensional tail when ALL operator assumptions are
//! independently established. The finite adapter below only checks a
//! finite leading/trailing matrix partition, not any infinite tail.
//!
//! All matrix-entry error budgets are *caller-supplied assumptions*.
//! Binary64 calculations here have no directed rounding and cannot be
//! promoted to validated interval-arithmetic or RH certificates.

use std::fmt;

use crate::weil_quadratic_matrix::FiniteWeilQuadraticMatrixAudit;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConditionalSchurVerdict {
    StrictlyPositiveUnderAssumptions,
    Undetermined,
}

/// General conditional Schur lemma independent of any specific matrix basis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConditionalSchurMargin {
    head_lower_bound: f64,
    tail_lower_bound: f64,
    coupling_norm_upper_bound: f64,
    tail_schur_lower_bound: Option<f64>,
    verdict: ConditionalSchurVerdict,
}

impl ConditionalSchurMargin {
    pub const fn head_lower_bound(self) -> f64 {
        self.head_lower_bound
    }

    pub const fn tail_lower_bound(self) -> f64 {
        self.tail_lower_bound
    }

    pub const fn coupling_norm_upper_bound(self) -> f64 {
        self.coupling_norm_upper_bound
    }

    /// Returns None if the available head lower bound is not positive.
    pub const fn tail_schur_lower_bound(self) -> Option<f64> {
        self.tail_schur_lower_bound
    }

    pub const fn verdict(self) -> ConditionalSchurVerdict {
        self.verdict
    }
}

/// Result for one declared leading/trailing split of a finite Weil matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteWeilSchurAudit {
    dimension: usize,
    head_dimension: usize,
    margin: ConditionalSchurMargin,
}

impl FiniteWeilSchurAudit {
    pub const fn dimension(self) -> usize {
        self.dimension
    }

    pub const fn head_dimension(self) -> usize {
        self.head_dimension
    }

    pub const fn tail_dimension(self) -> usize {
        self.dimension - self.head_dimension
    }

    pub const fn margin(self) -> ConditionalSchurMargin {
        self.margin
    }

    pub const fn verdict(self) -> ConditionalSchurVerdict {
        self.margin.verdict()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteWeilSchurError {
    EmptyMatrix,
    DimensionOverflow,
    InvalidSplit {
        head_dimension: usize,
        available: usize,
    },
    InvalidLength {
        field: &'static str,
        expected: usize,
        actual: usize,
    },
    InvalidValue {
        field: &'static str,
        index: usize,
        value: f64,
    },
    AsymmetricMatrix {
        row: usize,
        col: usize,
    },
    NonFiniteComputation {
        stage: &'static str,
    },
}

impl fmt::Display for FiniteWeilSchurError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for FiniteWeilSchurError {}

/// Conditional operator-level Schur implication.
///
/// * head_lower_bound: independently justified lower bound on A
/// * tail_lower_bound: independently justified lower bound on D
/// * coupling_norm_upper_bound: independently justified operator norm bound
///
/// This evaluates a finite-precision arithmetic *diagnostic*. It does NOT
/// establish that the caller's inputs are valid bounds, nor does it prove
/// that a global or semilocal Weil form admits this block decomposition.
pub fn audit_conditional_schur_margin(
    head_lower_bound: f64,
    tail_lower_bound: f64,
    coupling_norm_upper_bound: f64,
) -> Result<ConditionalSchurMargin, FiniteWeilSchurError> {
    for (index, value) in [head_lower_bound, tail_lower_bound, coupling_norm_upper_bound]
        .into_iter()
        .enumerate()
    {
        if !value.is_finite() || (index == 2 && value < 0.0) {
            return Err(FiniteWeilSchurError::InvalidValue {
                field: "conditional Schur margin",
                index,
                value,
            });
        }
    }
    let tail_schur_lower_bound = if head_lower_bound > 0.0 {
        let coupling_square = coupling_norm_upper_bound * coupling_norm_upper_bound;
        checked_finite("coupling norm squared", coupling_square)?;
        let reduction = coupling_square / head_lower_bound;
        checked_finite("Schur reduction", reduction)?;
        let lower = tail_lower_bound - reduction;
        checked_finite("Schur tail lower bound", lower)?;
        Some(lower)
    } else {
        None
    };

    let verdict = if tail_schur_lower_bound.is_some_and(|lower| lower > 0.0) {
        ConditionalSchurVerdict::StrictlyPositiveUnderAssumptions
    } else {
        ConditionalSchurVerdict::Undetermined
    };

    Ok(ConditionalSchurMargin {
        head_lower_bound,
        tail_lower_bound,
        coupling_norm_upper_bound,
        tail_schur_lower_bound,
        verdict,
    })
}

/// Evaluate a block Schur condition on the existing compact Weil basis.
///
/// Every entrywise budget E_ij must be an independently justified upper bound
/// on |M_exact(i,j) - M_numeric(i,j)|. It must cover the complete quadrature,
/// approximation, and floating-point error. A quadrature-refinement spread is
/// not such a bound. The returned verdict always remains CONDITIONAL.
pub fn audit_finite_weil_schur_gate(
    pairing: &FiniteWeilQuadraticMatrixAudit,
    entrywise_error_bounds: &[f64],
    head_dimension: usize,
) -> Result<FiniteWeilSchurAudit, FiniteWeilSchurError> {
    let n = pairing.dimension();
    let size = n
        .checked_mul(n)
        .ok_or(FiniteWeilSchurError::DimensionOverflow)?;
    let mut matrix = Vec::with_capacity(size);
    for i in 0..n {
        for j in 0..n {
            matrix.push(
                pairing
                    .entry(i, j)
                    .expect("finite Weil pairing indices are in range"),
            );
        }
    }
    audit_dense_schur_gate(n, &matrix, entrywise_error_bounds, head_dimension)
}

fn audit_dense_schur_gate(
    n: usize,
    matrix: &[f64],
    errors: &[f64],
    head_dimension: usize,
) -> Result<FiniteWeilSchurAudit, FiniteWeilSchurError> {
    if n == 0 {
        return Err(FiniteWeilSchurError::EmptyMatrix);
    }
    if head_dimension == 0 || head_dimension >= n {
        return Err(FiniteWeilSchurError::InvalidSplit {
            head_dimension,
            available: n,
        });
    }
    let count = n
        .checked_mul(n)
        .ok_or(FiniteWeilSchurError::DimensionOverflow)?;
    for (field, actual, expected) in [
        ("matrix", matrix.len(), count),
        ("entrywise_error_bounds", errors.len(), count),
    ] {
        if actual != expected {
            return Err(FiniteWeilSchurError::InvalidLength {
                field,
                expected,
                actual,
            });
        }
    }
    for (field, values, nonnegative) in [
        ("matrix", matrix, false),
        ("entrywise_error_bounds", errors, true),
    ] {
        for (index, &value) in values.iter().enumerate() {
            if !value.is_finite() || (nonnegative && value < 0.0) {
                return Err(FiniteWeilSchurError::InvalidValue {
                    field,
                    index,
                    value,
                });
            }
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if matrix[i * n + j] != matrix[j * n + i] {
                return Err(FiniteWeilSchurError::AsymmetricMatrix { row: i, col: j });
            }
        }
    }

    let head_lower = block_gershgorin_lower(n, matrix, errors, 0, head_dimension)?;
    let tail_lower = block_gershgorin_lower(n, matrix, errors, head_dimension, n)?;
    // Frobenius bounds spectral norm without any eigensolver. Each symmetric
    // cross entry is enclosed using the more conservative of its two budgets.
    let mut coupling_square = 0.0;
    for i in 0..head_dimension {
        for j in head_dimension..n {
            let error = errors[i * n + j].max(errors[j * n + i]);
            let magnitude = matrix[i * n + j].abs() + error;
            coupling_square += magnitude * magnitude;
        }
    }
    checked_finite("coupling Frobenius norm squared", coupling_square)?;
    let coupling_upper = coupling_square.sqrt();
    checked_finite("coupling Frobenius norm", coupling_upper)?;

    let margin = audit_conditional_schur_margin(
        head_lower,
        tail_lower,
        coupling_upper,
    )?;
    Ok(FiniteWeilSchurAudit {
        dimension: n,
        head_dimension,
        margin,
    })
}

fn block_gershgorin_lower(
    n: usize,
    matrix: &[f64],
    errors: &[f64],
    start: usize,
    end: usize,
) -> Result<f64, FiniteWeilSchurError> {
    let mut minimum = f64::INFINITY;
    for i in start..end {
        let mut radius = 0.0;
        for j in start..end {
            if j != i {
                let error = errors[i * n + j].max(errors[j * n + i]);
                radius += matrix[i * n + j].abs() + error;
            }
        }
        let lower = matrix[i * n + i] - errors[i * n + i] - radius;
        checked_finite("block Gershgorin lower bound", lower)?;
        minimum = minimum.min(lower);
    }
    Ok(minimum)
}

fn checked_finite(
    stage: &'static str,
    value: f64,
) -> Result<(), FiniteWeilSchurError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilSchurError::NonFiniteComputation { stage })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::{
        CompactArchimedeanBump, PositiveRational,
    };
    use crate::weil_quadratic_matrix::audit_finite_weil_quadratic_matrix;

    #[test]
    fn block_schur_passes_when_whole_matrix_is_not_diagonally_dominant() {
        // The final row fails whole-matrix Gershgorin: 3 - 1.6 - 1.6 < 0.
        // The leading 2x2 is stable enough to absorb the cross coupling.
        let a = [2.0, 0.2, 1.6, 0.2, 2.0, 1.6, 1.6, 1.6, 3.0];
        let result = audit_dense_schur_gate(3, &a, &[0.01; 9], 2).unwrap();
        assert_eq!(result.dimension(), 3);
        assert_eq!(result.tail_dimension(), 1);
        assert_eq!(
            result.verdict(),
            ConditionalSchurVerdict::StrictlyPositiveUnderAssumptions
        );
        assert!(result.margin().head_lower_bound() > 1.7);
        assert!(result.margin().tail_schur_lower_bound().unwrap() > 0.0);
    }

    #[test]
    fn too_wide_error_budget_is_undetermined() {
        let a = [2.0, 0.2, 1.6, 0.2, 2.0, 1.6, 1.6, 1.6, 3.0];
        let result = audit_dense_schur_gate(3, &a, &[0.1; 9], 2).unwrap();
        assert_eq!(result.verdict(), ConditionalSchurVerdict::Undetermined);
    }

    #[test]
    fn zero_coupling_and_strict_positive_blocks_are_positive() {
        let a = [0.8, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 2.0];
        let audit = audit_dense_schur_gate(3, &a, &[0.01; 9], 1).unwrap();
        assert_eq!(audit.verdict(), ConditionalSchurVerdict::StrictlyPositiveUnderAssumptions);
        assert!(audit.margin().tail_schur_lower_bound().unwrap() > 0.0);
    }

    #[test]
    fn failed_head_pivot_does_not_invert_a_nonpositive_lower_bound() {
        let audit = audit_conditional_schur_margin(-1.0, 2.0, 0.1).unwrap();
        assert_eq!(audit.verdict(), ConditionalSchurVerdict::Undetermined);
        assert_eq!(audit.tail_schur_lower_bound(), None);
    }

    #[test]
    fn conditional_operator_tail_contract_is_explicit() {
        let audit = audit_conditional_schur_margin(2.0, 3.0, 1.0).unwrap();
        assert_eq!(audit.tail_schur_lower_bound(), Some(2.5));
        assert_eq!(
            audit.verdict(),
            ConditionalSchurVerdict::StrictlyPositiveUnderAssumptions
        );
    }

    #[test]
    fn invalid_data_are_rejected_before_sign_decision() {
        let a = [1.0, 0.0, 0.0, 1.0];
        assert!(matches!(
            audit_dense_schur_gate(0, &[], &[], 0),
            Err(FiniteWeilSchurError::EmptyMatrix)
        ));
        assert!(matches!(
            audit_dense_schur_gate(2, &a, &[0.0; 4], 2),
            Err(FiniteWeilSchurError::InvalidSplit { .. })
        ));
        assert!(matches!(
            audit_dense_schur_gate(2, &a, &[0.0; 3], 1),
            Err(FiniteWeilSchurError::InvalidLength { .. })
        ));
        assert!(matches!(
            audit_dense_schur_gate(2, &a, &[-0.1, 0.0, 0.0, 0.0], 1),
            Err(FiniteWeilSchurError::InvalidValue { .. })
        ));
        assert!(matches!(
            audit_dense_schur_gate(2, &[1.0, 0.5, 0.1, 1.0], &[0.0; 4], 1),
            Err(FiniteWeilSchurError::AsymmetricMatrix { .. })
        ));
        assert!(matches!(
            audit_conditional_schur_margin(1.0, 1.0, -1.0),
            Err(FiniteWeilSchurError::InvalidValue { .. })
        ));
        assert!(matches!(
            audit_conditional_schur_margin(1.0, 1.0, f64::MAX),
            Err(FiniteWeilSchurError::NonFiniteComputation { .. })
        ));
    }

    #[test]
    fn existing_weil_matrix_is_accepted_without_claiming_error_certification() {
        let bump = CompactArchimedeanBump::new(
            PositiveRational::new(1, 2).unwrap(),
            PositiveRational::new(7, 2).unwrap(),
        )
        .unwrap();
        let matrix = audit_finite_weil_quadratic_matrix(bump, 2, 28, 28, 40)
            .unwrap();
        // E=0 here probes only the adapter; it is NOT a certified error budget.
        let audit = audit_finite_weil_schur_gate(&matrix, &[0.0; 4], 1).unwrap();
        assert_eq!(audit.head_dimension(), 1);
        assert!(audit.margin().coupling_norm_upper_bound().is_finite());
    }
}
