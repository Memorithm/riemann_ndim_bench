//! Conditional finite-span sign gate for the compact Riemann--Weil pairing.
//!
//! The pairing is evaluated numerically in the basis of
//! weil_quadratic_matrix. This module accepts *externally justified*
//! entrywise absolute error bounds and returns two useful implications:
//!
//! - strict Gershgorin dominance gives positivity on the declared finite span;
//! - a strictly negative upper bound on a Rayleigh numerator gives a
//!   finite-dimensional negative direction.
//!
//! Both implications are CONDITIONAL on the supplied bounds enclosing the
//! exact real-symmetric Weil pairing entries, including quadrature and
//! floating-point errors. Merely observing quadrature convergence does not
//! establish those bounds. Binary64 arithmetic here is not outward-rounded
//! interval arithmetic, so the output is not a formal certificate.
//! Neither outcome establishes a complete-space statement or RH.

use std::fmt;

use crate::weil_quadratic_matrix::FiniteWeilQuadraticMatrixAudit;

/// What the finite calculation would imply if the supplied entry bounds are valid.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConditionalFiniteWeilSign {
    /// The exact pairing would be positive definite on this finite span.
    StrictlyPositiveFiniteSpan,
    /// The declared nonzero coefficient direction would have negative pairing.
    StrictlyNegativeWitness,
    /// The available bounds decide neither question.
    Undetermined,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteWeilSignGateAudit {
    dimension: usize,
    direction_norm_squared: f64,
    direction_rayleigh: f64,
    direction_rayleigh_error: f64,
    gershgorin_lower_bound: f64,
    verdict: ConditionalFiniteWeilSign,
}

impl FiniteWeilSignGateAudit {
    pub const fn dimension(self) -> usize {
        self.dimension
    }

    pub const fn direction_norm_squared(self) -> f64 {
        self.direction_norm_squared
    }

    /// Approximate Rayleigh quotient using the symmetrized finite Weil matrix.
    pub const fn direction_rayleigh(self) -> f64 {
        self.direction_rayleigh
    }

    /// Bound on Rayleigh quotient error, conditional on valid entry bounds.
    pub const fn direction_rayleigh_error(self) -> f64 {
        self.direction_rayleigh_error
    }

    /// Conservative lower bound on the smallest eigenvalue of the exact
    /// finite pairing, conditional on valid entrywise bounds.
    pub const fn gershgorin_lower_bound(self) -> f64 {
        self.gershgorin_lower_bound
    }

    pub const fn verdict(self) -> ConditionalFiniteWeilSign {
        self.verdict
    }

    pub fn direction_interval(self) -> (f64, f64) {
        (
            self.direction_rayleigh - self.direction_rayleigh_error,
            self.direction_rayleigh + self.direction_rayleigh_error,
        )
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteWeilSignGateError {
    EmptyMatrix,
    MatrixDimensionOverflow,
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
    ZeroDirection,
    NonFiniteComputation {
        stage: &'static str,
    },
}

impl fmt::Display for FiniteWeilSignGateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for FiniteWeilSignGateError {}

/// Evaluate finite-span sign margins using the existing compact Riemann--Weil
/// quadratic matrix.
///
/// The caller MUST provide a valid absolute bound for EVERY matrix entry.
/// It must include all integration, model, and rounding errors. The bounds
/// are an explicit assumption and are not inferred from numerical refinement.
pub fn audit_finite_weil_sign_gate(
    matrix: &FiniteWeilQuadraticMatrixAudit,
    entrywise_error_bounds: &[f64],
    direction: &[f64],
) -> Result<FiniteWeilSignGateAudit, FiniteWeilSignGateError> {
    let n = matrix.dimension();
    let size = n
        .checked_mul(n)
        .ok_or(FiniteWeilSignGateError::MatrixDimensionOverflow)?;
    let mut values = Vec::with_capacity(size);
    for i in 0..n {
        for j in 0..n {
            values.push(
                matrix
                    .entry(i, j)
                    .expect("indices are inside finite Weil pairing"),
            );
        }
    }
    evaluate_dense_sign_gate(n, &values, entrywise_error_bounds, direction)
}

fn evaluate_dense_sign_gate(
    n: usize,
    matrix: &[f64],
    errors: &[f64],
    direction: &[f64],
) -> Result<FiniteWeilSignGateAudit, FiniteWeilSignGateError> {
    if n == 0 {
        return Err(FiniteWeilSignGateError::EmptyMatrix);
    }
    let size = n
        .checked_mul(n)
        .ok_or(FiniteWeilSignGateError::MatrixDimensionOverflow)?;
    for (field, actual, expected) in [
        ("matrix", matrix.len(), size),
        ("entrywise_error_bounds", errors.len(), size),
        ("direction", direction.len(), n),
    ] {
        if actual != expected {
            return Err(FiniteWeilSignGateError::InvalidLength {
                field,
                expected,
                actual,
            });
        }
    }

    for (field, values, nonnegative) in [
        ("matrix", matrix, false),
        ("entrywise_error_bounds", errors, true),
        ("direction", direction, false),
    ] {
        for (index, &value) in values.iter().enumerate() {
            if !value.is_finite() || (nonnegative && value < 0.0) {
                return Err(FiniteWeilSignGateError::InvalidValue {
                    field,
                    index,
                    value,
                });
            }
        }
    }

    let norm_squared = direction.iter().map(|v| v * v).sum::<f64>();
    finite("direction norm", norm_squared)?;
    if norm_squared == 0.0 {
        return Err(FiniteWeilSignGateError::ZeroDirection);
    }

    // A source Weil pairing is real-symmetric. Do not silently symmetrize
    // an arbitrary, asymmetric matrix supplied to this calculation.
    for i in 0..n {
        for j in (i + 1)..n {
            if matrix[i * n + j] != matrix[j * n + i] {
                return Err(FiniteWeilSignGateError::AsymmetricMatrix { row: i, col: j });
            }
        }
    }

    let mut gershgorin_lower = f64::INFINITY;
    let mut rayleigh_numerator = 0.0;
    let mut rayleigh_uncertainty_numerator = 0.0;
    for i in 0..n {
        let mut row_radius = 0.0;
        for j in 0..n {
            // An exact symmetric entry is enclosed by both declared bounds.
            // Taking the larger of the two makes mismatched budgets safe.
            let error = errors[i * n + j].max(errors[j * n + i]);
            rayleigh_numerator += direction[i] * matrix[i * n + j] * direction[j];
            rayleigh_uncertainty_numerator += direction[i].abs() * error * direction[j].abs();
            if i != j {
                row_radius += matrix[i * n + j].abs() + error;
            }
        }
        let row_lower = matrix[i * n + i] - errors[i * n + i] - row_radius;
        finite("Gershgorin row lower bound", row_lower)?;
        gershgorin_lower = gershgorin_lower.min(row_lower);
    }

    finite("Rayleigh numerator", rayleigh_numerator)?;
    finite("Rayleigh error numerator", rayleigh_uncertainty_numerator)?;
    let direction_rayleigh = rayleigh_numerator / norm_squared;
    let direction_rayleigh_error = rayleigh_uncertainty_numerator / norm_squared;
    finite("Rayleigh quotient", direction_rayleigh)?;
    finite("Rayleigh error", direction_rayleigh_error)?;

    let (direction_lower, direction_upper) = (
        direction_rayleigh - direction_rayleigh_error,
        direction_rayleigh + direction_rayleigh_error,
    );
    finite("Rayleigh lower bound", direction_lower)?;
    finite("Rayleigh upper bound", direction_upper)?;

    let verdict = if gershgorin_lower > 0.0 {
        ConditionalFiniteWeilSign::StrictlyPositiveFiniteSpan
    } else if direction_upper < 0.0 {
        ConditionalFiniteWeilSign::StrictlyNegativeWitness
    } else {
        ConditionalFiniteWeilSign::Undetermined
    };

    Ok(FiniteWeilSignGateAudit {
        dimension: n,
        direction_norm_squared: norm_squared,
        direction_rayleigh,
        direction_rayleigh_error,
        gershgorin_lower_bound: gershgorin_lower,
        verdict,
    })
}

fn finite(stage: &'static str, value: f64) -> Result<(), FiniteWeilSignGateError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilSignGateError::NonFiniteComputation { stage })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::{CompactArchimedeanBump, PositiveRational};
    use crate::weil_quadratic_matrix::audit_finite_weil_quadratic_matrix;

    #[test]
    fn strictly_positive_finite_span_with_entry_uncertainty() {
        let a = [2.0, 0.25, 0.25, 1.0];
        let errors = [0.1, 0.05, 0.05, 0.1];
        let audit = evaluate_dense_sign_gate(2, &a, &errors, &[1.0, -1.0]).unwrap();
        assert_eq!(
            audit.verdict(),
            ConditionalFiniteWeilSign::StrictlyPositiveFiniteSpan
        );
        assert!(audit.gershgorin_lower_bound() >= 0.6 - 1.0e-14);
        assert!(audit.direction_interval().0 > 0.0);
    }

    #[test]
    fn negative_witness_survives_declared_error_budget() {
        let a = [-1.0, 0.0, 0.0, 2.0];
        let errors = [0.1, 0.0, 0.0, 0.1];
        let audit = evaluate_dense_sign_gate(2, &a, &errors, &[3.0, 0.0]).unwrap();
        assert_eq!(
            audit.verdict(),
            ConditionalFiniteWeilSign::StrictlyNegativeWitness
        );
        assert!(audit.direction_interval().1 < 0.0);
        assert!((audit.direction_rayleigh() + 1.0).abs() < 1.0e-14);
    }

    #[test]
    fn undecided_sign_remains_undetermined() {
        let a = [0.0, 1.0, 1.0, 0.0];
        let audit = evaluate_dense_sign_gate(2, &a, &[0.2, 0.2, 0.2, 0.2], &[1.0, 1.0]).unwrap();
        assert_eq!(audit.verdict(), ConditionalFiniteWeilSign::Undetermined);
    }

    #[test]
    fn inconsistent_and_nonfinite_inputs_fail_closed() {
        let a = [1.0, 0.0, 0.0, 1.0];
        assert!(matches!(
            evaluate_dense_sign_gate(2, &a, &[0.0; 4], &[0.0, 0.0]),
            Err(FiniteWeilSignGateError::ZeroDirection)
        ));
        assert!(matches!(
            evaluate_dense_sign_gate(2, &a, &[0.0; 3], &[1.0, 0.0]),
            Err(FiniteWeilSignGateError::InvalidLength { .. })
        ));
        assert!(matches!(
            evaluate_dense_sign_gate(2, &a, &[0.0, -0.1, 0.0, 0.0], &[1.0, 0.0]),
            Err(FiniteWeilSignGateError::InvalidValue { .. })
        ));
        assert!(matches!(
            evaluate_dense_sign_gate(2, &[1.0, 0.5, 0.25, 1.0], &[0.0; 4], &[1.0, 0.0]),
            Err(FiniteWeilSignGateError::AsymmetricMatrix { .. })
        ));
        assert!(matches!(
            evaluate_dense_sign_gate(2, &a, &[0.0; 4], &[f64::NAN, 0.0]),
            Err(FiniteWeilSignGateError::InvalidValue { .. })
        ));
        assert!(matches!(
            evaluate_dense_sign_gate(2, &a, &[0.0; 4], &[f64::MAX, f64::MAX]),
            Err(FiniteWeilSignGateError::NonFiniteComputation { .. })
        ));
    }

    #[test]
    fn existing_weil_matrix_feeds_the_gate_without_auto_claims() {
        let bump = CompactArchimedeanBump::new(
            PositiveRational::new(1, 2).unwrap(),
            PositiveRational::new(7, 2).unwrap(),
        )
        .unwrap();
        let matrix = audit_finite_weil_quadratic_matrix(bump, 1, 28, 28, 40).unwrap();
        // Zero errors here are a *synthetic input* for API integration, not
        // an assertion that the numerical quadrature is error-free.
        let audit = audit_finite_weil_sign_gate(&matrix, &[0.0], &[1.0]).unwrap();
        assert_eq!(audit.dimension(), 1);
        assert!(audit.direction_rayleigh().is_finite());
    }
}
