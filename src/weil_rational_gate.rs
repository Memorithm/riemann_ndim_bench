//! Exact integer-arithmetic checks for finite Riemann--Weil sign obligations.
//!
//! Each matrix entry is enclosed in a rational interval [lo/S, hi/S]
//! where S is the SAME positive integer for every entry. The caller must
//! establish these enclosures independently (e.g. from validated quadrature).
//! All comparisons here use checked i128 operations; an overflow is an error,
//! never an accepted sign. No floating-point conversion occurs.
//!
//! For symmetric M with these bounds:
//! * Gershgorin proves strict positivity if every diagonal lower endpoint
//!   exceeds the sum of absolute off-diagonal upper magnitudes.
//! * A nonzero integer vector v proves a negative direction if the EXACT upper
//!   endpoint of v^T M v is negative.
//! * For a declared leading/trailing split, the block Schur criterion uses
//!   alpha*delta > ||B||_F^2. Both sides have common denominator S^2,
//!   so the comparison uses only exact integer multiplication.
//!
//! These are rigorously checked algebraic implications CONDITIONAL on the
//! input intervals containing the exact compact Weil pairing entries. They
//! neither generate validated integration enclosures nor address the missing
//! infinite-dimensional semilocal-to-Weil bridge or RH.

use std::fmt;

/// One closed interval with endpoints expressed in common units of 1/scale.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WeilRationalInterval {
    lower: i128,
    upper: i128,
}

impl WeilRationalInterval {
    pub fn new(lower: i128, upper: i128) -> Result<Self, ExactWeilSignError> {
        if lower > upper {
            return Err(ExactWeilSignError::ReversedInterval { lower, upper });
        }
        Ok(Self { lower, upper })
    }

    pub const fn lower(self) -> i128 {
        self.lower
    }

    pub const fn upper(self) -> i128 {
        self.upper
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ExactFiniteWeilSign {
    /// A strictly positive finite form, conditional on valid input enclosures.
    StrictlyPositiveFiniteSpan,
    /// A strictly negative quadratic direction on the finite span.
    StrictlyNegativeWitness,
    /// Neither strict inequality is established.
    Undetermined,
}

/// Exact numerators. Gershgorin and Rayleigh values are in 1/scale units.
/// The Schur discriminant is in 1/scale^2 units.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ExactWeilSignAudit {
    dimension: usize,
    scale: i128,
    gershgorin_lower_numerator: i128,
    rayleigh_upper_numerator: i128,
    head_dimension: Option<usize>,
    head_lower_numerator: Option<i128>,
    tail_lower_numerator: Option<i128>,
    coupling_frobenius_squared_numerator: Option<i128>,
    schur_discriminant_numerator: Option<i128>,
    verdict: ExactFiniteWeilSign,
}

impl ExactWeilSignAudit {
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    pub const fn scale(&self) -> i128 {
        self.scale
    }

    pub const fn gershgorin_lower_numerator(&self) -> i128 {
        self.gershgorin_lower_numerator
    }

    /// Upper bound for v^T M v, without division by v^T v.
    pub const fn rayleigh_upper_numerator(&self) -> i128 {
        self.rayleigh_upper_numerator
    }

    pub const fn head_dimension(&self) -> Option<usize> {
        self.head_dimension
    }

    pub const fn head_lower_numerator(&self) -> Option<i128> {
        self.head_lower_numerator
    }

    pub const fn tail_lower_numerator(&self) -> Option<i128> {
        self.tail_lower_numerator
    }

    pub const fn coupling_frobenius_squared_numerator(&self) -> Option<i128> {
        self.coupling_frobenius_squared_numerator
    }

    /// Some(alpha*delta - beta_squared), only if alpha > 0.
    /// A strictly positive numerator implies positive Schur complement.
    pub const fn schur_discriminant_numerator(&self) -> Option<i128> {
        self.schur_discriminant_numerator
    }

    pub const fn verdict(&self) -> ExactFiniteWeilSign {
        self.verdict
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ExactWeilSignError {
    EmptyMatrix,
    InvalidScale { scale: i128 },
    DimensionOverflow,
    InvalidLength {
        field: &'static str,
        expected: usize,
        actual: usize,
    },
    ReversedInterval { lower: i128, upper: i128 },
    AsymmetricInterval { row: usize, col: usize },
    InvalidSplit { head: usize, dimension: usize },
    ZeroDirection,
    ArithmeticOverflow { stage: &'static str },
}

impl fmt::Display for ExactWeilSignError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for ExactWeilSignError {}

/// Check interval finite-form sign statements without floating-point arithmetic.
///
/// Entries are row-major and MUST be symmetric as closed intervals.
/// No interval is inferred from a floating-point center or an empirical
/// quadrature difference. A caller that lacks genuine enclosure data
/// must not interpret the result as a Weil pairing certificate.
///
/// The optional split is a FINITE leading/trailing partition, NOT a claim
/// about an infinite Sonin complement.
pub fn audit_exact_weil_interval_sign(
    dimension: usize,
    scale: i128,
    matrix: &[WeilRationalInterval],
    witness: &[i128],
    head_dimension: Option<usize>,
) -> Result<ExactWeilSignAudit, ExactWeilSignError> {
    if dimension == 0 {
        return Err(ExactWeilSignError::EmptyMatrix);
    }
    if scale <= 0 {
        return Err(ExactWeilSignError::InvalidScale { scale });
    }
    let entries = dimension
        .checked_mul(dimension)
        .ok_or(ExactWeilSignError::DimensionOverflow)?;
    for (field, actual, expected) in [
        ("matrix", matrix.len(), entries),
        ("witness", witness.len(), dimension),
    ] {
        if actual != expected {
            return Err(ExactWeilSignError::InvalidLength {
                field,
                expected,
                actual,
            });
        }
    }
    if let Some(head) = head_dimension
        && (head == 0 || head >= dimension)
    {
        return Err(ExactWeilSignError::InvalidSplit { head, dimension });
    }
    if witness.iter().all(|&value| value == 0) {
        return Err(ExactWeilSignError::ZeroDirection);
    }
    for i in 0..dimension {
        for j in (i + 1)..dimension {
            if matrix[i * dimension + j] != matrix[j * dimension + i] {
                return Err(ExactWeilSignError::AsymmetricInterval { row: i, col: j });
            }
        }
    }

    let gershgorin_lower_numerator =
        gershgorin_block_lower(dimension, matrix, 0, dimension)?;

    // For each term v_i * M_ij * v_j, choose the upper interval endpoint
    // according to the exact sign of the integer coefficient v_i * v_j.
    let mut rayleigh_upper_numerator = 0_i128;
    for i in 0..dimension {
        for j in 0..dimension {
            let coefficient = mul(witness[i], witness[j], "Rayleigh coefficient")?;
            let endpoint = if coefficient >= 0 {
                matrix[i * dimension + j].upper
            } else {
                matrix[i * dimension + j].lower
            };
            let term = mul(coefficient, endpoint, "Rayleigh term")?;
            rayleigh_upper_numerator =
                add(rayleigh_upper_numerator, term, "Rayleigh upper sum")?;
        }
    }

    let mut head_lower_numerator = None;
    let mut tail_lower_numerator = None;
    let mut coupling_frobenius_squared_numerator = None;
    let mut schur_discriminant_numerator = None;

    if let Some(head) = head_dimension {
        let head_lower = gershgorin_block_lower(dimension, matrix, 0, head)?;
        let tail_lower =
            gershgorin_block_lower(dimension, matrix, head, dimension)?;
        let mut coupling_squared = 0_i128;
        for i in 0..head {
            for j in head..dimension {
                let coupling = max_abs(matrix[i * dimension + j])?;
                let square = mul(coupling, coupling, "cross-block coupling square")?;
                coupling_squared =
                    add(coupling_squared, square, "cross-block Frobenius square")?;
            }
        }
        head_lower_numerator = Some(head_lower);
        tail_lower_numerator = Some(tail_lower);
        coupling_frobenius_squared_numerator = Some(coupling_squared);
        if head_lower > 0 {
            // Since all denominators are S^2>0, compare the integers
            // alpha_num * delta_num and beta_squared_num exactly.
            let product = mul(head_lower, tail_lower, "Schur block product")?;
            schur_discriminant_numerator =
                Some(sub(product, coupling_squared, "Schur discriminant")?);
        }
    }

    let positive_by_schur =
        schur_discriminant_numerator.is_some_and(|discriminant| discriminant > 0);
    let verdict = if gershgorin_lower_numerator > 0 || positive_by_schur {
        ExactFiniteWeilSign::StrictlyPositiveFiniteSpan
    } else if rayleigh_upper_numerator < 0 {
        ExactFiniteWeilSign::StrictlyNegativeWitness
    } else {
        ExactFiniteWeilSign::Undetermined
    };

    Ok(ExactWeilSignAudit {
        dimension,
        scale,
        gershgorin_lower_numerator,
        rayleigh_upper_numerator,
        head_dimension,
        head_lower_numerator,
        tail_lower_numerator,
        coupling_frobenius_squared_numerator,
        schur_discriminant_numerator,
        verdict,
    })
}

fn gershgorin_block_lower(
    dimension: usize,
    matrix: &[WeilRationalInterval],
    start: usize,
    end: usize,
) -> Result<i128, ExactWeilSignError> {
    let mut lower = i128::MAX;
    for i in start..end {
        let mut radius = 0_i128;
        for j in start..end {
            if i != j {
                radius = add(
                    radius,
                    max_abs(matrix[i * dimension + j])?,
                    "Gershgorin radius",
                )?;
            }
        }
        let row = sub(
            matrix[i * dimension + i].lower,
            radius,
            "Gershgorin row",
        )?;
        lower = lower.min(row);
    }
    Ok(lower)
}

fn max_abs(interval: WeilRationalInterval) -> Result<i128, ExactWeilSignError> {
    let a = interval
        .lower
        .checked_abs()
        .ok_or(ExactWeilSignError::ArithmeticOverflow {
            stage: "interval lower magnitude",
        })?;
    let b = interval
        .upper
        .checked_abs()
        .ok_or(ExactWeilSignError::ArithmeticOverflow {
            stage: "interval upper magnitude",
        })?;
    Ok(a.max(b))
}

fn add(left: i128, right: i128, stage: &'static str) -> Result<i128, ExactWeilSignError> {
    left.checked_add(right)
        .ok_or(ExactWeilSignError::ArithmeticOverflow { stage })
}

fn sub(left: i128, right: i128, stage: &'static str) -> Result<i128, ExactWeilSignError> {
    left.checked_sub(right)
        .ok_or(ExactWeilSignError::ArithmeticOverflow { stage })
}

fn mul(left: i128, right: i128, stage: &'static str) -> Result<i128, ExactWeilSignError> {
    left.checked_mul(right)
        .ok_or(ExactWeilSignError::ArithmeticOverflow { stage })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interval(lo: i128, hi: i128) -> WeilRationalInterval {
        WeilRationalInterval::new(lo, hi).unwrap()
    }

    #[test]
    fn exact_gershgorin_positive_with_nonzero_intervals() {
        let matrix = [
            interval(190, 210), interval(-20, 20),
            interval(-20, 20), interval(90, 110),
        ];
        let result = audit_exact_weil_interval_sign(2, 100, &matrix, &[1, -1], None)
            .unwrap();
        assert_eq!(result.gershgorin_lower_numerator(), 70);
        assert_eq!(result.verdict(), ExactFiniteWeilSign::StrictlyPositiveFiniteSpan);
        assert!(result.rayleigh_upper_numerator() > 0);
    }

    #[test]
    fn exact_schur_can_succeed_when_full_gershgorin_fails() {
        // At scale 100, M approximately equals
        // [2.0 0.2 1.6; 0.2 2.0 1.6; 1.6 1.6 3.0].
        // Entrywise uncertainty 0.01 is encoded exactly as +/-1 unit.
        let centers = [
            200, 20, 160,
            20, 200, 160,
            160, 160, 300,
        ];
        let matrix = centers.map(|value| interval(value - 1, value + 1));
        let result = audit_exact_weil_interval_sign(3, 100, &matrix, &[1, 0, 0], Some(2))
            .unwrap();
        assert_eq!(result.gershgorin_lower_numerator(), -23);
        assert_eq!(result.head_lower_numerator(), Some(178));
        assert_eq!(result.tail_lower_numerator(), Some(299));
        assert_eq!(result.coupling_frobenius_squared_numerator(), Some(51_842));
        assert_eq!(result.schur_discriminant_numerator(), Some(1_380));
        assert_eq!(result.verdict(), ExactFiniteWeilSign::StrictlyPositiveFiniteSpan);
    }

    #[test]
    fn wide_enclosure_invalidates_schur_decision() {
        let centers = [
            200, 20, 160,
            20, 200, 160,
            160, 160, 300,
        ];
        let matrix = centers.map(|value| interval(value - 10, value + 10));
        let result = audit_exact_weil_interval_sign(3, 100, &matrix, &[1, 0, 0], Some(2))
            .unwrap();
        assert_eq!(result.verdict(), ExactFiniteWeilSign::Undetermined);
        assert!(result.schur_discriminant_numerator().unwrap() < 0);
    }

    #[test]
    fn exact_negative_rayleigh_witness_requires_strict_upper_sign() {
        let matrix = [
            interval(-101, -99), interval(0, 0),
            interval(0, 0), interval(199, 201),
        ];
        let result = audit_exact_weil_interval_sign(2, 100, &matrix, &[3, 0], Some(1))
            .unwrap();
        assert_eq!(result.rayleigh_upper_numerator(), -891);
        assert_eq!(result.verdict(), ExactFiniteWeilSign::StrictlyNegativeWitness);
        assert_eq!(result.schur_discriminant_numerator(), None);
    }

    #[test]
    fn zero_boundary_is_not_misclassified() {
        let matrix = [interval(-1, 1)];
        let result = audit_exact_weil_interval_sign(1, 1, &matrix, &[1], None).unwrap();
        assert_eq!(result.verdict(), ExactFiniteWeilSign::Undetermined);
    }

    #[test]
    fn invalid_intervals_symmetry_and_dimensions_fail_closed() {
        assert!(matches!(
            WeilRationalInterval::new(5, 4),
            Err(ExactWeilSignError::ReversedInterval { .. })
        ));
        assert!(matches!(
            audit_exact_weil_interval_sign(0, 1, &[], &[], None),
            Err(ExactWeilSignError::EmptyMatrix)
        ));
        let matrix = [
            interval(1, 1), interval(0, 0),
            interval(0, 1), interval(1, 1),
        ];
        assert!(matches!(
            audit_exact_weil_interval_sign(2, 1, &matrix, &[1, 0], None),
            Err(ExactWeilSignError::AsymmetricInterval { .. })
        ));
        assert!(matches!(
            audit_exact_weil_interval_sign(2, 0, &matrix, &[1, 0], None),
            Err(ExactWeilSignError::InvalidScale { .. })
        ));
        assert!(matches!(
            audit_exact_weil_interval_sign(2, 1, &matrix[..3], &[1, 0], None),
            Err(ExactWeilSignError::InvalidLength { .. })
        ));
        assert!(matches!(
            audit_exact_weil_interval_sign(2, 1, &matrix, &[0, 0], None),
            Err(ExactWeilSignError::ZeroDirection)
        ));
        assert!(matches!(
            audit_exact_weil_interval_sign(2, 1, &matrix, &[1, 0], Some(2)),
            Err(ExactWeilSignError::InvalidSplit { .. })
        ));
    }

    #[test]
    fn integer_overflow_is_rejected_not_promoted_to_a_sign() {
        let minimum = [interval(i128::MIN, i128::MIN)];
        assert!(matches!(
            audit_exact_weil_interval_sign(1, 1, &minimum, &[1], None),
            Err(ExactWeilSignError::ArithmeticOverflow { .. })
                | Ok(ExactWeilSignAudit {
                    verdict: ExactFiniteWeilSign::StrictlyNegativeWitness,
                    ..
                })
        ));
        // The following 2x2 deliberately overflows the exact Schur product:
        let huge = i128::MAX / 3;
        let matrix = [
            interval(huge, huge), interval(0, 0),
            interval(0, 0), interval(huge, huge),
        ];
        assert!(matches!(
            audit_exact_weil_interval_sign(2, 1, &matrix, &[1, 0], Some(1)),
            Err(ExactWeilSignError::ArithmeticOverflow { .. })
        ));
        assert!(matches!(
            audit_exact_weil_interval_sign(1, 1, &[interval(1, 1)], &[i128::MAX], None),
            Err(ExactWeilSignError::ArithmeticOverflow { .. })
        ));
    }

    #[test]
    fn exact_result_is_unchanged_under_common_scale_multiplication() {
        let a = [
            interval(20, 20), interval(-1, 1),
            interval(-1, 1), interval(10, 10),
        ];
        let b = [
            interval(200, 200), interval(-10, 10),
            interval(-10, 10), interval(100, 100),
        ];
        let small = audit_exact_weil_interval_sign(2, 10, &a, &[1, 1], Some(1))
            .unwrap();
        let large = audit_exact_weil_interval_sign(2, 100, &b, &[1, 1], Some(1))
            .unwrap();
        assert_eq!(small.verdict(), large.verdict());
        assert_eq!(small.schur_discriminant_numerator().unwrap() * 100,
                   large.schur_discriminant_numerator().unwrap());
    }
}
