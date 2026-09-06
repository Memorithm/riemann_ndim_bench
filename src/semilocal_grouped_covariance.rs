//! Manufactured finite `Q_S -> M_S` regrouping audit built from the exact
//! orbit grouping and certified unit-covariance fixtures.
//!
//! The ungrouped side evaluates the original factorizable lattice contribution
//! for explicit `Q_S` samples in input order. The grouped side independently
//! re-evaluates each sample after its exact `Q_S^*` unit transport, organizes the
//! terms by the unique `m in M_S` representative, and applies one explicit
//! caller-supplied weight per representative.
//!
//! This is a finite manufactured regression. It does not justify rearranging an
//! infinite `Q_S` sum, construct `X_S`, or establish the Hilbert-space Poisson
//! map `E`.

use std::collections::BTreeMap;
use std::fmt;

use crate::semilocal_factorizable_poisson::LocalBallSpec;
use crate::semilocal_orbit_grouping::{OrbitGroupingError, group_qs_samples_by_m};
use crate::semilocal_qs::QsRational;
use crate::semilocal_trace_contract::FinitePlaceSet;
use crate::semilocal_unit_covariance::{UnitCovarianceError, compare_unit_covariance};

/// One `M_S`-group contribution in the manufactured covariance audit.
#[derive(Clone, Debug, PartialEq)]
pub struct GroupedCovarianceGroupAudit {
    monoid_representative: u64,
    weight: f64,
    member_count: usize,
    direct_total: f64,
    transported_total: f64,
    combined_tail_bound: f64,
}

impl GroupedCovarianceGroupAudit {
    #[inline]
    pub fn monoid_representative(&self) -> u64 {
        self.monoid_representative
    }

    #[inline]
    pub fn weight(&self) -> f64 {
        self.weight
    }

    #[inline]
    pub fn member_count(&self) -> usize {
        self.member_count
    }

    #[inline]
    pub fn direct_total(&self) -> f64 {
        self.direct_total
    }

    #[inline]
    pub fn transported_total(&self) -> f64 {
        self.transported_total
    }

    #[inline]
    pub fn residual(&self) -> f64 {
        self.direct_total - self.transported_total
    }

    #[inline]
    pub fn combined_tail_bound(&self) -> f64 {
        self.combined_tail_bound
    }
}

/// Auditable finite comparison between input-order direct evaluation and
/// `M_S`-grouped unit-transported evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGroupedCovarianceAudit {
    groups: Vec<GroupedCovarianceGroupAudit>,
    ungrouped_direct_total: f64,
    grouped_transported_total: f64,
    combined_tail_bound: f64,
}

impl FiniteGroupedCovarianceAudit {
    #[inline]
    pub fn groups(&self) -> &[GroupedCovarianceGroupAudit] {
        &self.groups
    }

    #[inline]
    pub fn ungrouped_direct_total(&self) -> f64 {
        self.ungrouped_direct_total
    }

    #[inline]
    pub fn grouped_transported_total(&self) -> f64 {
        self.grouped_transported_total
    }

    #[inline]
    pub fn residual(&self) -> f64 {
        self.ungrouped_direct_total - self.grouped_transported_total
    }

    /// Rigorous Gaussian truncation contribution from both independently
    /// evaluated sides. Floating-point roundoff is intentionally excluded.
    #[inline]
    pub fn combined_tail_bound(&self) -> f64 {
        self.combined_tail_bound
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GroupedCovarianceError {
    Grouping(OrbitGroupingError),
    Covariance(UnitCovarianceError),
    InvalidWeight {
        monoid_representative: u64,
        weight: f64,
    },
    InconsistentOrbitData {
        monoid_representative: u64,
    },
}

impl fmt::Display for GroupedCovarianceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Grouping(error) => write!(f, "Q_S orbit grouping failed: {error}"),
            Self::Covariance(error) => write!(f, "Q_S unit covariance audit failed: {error}"),
            Self::InvalidWeight {
                monoid_representative,
                weight,
            } => write!(
                f,
                "manufactured group weight must be finite for m={monoid_representative}: {weight}"
            ),
            Self::InconsistentOrbitData {
                monoid_representative,
            } => write!(
                f,
                "independent unit covariance data disagrees with exact grouping for m={monoid_representative}"
            ),
        }
    }
}

impl std::error::Error for GroupedCovarianceError {}

impl From<OrbitGroupingError> for GroupedCovarianceError {
    fn from(value: OrbitGroupingError) -> Self {
        Self::Grouping(value)
    }
}

impl From<UnitCovarianceError> for GroupedCovarianceError {
    fn from(value: UnitCovarianceError) -> Self {
        Self::Covariance(value)
    }
}

/// Compare a finite direct `Q_S` total with the independently transported total
/// after exact grouping by `m in M_S`.
///
/// `weight_for_m` is evaluated exactly once for every represented `m`; this
/// keeps the grouping coefficient explicit and prevents a stateful callback
/// from assigning different weights to the two sides of the comparison.
pub fn compare_grouped_unit_covariance(
    samples: &[QsRational],
    places: &FinitePlaceSet,
    original_balls: &[LocalBallSpec],
    max_abs_n: u64,
    mut weight_for_m: impl FnMut(u64) -> f64,
) -> Result<FiniteGroupedCovarianceAudit, GroupedCovarianceError> {
    let exact_groups = group_qs_samples_by_m(samples, places, original_balls)?;
    let mut weights = BTreeMap::new();

    for group in &exact_groups {
        let representative = group.monoid_representative();
        let weight = weight_for_m(representative);
        if !weight.is_finite() {
            return Err(GroupedCovarianceError::InvalidWeight {
                monoid_representative: representative,
                weight,
            });
        }
        weights.insert(representative, weight);
    }

    let mut ungrouped_direct_total = 0.0_f64;
    let mut ungrouped_direct_tail = 0.0_f64;
    for &sample in samples {
        let comparison = compare_unit_covariance(sample, places, original_balls, max_abs_n)?;
        let representative = comparison.monoid_representative();
        let weight = *weights
            .get(&representative)
            .ok_or(GroupedCovarianceError::InconsistentOrbitData {
                monoid_representative: representative,
            })?;
        ungrouped_direct_total += weight * comparison.original().original_sum().value();
        ungrouped_direct_tail +=
            weight.abs() * comparison.original().original_sum().absolute_tail_bound();
    }

    let mut group_audits = Vec::with_capacity(exact_groups.len());
    let mut grouped_transported_total = 0.0_f64;
    let mut grouped_transported_tail = 0.0_f64;

    for group in exact_groups {
        let representative = group.monoid_representative();
        let weight = weights[&representative];
        let mut direct_total = 0.0_f64;
        let mut transported_total = 0.0_f64;
        let mut combined_tail_bound = 0.0_f64;

        for term in group.terms() {
            let comparison =
                compare_unit_covariance(term.sample(), places, original_balls, max_abs_n)?;
            if comparison.monoid_representative() != representative
                || comparison.transported_balls() != term.transported_balls()
            {
                return Err(GroupedCovarianceError::InconsistentOrbitData {
                    monoid_representative: representative,
                });
            }

            direct_total += weight * comparison.original().original_sum().value();
            transported_total += weight * comparison.rescaled_transported_sum().value();
            combined_tail_bound += weight.abs() * comparison.combined_tail_bound();
            grouped_transported_tail +=
                weight.abs() * comparison.rescaled_transported_sum().absolute_tail_bound();
        }

        grouped_transported_total += transported_total;
        group_audits.push(GroupedCovarianceGroupAudit {
            monoid_representative: representative,
            weight,
            member_count: group.len(),
            direct_total,
            transported_total,
            combined_tail_bound,
        });
    }

    Ok(FiniteGroupedCovarianceAudit {
        groups: group_audits,
        ungrouped_direct_total,
        grouped_transported_total,
        combined_tail_bound: ungrouped_direct_tail + grouped_transported_tail,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_audit_close(audit: &FiniteGroupedCovarianceAudit, relative_roundoff: f64) {
        let scale = audit
            .ungrouped_direct_total()
            .abs()
            .max(audit.grouped_transported_total().abs())
            .max(1.0);
        let roundoff = relative_roundoff * scale;
        assert!(
            audit.residual().abs() <= audit.combined_tail_bound() + roundoff,
            "residual={:.3e} tails={:.3e} roundoff={:.3e}",
            audit.residual(),
            audit.combined_tail_bound(),
            roundoff
        );
    }

    #[test]
    fn dyadic_grouped_total_matches_input_order_direct_total() {
        let places = FinitePlaceSet::new(vec![2]).unwrap();
        let balls = [LocalBallSpec::new(2, 1)];
        let samples = [
            QsRational::new(3, 1, &places).unwrap(),
            QsRational::new(6, 1, &places).unwrap(),
            QsRational::new(3, 2, &places).unwrap(),
            QsRational::new(3, 8, &places).unwrap(),
            QsRational::new(-12, 1, &places).unwrap(),
            QsRational::new(5, 1, &places).unwrap(),
            QsRational::new(20, 1, &places).unwrap(),
            QsRational::new(5, 4, &places).unwrap(),
        ];

        let audit = compare_grouped_unit_covariance(
            &samples,
            &places,
            &balls,
            256,
            |m| match m {
                3 => 1.0,
                5 => 2.0,
                _ => unreachable!(),
            },
        )
        .unwrap();

        assert_eq!(audit.groups().len(), 2);
        assert_eq!(audit.groups()[0].monoid_representative(), 3);
        assert_eq!(audit.groups()[0].member_count(), 5);
        assert_eq!(audit.groups()[1].monoid_representative(), 5);
        assert_eq!(audit.groups()[1].member_count(), 3);
        assert_audit_close(&audit, 2.0e-15);
    }

    #[test]
    fn multi_prime_grouped_total_survives_float_audit() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let balls = [LocalBallSpec::new(2, 1), LocalBallSpec::new(3, -2)];
        let samples = [
            QsRational::new(5, 1, &places).unwrap(),
            QsRational::new(10, 1, &places).unwrap(),
            QsRational::new(45, 8, &places).unwrap(),
            QsRational::new(-5, 1, &places).unwrap(),
            QsRational::new(7, 1, &places).unwrap(),
            QsRational::new(14, 1, &places).unwrap(),
        ];

        let audit = compare_grouped_unit_covariance(
            &samples,
            &places,
            &balls,
            256,
            |m| 1.0 / m as f64,
        )
        .unwrap();

        assert_eq!(audit.groups().len(), 2);
        assert_audit_close(&audit, 1.0e-12);
    }

    #[test]
    fn non_finite_group_weight_is_rejected() {
        let places = FinitePlaceSet::new(vec![2]).unwrap();
        let sample = QsRational::new(3, 2, &places).unwrap();
        let error = compare_grouped_unit_covariance(
            &[sample],
            &places,
            &[LocalBallSpec::new(2, 0)],
            64,
            |_| f64::NAN,
        )
        .unwrap_err();

        assert!(matches!(
            error,
            GroupedCovarianceError::InvalidWeight {
                monoid_representative: 3,
                weight
            } if weight.is_nan()
        ));
    }
}
