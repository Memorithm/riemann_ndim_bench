//! Finite manufactured bridge from exact `Q_S` orbit descent to the source
//! semilocal `E`-sum index `M_S`.
//!
//! For an elementary local product `prod_p 1_{p^{k_p} Z_p}`, every explicit
//! non-zero sample `q = u m` is checked in two independent coordinates:
//!
//! 1. directly, by evaluating the original local balls on the diagonal `q`;
//! 2. after exact unit transport, by evaluating the transported local balls on
//!    the monoid representative `m`.
//!
//! The two exact indicators must agree term by term. The transported
//! contributions are then accumulated by `m` and passed through
//! [`SemilocalPoissonMonoid::finite_e_sum`]. This is only a finite manufactured
//! regression; it does not construct a general Bruhat--Schwartz function on
//! `A_S`, the quotient `X_S`, or the Hilbert-space Poisson map `E`.

use std::collections::BTreeMap;
use std::fmt;

use crate::semilocal_factorizable_poisson::LocalBallSpec;
use crate::semilocal_orbit_grouping::{OrbitGroupingError, group_qs_samples_by_m};
use crate::semilocal_padic_fourier::{PadicBall, PadicFourierError};
use crate::semilocal_poisson::{FiniteESum, SemilocalPoissonError, SemilocalPoissonMonoid};
use crate::semilocal_qs::{QsArithmeticError, QsRational};
use crate::semilocal_trace_contract::FinitePlaceSet;

/// Exact elementary finite factor `prod_p 1_{p^{k_p} Z_p}` for one complete
/// declared finite place set.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FiniteLocalBallProduct {
    places: FinitePlaceSet,
    specs: Vec<LocalBallSpec>,
    balls: Vec<PadicBall>,
}

impl FiniteLocalBallProduct {
    /// Build a complete local product, requiring exactly one ball for every
    /// finite place in `places`.
    pub fn new(places: FinitePlaceSet, specs: &[LocalBallSpec]) -> Result<Self, FiniteBruhatError> {
        if specs.len() != places.finite_primes().len() {
            return Err(FiniteBruhatError::IncompleteLocalProduct);
        }

        let mut sorted = specs.to_vec();
        sorted.sort_unstable_by_key(|spec| spec.prime());
        let mut balls = Vec::with_capacity(sorted.len());

        for (&expected_prime, spec) in places.finite_primes().iter().zip(sorted.iter()) {
            if spec.prime() != expected_prime {
                return Err(FiniteBruhatError::PlaceSetMismatch {
                    expected_prime,
                    actual_prime: spec.prime(),
                });
            }
            balls.push(PadicBall::new(spec.prime(), spec.exponent(), &places)?);
        }

        Ok(Self {
            places,
            specs: sorted,
            balls,
        })
    }

    #[inline]
    pub fn places(&self) -> &FinitePlaceSet {
        &self.places
    }

    #[inline]
    pub fn specs(&self) -> &[LocalBallSpec] {
        &self.specs
    }

    /// Exact diagonal evaluation of `prod_p 1_{p^{k_p} Z_p}`.
    pub fn evaluate_diagonal(&self, sample: QsRational) -> u8 {
        self.balls.iter().all(|ball| ball.contains_diagonal(sample)) as u8
    }

    /// Evaluate the product on `m in M_S` without coercing `m` through a signed
    /// scalar representation. Since `m` is prime to every finite place,
    /// `v_p(m)=0`, so membership is exactly `0 >= k_p` at each place.
    pub fn evaluate_monoid_representative(&self, m: u64) -> Result<u8, FiniteBruhatError> {
        if m == 0
            || self
                .places
                .finite_primes()
                .iter()
                .any(|&prime| m.is_multiple_of(prime))
        {
            return Err(FiniteBruhatError::NotMonoidElement { value: m });
        }

        Ok(self.balls.iter().all(|ball| ball.exponent() <= 0) as u8)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FiniteBruhatError {
    IncompleteLocalProduct,
    PlaceSetMismatch {
        expected_prime: u64,
        actual_prime: u64,
    },
    NotMonoidElement {
        value: u64,
    },
    Padic(PadicFourierError),
}

impl fmt::Display for FiniteBruhatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IncompleteLocalProduct => write!(
                f,
                "finite Bruhat product must contain exactly one ball for every finite place"
            ),
            Self::PlaceSetMismatch {
                expected_prime,
                actual_prime,
            } => write!(
                f,
                "finite Bruhat product place mismatch: expected p={expected_prime}, got p={actual_prime}"
            ),
            Self::NotMonoidElement { value } => {
                write!(f, "value is not an element of the declared M_S: {value}")
            }
            Self::Padic(error) => write!(f, "p-adic local ball construction failed: {error}"),
        }
    }
}

impl std::error::Error for FiniteBruhatError {}

impl From<PadicFourierError> for FiniteBruhatError {
    fn from(value: PadicFourierError) -> Self {
        Self::Padic(value)
    }
}

/// One represented `m in M_S` in the finite manufactured bridge.
#[derive(Clone, Debug, PartialEq)]
pub struct BruhatEGroupAudit {
    monoid_representative: u64,
    member_count: usize,
    archimedean_value: f64,
    direct_indicator_sum: usize,
    transported_indicator_sum: usize,
    grouped_raw_contribution: f64,
}

impl BruhatEGroupAudit {
    #[inline]
    pub fn monoid_representative(&self) -> u64 {
        self.monoid_representative
    }

    #[inline]
    pub fn member_count(&self) -> usize {
        self.member_count
    }

    #[inline]
    pub fn archimedean_value(&self) -> f64 {
        self.archimedean_value
    }

    #[inline]
    pub fn direct_indicator_sum(&self) -> usize {
        self.direct_indicator_sum
    }

    #[inline]
    pub fn transported_indicator_sum(&self) -> usize {
        self.transported_indicator_sum
    }

    #[inline]
    pub fn grouped_raw_contribution(&self) -> f64 {
        self.grouped_raw_contribution
    }
}

/// Audit joining direct finite `Q_S` evaluation, exact orbit descent, and the
/// finite callback-driven source `E` sum.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteBruhatEBridgeAudit {
    groups: Vec<BruhatEGroupAudit>,
    modulus: f64,
    direct_raw_total: f64,
    grouped_raw_total: f64,
    e_sum: FiniteESum,
}

impl FiniteBruhatEBridgeAudit {
    #[inline]
    pub fn groups(&self) -> &[BruhatEGroupAudit] {
        &self.groups
    }

    #[inline]
    pub fn modulus(&self) -> f64 {
        self.modulus
    }

    #[inline]
    pub fn direct_raw_total(&self) -> f64 {
        self.direct_raw_total
    }

    #[inline]
    pub fn grouped_raw_total(&self) -> f64 {
        self.grouped_raw_total
    }

    #[inline]
    pub fn raw_regrouping_residual(&self) -> f64 {
        self.direct_raw_total - self.grouped_raw_total
    }

    #[inline]
    pub fn finite_e_sum(&self) -> FiniteESum {
        self.e_sum
    }

    #[inline]
    pub fn direct_half_density_value(&self) -> f64 {
        self.modulus.sqrt() * self.direct_raw_total
    }

    #[inline]
    pub fn e_bridge_residual(&self) -> f64 {
        self.direct_half_density_value() - self.e_sum.value()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BruhatEBridgeError {
    Grouping(OrbitGroupingError),
    Qs(QsArithmeticError),
    Bruhat(FiniteBruhatError),
    Poisson(SemilocalPoissonError),
    MaxMonoidBoundTooSmall {
        max_m: u64,
        required_m: u64,
    },
    InvalidArchimedeanValue {
        monoid_representative: u64,
        value: f64,
    },
    InconsistentOrbitData {
        monoid_representative: u64,
    },
    LocalTransportMismatch {
        monoid_representative: u64,
        direct_indicator: u8,
        transported_indicator: u8,
    },
}

impl fmt::Display for BruhatEBridgeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Grouping(error) => write!(f, "Q_S orbit grouping failed: {error}"),
            Self::Qs(error) => write!(f, "Q_S arithmetic failed: {error}"),
            Self::Bruhat(error) => write!(f, "finite Bruhat evaluation failed: {error}"),
            Self::Poisson(error) => write!(f, "finite E sum failed: {error}"),
            Self::MaxMonoidBoundTooSmall { max_m, required_m } => write!(
                f,
                "finite E bound max_m={max_m} omits represented monoid element m={required_m}"
            ),
            Self::InvalidArchimedeanValue {
                monoid_representative,
                value,
            } => write!(
                f,
                "manufactured archimedean value must be finite for m={monoid_representative}: {value}"
            ),
            Self::InconsistentOrbitData {
                monoid_representative,
            } => write!(
                f,
                "independent Q_S decomposition disagrees with grouped orbit data for m={monoid_representative}"
            ),
            Self::LocalTransportMismatch {
                monoid_representative,
                direct_indicator,
                transported_indicator,
            } => write!(
                f,
                "local Bruhat indicator changed under exact unit descent for m={monoid_representative}: direct={direct_indicator}, transported={transported_indicator}"
            ),
        }
    }
}

impl std::error::Error for BruhatEBridgeError {}

impl From<OrbitGroupingError> for BruhatEBridgeError {
    fn from(value: OrbitGroupingError) -> Self {
        Self::Grouping(value)
    }
}

impl From<QsArithmeticError> for BruhatEBridgeError {
    fn from(value: QsArithmeticError) -> Self {
        Self::Qs(value)
    }
}

impl From<FiniteBruhatError> for BruhatEBridgeError {
    fn from(value: FiniteBruhatError) -> Self {
        Self::Bruhat(value)
    }
}

impl From<SemilocalPoissonError> for BruhatEBridgeError {
    fn from(value: SemilocalPoissonError) -> Self {
        Self::Poisson(value)
    }
}

/// Compare a finite explicit diagonal `Q_S` fixture with the same data after
/// exact unit descent and enumeration through `SemilocalPoissonMonoid`.
///
/// `archimedean_at_m` is evaluated exactly once for every represented `m` and
/// defines the manufactured archimedean profile on that finite support. Values
/// at unrepresented `m <= max_m` are defined to be zero for this fixture.
pub fn compare_finite_bruhat_e_bridge(
    samples: &[QsRational],
    places: &FinitePlaceSet,
    original_balls: &[LocalBallSpec],
    modulus: f64,
    max_m: u64,
    mut archimedean_at_m: impl FnMut(u64) -> f64,
) -> Result<FiniteBruhatEBridgeAudit, BruhatEBridgeError> {
    let groups = group_qs_samples_by_m(samples, places, original_balls)?;
    let original_product = FiniteLocalBallProduct::new(places.clone(), original_balls)?;
    let monoid = SemilocalPoissonMonoid::new(places.clone());

    if let Some(required_m) = groups.last().map(|group| group.monoid_representative())
        && required_m > max_m
    {
        return Err(BruhatEBridgeError::MaxMonoidBoundTooSmall { max_m, required_m });
    }

    let mut archimedean_values = BTreeMap::new();
    for group in &groups {
        let m = group.monoid_representative();
        let value = archimedean_at_m(m);
        if !value.is_finite() {
            return Err(BruhatEBridgeError::InvalidArchimedeanValue {
                monoid_representative: m,
                value,
            });
        }
        archimedean_values.insert(m, value);
    }

    // Direct side: preserve the caller's original Q_S sample order and
    // independently recover each representative from exact Q_S arithmetic.
    let mut direct_raw_total = 0.0_f64;
    for &sample in samples {
        let decomposition = sample.unit_monoid_decomposition(places)?;
        let m = decomposition.monoid_element();
        let archimedean_value =
            *archimedean_values
                .get(&m)
                .ok_or(BruhatEBridgeError::InconsistentOrbitData {
                    monoid_representative: m,
                })?;
        direct_raw_total +=
            archimedean_value * f64::from(original_product.evaluate_diagonal(sample));
    }

    // Grouped side: evaluate the transported local coordinates on m itself.
    let mut group_audits = Vec::with_capacity(groups.len());
    let mut grouped_raw_total = 0.0_f64;
    let mut grouped_contributions = BTreeMap::new();

    for group in groups {
        let m = group.monoid_representative();
        let archimedean_value = archimedean_values[&m];
        let mut direct_indicator_sum = 0_usize;
        let mut transported_indicator_sum = 0_usize;

        for term in group.terms() {
            let direct_indicator = original_product.evaluate_diagonal(term.sample());
            let transported_product =
                FiniteLocalBallProduct::new(places.clone(), term.transported_balls())?;
            let transported_indicator = transported_product.evaluate_monoid_representative(m)?;

            if direct_indicator != transported_indicator {
                return Err(BruhatEBridgeError::LocalTransportMismatch {
                    monoid_representative: m,
                    direct_indicator,
                    transported_indicator,
                });
            }

            direct_indicator_sum += usize::from(direct_indicator);
            transported_indicator_sum += usize::from(transported_indicator);
        }

        let grouped_raw_contribution = archimedean_value * transported_indicator_sum as f64;
        grouped_raw_total += grouped_raw_contribution;
        grouped_contributions.insert(m, grouped_raw_contribution);
        group_audits.push(BruhatEGroupAudit {
            monoid_representative: m,
            member_count: group.len(),
            archimedean_value,
            direct_indicator_sum,
            transported_indicator_sum,
            grouped_raw_contribution,
        });
    }

    let e_sum = monoid.finite_e_sum(modulus, max_m, |m| {
        grouped_contributions.get(&m).copied().unwrap_or(0.0)
    })?;

    Ok(FiniteBruhatEBridgeAudit {
        groups: group_audits,
        modulus,
        direct_raw_total,
        grouped_raw_total,
        e_sum,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn local_unit_descent_changes_coordinates_without_changing_membership() {
        let places = FinitePlaceSet::new(vec![2]).unwrap();
        let balls = [LocalBallSpec::new(2, 0)];
        let samples = [
            QsRational::new(3, 1, &places).unwrap(),
            QsRational::new(6, 1, &places).unwrap(),
            QsRational::new(3, 2, &places).unwrap(),
            QsRational::new(-12, 1, &places).unwrap(),
            QsRational::new(5, 1, &places).unwrap(),
            QsRational::new(5, 4, &places).unwrap(),
        ];

        let audit =
            compare_finite_bruhat_e_bridge(&samples, &places, &balls, 4.0, 5, |m| match m {
                3 => 2.0,
                5 => 4.0,
                _ => unreachable!(),
            })
            .unwrap();

        assert_eq!(audit.groups().len(), 2);
        assert_eq!(audit.groups()[0].monoid_representative(), 3);
        assert_eq!(audit.groups()[0].member_count(), 4);
        assert_eq!(audit.groups()[0].direct_indicator_sum(), 3);
        assert_eq!(audit.groups()[0].transported_indicator_sum(), 3);
        assert_eq!(audit.groups()[1].monoid_representative(), 5);
        assert_eq!(audit.groups()[1].member_count(), 2);
        assert_eq!(audit.groups()[1].direct_indicator_sum(), 1);
        assert_eq!(audit.groups()[1].transported_indicator_sum(), 1);

        assert_eq!(audit.direct_raw_total(), 10.0);
        assert_eq!(audit.grouped_raw_total(), 10.0);
        assert_eq!(audit.raw_regrouping_residual(), 0.0);
        assert_eq!(audit.finite_e_sum().raw_sum(), 10.0);
        assert_eq!(audit.finite_e_sum().value(), 20.0);
        assert_eq!(audit.direct_half_density_value(), 20.0);
        assert_eq!(audit.e_bridge_residual(), 0.0);
        // M_{ {2} } through 5 is {1,3,5}; the unrepresented m=1 term is
        // deliberately zero in this manufactured finite support.
        assert_eq!(audit.finite_e_sum().term_count(), 3);
    }

    #[test]
    fn multi_prime_fixture_preserves_local_membership_after_transport() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let balls = [LocalBallSpec::new(2, 0), LocalBallSpec::new(3, 0)];
        let samples = [
            QsRational::new(5, 1, &places).unwrap(),
            QsRational::new(10, 1, &places).unwrap(),
            QsRational::new(15, 1, &places).unwrap(),
            QsRational::new(5, 2, &places).unwrap(),
            QsRational::new(5, 3, &places).unwrap(),
            QsRational::new(45, 8, &places).unwrap(),
            QsRational::new(7, 1, &places).unwrap(),
            QsRational::new(14, 1, &places).unwrap(),
            QsRational::new(7, 3, &places).unwrap(),
        ];

        let audit =
            compare_finite_bruhat_e_bridge(&samples, &places, &balls, 9.0, 7, |m| 1.0 / m as f64)
                .unwrap();

        assert_eq!(audit.groups().len(), 2);
        assert_eq!(audit.groups()[0].monoid_representative(), 5);
        assert_eq!(audit.groups()[0].direct_indicator_sum(), 3);
        assert_eq!(audit.groups()[0].transported_indicator_sum(), 3);
        assert_eq!(audit.groups()[1].monoid_representative(), 7);
        assert_eq!(audit.groups()[1].direct_indicator_sum(), 2);
        assert_eq!(audit.groups()[1].transported_indicator_sum(), 2);

        let scale = audit
            .direct_half_density_value()
            .abs()
            .max(audit.finite_e_sum().value().abs())
            .max(1.0);
        assert!(audit.raw_regrouping_residual().abs() <= 2.0e-15 * scale);
        assert!(audit.e_bridge_residual().abs() <= 2.0e-15 * scale);
    }

    #[test]
    fn finite_e_bound_must_cover_every_represented_m() {
        let places = FinitePlaceSet::new(vec![2]).unwrap();
        let sample = QsRational::new(5, 1, &places).unwrap();
        let error = compare_finite_bruhat_e_bridge(
            &[sample],
            &places,
            &[LocalBallSpec::new(2, 0)],
            1.0,
            3,
            |_| 1.0,
        )
        .unwrap_err();

        assert_eq!(
            error,
            BruhatEBridgeError::MaxMonoidBoundTooSmall {
                max_m: 3,
                required_m: 5,
            }
        );
    }

    #[test]
    fn non_finite_archimedean_profile_is_rejected() {
        let places = FinitePlaceSet::new(vec![2]).unwrap();
        let sample = QsRational::new(3, 1, &places).unwrap();
        let error = compare_finite_bruhat_e_bridge(
            &[sample],
            &places,
            &[LocalBallSpec::new(2, 0)],
            1.0,
            3,
            |_| f64::NAN,
        )
        .unwrap_err();

        assert!(matches!(
            error,
            BruhatEBridgeError::InvalidArchimedeanValue {
                monoid_representative: 3,
                value,
            } if value.is_nan()
        ));
    }
}
