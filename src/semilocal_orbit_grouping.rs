//! Exact finite grouping of non-zero `Q_S` samples by their unique
//! `Q_S^*`-orbit representative in `M_S`.
//!
//! The source decomposition immediately before equation (4.6) writes every
//! non-zero `q in Q_S` uniquely as
//!
//! `q = u m`, `u in Q_S^*`, `m in M_S`.
//!
//! This module makes that finite grouping executable while retaining the full
//! unit decomposition and the induced transport of the finite local ball
//! factors.  It deliberately does not replace an infinite `Q_S` sum by a
//! finite `M_S` sum, and it does not claim the Hilbert-space quotient map `E`.

use std::collections::BTreeMap;
use std::fmt;

use crate::semilocal_factorizable_poisson::LocalBallSpec;
use crate::semilocal_qs::{QsArithmeticError, QsRational, QsUnitMonoidDecomposition};
use crate::semilocal_trace_contract::FinitePlaceSet;
use crate::semilocal_unit_orbit::{SemilocalUnitOrbitTransport, UnitOrbitError};

/// One explicit non-zero diagonal sample together with its exact orbit data.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QsOrbitTerm {
    sample: QsRational,
    decomposition: QsUnitMonoidDecomposition,
    transported_balls: Vec<LocalBallSpec>,
}

impl QsOrbitTerm {
    #[inline]
    pub fn sample(&self) -> QsRational {
        self.sample
    }

    #[inline]
    pub fn decomposition(&self) -> &QsUnitMonoidDecomposition {
        &self.decomposition
    }

    #[inline]
    pub fn monoid_representative(&self) -> u64 {
        self.decomposition.monoid_element()
    }

    #[inline]
    pub fn unit_sign(&self) -> i8 {
        self.decomposition.unit_sign()
    }

    #[inline]
    pub fn unit_exponents(&self) -> &[(u64, i32)] {
        self.decomposition.unit_exponents()
    }

    #[inline]
    pub fn transported_balls(&self) -> &[LocalBallSpec] {
        &self.transported_balls
    }

    /// Recompose the original exact rational as an audit of the grouping data.
    #[inline]
    pub fn recomposed_sample(&self) -> QsRational {
        self.decomposition.recompose()
    }
}

/// Finite set of explicit samples sharing one unique `m in M_S` representative.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QsOrbitGroup {
    monoid_representative: u64,
    terms: Vec<QsOrbitTerm>,
}

impl QsOrbitGroup {
    #[inline]
    pub fn monoid_representative(&self) -> u64 {
        self.monoid_representative
    }

    #[inline]
    pub fn terms(&self) -> &[QsOrbitTerm] {
        &self.terms
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OrbitGroupingError {
    Qs(QsArithmeticError),
    UnitOrbit(UnitOrbitError),
}

impl fmt::Display for OrbitGroupingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Qs(error) => write!(f, "Q_S orbit decomposition failed: {error}"),
            Self::UnitOrbit(error) => write!(f, "Q_S unit transport failed: {error}"),
        }
    }
}

impl std::error::Error for OrbitGroupingError {}

impl From<QsArithmeticError> for OrbitGroupingError {
    fn from(value: QsArithmeticError) -> Self {
        Self::Qs(value)
    }
}

impl From<UnitOrbitError> for OrbitGroupingError {
    fn from(value: UnitOrbitError) -> Self {
        Self::UnitOrbit(value)
    }
}

/// Group explicit non-zero `Q_S` samples by their unique `M_S` representative.
///
/// The returned groups are sorted by increasing `m`. Terms within one group
/// preserve the input order. Every term retains the exact `Q_S^*` unit data and
/// the transported finite local balls; no local coordinate is projected away.
pub fn group_qs_samples_by_m(
    samples: &[QsRational],
    places: &FinitePlaceSet,
    local_balls: &[LocalBallSpec],
) -> Result<Vec<QsOrbitGroup>, OrbitGroupingError> {
    let mut grouped: BTreeMap<u64, Vec<QsOrbitTerm>> = BTreeMap::new();

    for &sample in samples {
        let decomposition = sample.unit_monoid_decomposition(places)?;
        let action = SemilocalUnitOrbitTransport::from_decomposition(&decomposition, places)?;
        let transported_balls = action.transport_complete_product(local_balls)?;
        let representative = decomposition.monoid_element();

        grouped.entry(representative).or_default().push(QsOrbitTerm {
            sample,
            decomposition,
            transported_balls,
        });
    }

    Ok(grouped
        .into_iter()
        .map(|(monoid_representative, terms)| QsOrbitGroup {
            monoid_representative,
            terms,
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multiple_units_group_under_the_same_m_representative() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let balls = [LocalBallSpec::new(2, 1), LocalBallSpec::new(3, -2)];
        let samples = [
            QsRational::new(5, 1, &places).unwrap(),
            QsRational::new(10, 1, &places).unwrap(),
            QsRational::new(5, 2, &places).unwrap(),
            QsRational::new(15, 1, &places).unwrap(),
            QsRational::new(45, 8, &places).unwrap(),
            QsRational::new(-5, 1, &places).unwrap(),
            QsRational::new(7, 1, &places).unwrap(),
            QsRational::new(14, 1, &places).unwrap(),
        ];

        let groups = group_qs_samples_by_m(&samples, &places, &balls).unwrap();
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].monoid_representative(), 5);
        assert_eq!(groups[0].len(), 6);
        assert_eq!(groups[1].monoid_representative(), 7);
        assert_eq!(groups[1].len(), 2);

        for group in &groups {
            for term in group.terms() {
                assert_eq!(term.monoid_representative(), group.monoid_representative());
                assert_eq!(term.recomposed_sample(), term.sample());
            }
        }
    }

    #[test]
    fn transported_local_coordinates_are_retained_per_orbit_term() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let balls = [LocalBallSpec::new(2, 1), LocalBallSpec::new(3, -2)];
        let samples = [
            QsRational::new(5, 1, &places).unwrap(),
            QsRational::new(45, 8, &places).unwrap(),
        ];

        let groups = group_qs_samples_by_m(&samples, &places, &balls).unwrap();
        let terms = groups[0].terms();

        assert_eq!(terms[0].unit_exponents(), &[(2, 0), (3, 0)]);
        assert_eq!(terms[0].transported_balls(), &balls);

        assert_eq!(terms[1].unit_exponents(), &[(2, -3), (3, 2)]);
        assert_eq!(
            terms[1].transported_balls(),
            &[LocalBallSpec::new(2, 4), LocalBallSpec::new(3, -4)]
        );
    }

    #[test]
    fn zero_is_rejected_instead_of_being_silently_assigned_an_orbit() {
        let places = FinitePlaceSet::new(vec![2]).unwrap();
        let zero = QsRational::new(0, 8, &places).unwrap();
        let error = group_qs_samples_by_m(&[zero], &places, &[LocalBallSpec::new(2, 0)])
            .unwrap_err();
        assert_eq!(
            error,
            OrbitGroupingError::Qs(QsArithmeticError::ZeroHasNoUnitMonoidDecomposition)
        );
    }
}
