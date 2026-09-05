use riemann_ndim_bench::semilocal_poisson::SemilocalPoissonMonoid;
use riemann_ndim_bench::semilocal_qs::{QsArithmeticError, QsRational};
use riemann_ndim_bench::semilocal_trace_contract::FinitePlaceSet;

#[test]
fn source_decomposition_lands_in_the_existing_m_s_monoid() {
    let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
    let monoid = SemilocalPoissonMonoid::new(places.clone());

    for (numerator, denominator, expected_m) in [
        (45_i64, 8_u64, 5_u64),
        (-150, 72, 25),
        (7, 12, 7),
        (1000, 9, 125),
    ] {
        let rational = QsRational::new(numerator, denominator, &places).unwrap();
        let decomposition = rational.unit_monoid_decomposition(&places).unwrap();
        assert_eq!(decomposition.monoid_element(), expected_m);
        assert!(monoid.contains(decomposition.monoid_element()));
        assert_eq!(decomposition.recompose(), rational);
    }
}

#[test]
fn equivalent_rational_presentations_have_the_same_canonical_split() {
    let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
    let first = QsRational::new(45, 8, &places).unwrap();
    let second = QsRational::new(90, 16, &places).unwrap();

    assert_eq!(first, second);
    assert_eq!(
        first.unit_monoid_decomposition(&places).unwrap(),
        second.unit_monoid_decomposition(&places).unwrap()
    );
}

#[test]
fn archimedean_only_q_s_reduces_to_the_integers() {
    let places = FinitePlaceSet::new(vec![]).unwrap();
    let integer = QsRational::new(-35, 1, &places).unwrap();
    let decomposition = integer.unit_monoid_decomposition(&places).unwrap();

    assert_eq!(decomposition.unit_sign(), -1);
    assert!(decomposition.unit_exponents().is_empty());
    assert_eq!(decomposition.monoid_element(), 35);
    assert_eq!(decomposition.recompose(), integer);

    assert_eq!(
        QsRational::new(1, 2, &places).unwrap_err(),
        QsArithmeticError::DenominatorOutsidePlaceSet { residual: 2 }
    );
}

#[test]
fn finite_place_exponents_are_the_unique_unit_signature() {
    let places = FinitePlaceSet::new(vec![2, 3, 5]).unwrap();
    let rational = QsRational::new(-1350, 32, &places).unwrap();
    let decomposition = rational.unit_monoid_decomposition(&places).unwrap();

    assert_eq!(decomposition.unit_sign(), -1);
    assert_eq!(decomposition.exponent_for(2), Some(-4));
    assert_eq!(decomposition.exponent_for(3), Some(3));
    assert_eq!(decomposition.exponent_for(5), Some(2));
    assert_eq!(decomposition.monoid_element(), 1);
    assert_eq!(decomposition.recompose(), rational);
}
