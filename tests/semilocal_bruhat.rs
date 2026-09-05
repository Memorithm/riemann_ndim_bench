use riemann_ndim_bench::semilocal_bruhat::ElementaryFiniteBruhatFactor;
use riemann_ndim_bench::semilocal_poisson::SemilocalPoissonMonoid;
use riemann_ndim_bench::semilocal_qs::QsRational;
use riemann_ndim_bench::semilocal_trace_contract::FinitePlaceSet;

#[test]
fn monoid_scaling_is_invisible_to_declared_finite_unit_ball_factors() {
    let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
    let monoid = SemilocalPoissonMonoid::new(places.clone());
    let factor = ElementaryFiniteBruhatFactor::new(places.clone());

    let base = QsRational::new(3, 4, &places).unwrap();
    for m in monoid.elements_through(31) {
        let numerator = 3_i64 * i64::try_from(m).unwrap();
        let scaled = QsRational::new(numerator, 4, &places).unwrap();

        assert_eq!(
            factor.diagonal_valuation(base, 2),
            factor.diagonal_valuation(scaled, 2)
        );
        assert_eq!(
            factor.diagonal_valuation(base, 3),
            factor.diagonal_valuation(scaled, 3)
        );
        assert_eq!(
            factor.evaluate_diagonal(base),
            factor.evaluate_diagonal(scaled)
        );
    }
}

#[test]
fn finite_factor_is_nontrivial_on_qs_diagonal_samples() {
    let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
    let factor = ElementaryFiniteBruhatFactor::new(places.clone());

    let samples = [
        (QsRational::new(1, 1, &places).unwrap(), 1_u8),
        (QsRational::new(6, 1, &places).unwrap(), 1_u8),
        (QsRational::new(1, 2, &places).unwrap(), 0_u8),
        (QsRational::new(1, 3, &places).unwrap(), 0_u8),
        (QsRational::new(5, 6, &places).unwrap(), 0_u8),
    ];

    for (q, expected) in samples {
        assert_eq!(factor.evaluate_diagonal(q), expected);
    }
}

#[test]
fn factorizable_diagonal_value_keeps_archimedean_and_finite_roles_separate() {
    let places = FinitePlaceSet::new(vec![2]).unwrap();
    let factor = ElementaryFiniteBruhatFactor::new(places.clone());
    let archimedean_value = -2.75;

    let integral = QsRational::new(3, 1, &places).unwrap();
    let non_integral = QsRational::new(3, 2, &places).unwrap();

    assert_eq!(
        factor.evaluate_factorizable_diagonal(integral, archimedean_value),
        archimedean_value
    );
    assert_eq!(
        factor.evaluate_factorizable_diagonal(non_integral, archimedean_value),
        0.0
    );
}
