use riemann_ndim_bench::semilocal_factorizable_poisson::LocalBallSpec;
use riemann_ndim_bench::semilocal_qs::QsRational;
use riemann_ndim_bench::semilocal_trace_contract::FinitePlaceSet;
use riemann_ndim_bench::semilocal_unit_orbit::{SemilocalUnitOrbitTransport, UnitOrbitError};

#[test]
fn qs_unit_action_preserves_effective_local_exponents() {
    let places = FinitePlaceSet::new(vec![2, 3, 5]).unwrap();
    let q = QsRational::new(-675, 32, &places).unwrap();
    let decomposition = q.unit_monoid_decomposition(&places).unwrap();
    let action = SemilocalUnitOrbitTransport::from_decomposition(&decomposition, &places).unwrap();

    // -675/32 = -(2^-5)(3^3)(5^2) * 1.
    assert_eq!(action.unit_sign(), -1);
    assert_eq!(action.unit_exponents(), &[(2, -5), (3, 3), (5, 2)]);
    assert_eq!(decomposition.monoid_element(), 1);

    let original = [
        LocalBallSpec::new(2, -2),
        LocalBallSpec::new(3, 4),
        LocalBallSpec::new(5, 0),
    ];
    let transported = action.transport_complete_product(&original).unwrap();
    assert_eq!(
        transported,
        vec![
            LocalBallSpec::new(2, 3),
            LocalBallSpec::new(3, 1),
            LocalBallSpec::new(5, -2),
        ]
    );

    for ball in original {
        assert_eq!(action.compensated_exponent(ball).unwrap(), ball.exponent());
    }
}

#[test]
fn unit_transport_rejects_incomplete_or_mismatched_local_products() {
    let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
    let q = QsRational::new(45, 8, &places).unwrap();
    let decomposition = q.unit_monoid_decomposition(&places).unwrap();
    let action = SemilocalUnitOrbitTransport::from_decomposition(&decomposition, &places).unwrap();

    assert_eq!(
        action
            .transport_complete_product(&[LocalBallSpec::new(2, 0)])
            .unwrap_err(),
        UnitOrbitError::IncompleteLocalProduct
    );
    assert_eq!(
        action.transport_ball(LocalBallSpec::new(5, 0)).unwrap_err(),
        UnitOrbitError::PrimeOutsideUnitPlaces { prime: 5 }
    );
}

#[test]
fn declared_place_set_must_match_the_decomposition_origin() {
    let source_places = FinitePlaceSet::new(vec![2, 3]).unwrap();
    let q = QsRational::new(45, 8, &source_places).unwrap();
    let decomposition = q.unit_monoid_decomposition(&source_places).unwrap();
    let wrong_places = FinitePlaceSet::new(vec![2, 5]).unwrap();

    assert_eq!(
        SemilocalUnitOrbitTransport::from_decomposition(&decomposition, &wrong_places).unwrap_err(),
        UnitOrbitError::PlaceSetMismatch
    );
}
