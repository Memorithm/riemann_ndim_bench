use riemann_ndim_bench::semilocal_padic_fourier::{PadicBall, PadicFourierError};
use riemann_ndim_bench::semilocal_qs::QsRational;
use riemann_ndim_bench::semilocal_trace_contract::FinitePlaceSet;

#[test]
fn undeclared_finite_place_is_rejected() {
    let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
    assert_eq!(
        PadicBall::new(5, 0, &places).unwrap_err(),
        PadicFourierError::PrimeOutsidePlaceSet { prime: 5 }
    );
}

#[test]
fn z_p_is_symbolically_self_dual_under_the_source_normalization() {
    let places = FinitePlaceSet::new(vec![2, 3]).unwrap();

    for prime in [2_u64, 3] {
        let unit_ball = PadicBall::new(prime, 0, &places).unwrap();
        let transformed = unit_ball.fourier_transform();

        assert_eq!(transformed.ball(), unit_ball);
        assert!(transformed.scale().is_one());
        assert_eq!(transformed.scale().prime(), prime);
    }
}

#[test]
fn general_ball_fourier_image_has_exact_dual_exponent_and_scale() {
    let places = FinitePlaceSet::new(vec![2]).unwrap();

    for exponent in [-5_i32, -2, -1, 0, 1, 2, 5] {
        let ball = PadicBall::new(2, exponent, &places).unwrap();
        let transformed = ball.fourier_transform();

        assert_eq!(transformed.ball().prime(), 2);
        assert_eq!(transformed.ball().exponent(), -i64::from(exponent));
        assert_eq!(transformed.scale().prime(), 2);
        assert_eq!(transformed.scale().exponent(), -i64::from(exponent));
    }
}

#[test]
fn applying_local_fourier_twice_restores_ball_and_unit_scale_exactly() {
    let places = FinitePlaceSet::new(vec![2, 5]).unwrap();

    for prime in [2_u64, 5] {
        for exponent in [-8_i32, -3, -1, 0, 1, 4, 9] {
            let ball = PadicBall::new(prime, exponent, &places).unwrap();
            let twice = ball.fourier_transform().fourier_transform().unwrap();

            assert_eq!(twice.ball(), ball);
            assert!(twice.scale().is_one());
            assert_eq!(twice.scale().prime(), prime);
        }
    }
}

#[test]
fn diagonal_membership_matches_exact_p_adic_valuation_thresholds() {
    let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
    let q = QsRational::new(45, 8, &places).unwrap(); // v_2=-3, v_3=2
    let zero = QsRational::new(0, 1, &places).unwrap();

    assert!(PadicBall::new(2, -3, &places).unwrap().contains_diagonal(q));
    assert!(!PadicBall::new(2, -2, &places).unwrap().contains_diagonal(q));
    assert!(PadicBall::new(3, 2, &places).unwrap().contains_diagonal(q));
    assert!(!PadicBall::new(3, 3, &places).unwrap().contains_diagonal(q));

    for exponent in [-10_i32, 0, 10] {
        assert!(
            PadicBall::new(2, exponent, &places)
                .unwrap()
                .contains_diagonal(zero)
        );
    }
}
