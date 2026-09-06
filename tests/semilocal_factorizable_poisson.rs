use riemann_ndim_bench::semilocal_factorizable_poisson::{
    LocalBallSpec, compare_factorizable_ball_poisson,
};

fn assert_certified(specs: &[LocalBallSpec]) {
    let comparison = compare_factorizable_ball_poisson(specs, 256).unwrap();
    let scale = comparison
        .left_value()
        .abs()
        .max(comparison.right_value().abs())
        .max(1.0);
    let roundoff = 5.0e-13 * scale;
    assert!(
        comparison.residual().abs() <= comparison.combined_tail_bound() + roundoff,
        "specs={specs:?} residual={:.3e} tail={:.3e} roundoff={:.3e}",
        comparison.residual(),
        comparison.combined_tail_bound(),
        roundoff
    );
}

#[test]
fn certified_semilocal_product_poisson_holds_for_nontrivial_local_balls() {
    assert_certified(&[LocalBallSpec::new(2, 1)]);
    assert_certified(&[LocalBallSpec::new(2, -2)]);
    assert_certified(&[LocalBallSpec::new(3, 2)]);
    assert_certified(&[LocalBallSpec::new(5, -1)]);
}

#[test]
fn certified_semilocal_product_poisson_composes_across_places() {
    assert_certified(&[
        LocalBallSpec::new(2, 1),
        LocalBallSpec::new(3, -1),
        LocalBallSpec::new(5, 1),
    ]);
}
