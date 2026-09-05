use riemann_ndim_bench::archimedean_poisson::{
    compare_source_poisson_identity, source_poisson_fixture_at_zero,
    source_poisson_fixture_fourier_at_zero,
};

#[test]
fn manufactured_fixture_meets_the_two_source_zero_conditions() {
    assert!(source_poisson_fixture_at_zero().abs() < 2.0e-16);
    assert!(source_poisson_fixture_fourier_at_zero().abs() < 2.0e-16);
}

#[test]
fn source_poisson_identity_is_certified_across_reciprocal_scales() {
    for x in [0.5_f64, 0.75, 1.0, 1.5, 2.0] {
        let comparison = compare_source_poisson_identity(x, 32).unwrap();
        let truncation = comparison.combined_tail_bound();
        let roundoff = 64.0 * f64::EPSILON
            * comparison
                .left()
                .value()
                .abs()
                .max(comparison.right().value().abs())
                .max(1.0);
        let tolerance = truncation + roundoff;

        assert!(
            comparison.residual().abs() <= tolerance,
            "x={x}: left={:.16e} right={:.16e} residual={:.3e} tail={truncation:.3e} tolerance={tolerance:.3e}",
            comparison.left().value(),
            comparison.right().value(),
            comparison.residual().abs(),
        );
        assert_eq!(comparison.x(), x);
        assert_eq!(comparison.left().max_n(), 32);
        assert_eq!(comparison.right().max_n(), 32);
    }
}

#[test]
fn certified_gaussian_tails_are_negligible_on_the_test_window() {
    for x in [0.5_f64, 0.75, 1.0, 1.5, 2.0] {
        let comparison = compare_source_poisson_identity(x, 32).unwrap();
        assert!(
            comparison.combined_tail_bound() < 1.0e-35,
            "x={x}: tail={:.3e}",
            comparison.combined_tail_bound()
        );
    }
}
