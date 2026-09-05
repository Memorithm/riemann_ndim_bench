use std::f64::consts::PI;

use num_complex::Complex64;
use riemann_ndim_bench::semilocal_fourier_multiplier::{
    convergent_dirichlet_prefix, critical_line_dirichlet_exponent,
    finite_euler_deletion_factor, semilocal_multiplier_from_zeta,
};
use riemann_ndim_bench::semilocal_poisson::SemilocalPoissonMonoid;
use riemann_ndim_bench::semilocal_trace_contract::FinitePlaceSet;

fn assert_certified_target(prefix_value: f64, tail: f64, target: f64, label: &str) {
    let error = (prefix_value - target).abs();
    assert!(
        error <= tail + 8.0 * f64::EPSILON * target.abs().max(1.0),
        "{label}: value={prefix_value:.16e} target={target:.16e} error={error:.3e} tail={tail:.3e}"
    );
}

#[test]
fn convergent_semilocal_dirichlet_series_matches_known_euler_deletions() {
    let cases = [
        (vec![], PI * PI / 6.0, "S={infinity}"),
        (vec![2], PI * PI / 8.0, "S={infinity,2}"),
        (vec![2, 3], PI * PI / 9.0, "S={infinity,2,3}"),
    ];

    for (finite_primes, target, label) in cases {
        let monoid =
            SemilocalPoissonMonoid::new(FinitePlaceSet::new(finite_primes).unwrap());
        let prefix = convergent_dirichlet_prefix(&monoid, 2.0, 50_000).unwrap();
        assert_certified_target(prefix.value(), prefix.tail_upper_bound(), target, label);
        assert!(prefix.term_count() >= 1);
        assert_eq!(prefix.max_m(), 50_000);
    }
}

#[test]
fn finite_euler_factor_multiplies_zeta_two_to_the_same_targets() {
    let zeta_two = Complex64::new(PI * PI / 6.0, 0.0);
    let exponent = Complex64::new(2.0, 0.0);

    for (finite_primes, target) in [
        (vec![], PI * PI / 6.0),
        (vec![2], PI * PI / 8.0),
        (vec![2, 3], PI * PI / 9.0),
    ] {
        let places = FinitePlaceSet::new(finite_primes).unwrap();
        let multiplier =
            semilocal_multiplier_from_zeta(zeta_two, &places, exponent).unwrap();
        assert!((multiplier.re - target).abs() < 4.0e-15);
        assert!(multiplier.im.abs() < 4.0e-15);
    }
}

#[test]
fn critical_line_path_exposes_only_the_finite_source_factor() {
    let places = FinitePlaceSet::new(vec![2, 3, 5]).unwrap();
    let exponent = critical_line_dirichlet_exponent(14.134_725).unwrap();
    let factor = finite_euler_deletion_factor(&places, exponent).unwrap();

    assert!(factor.re.is_finite());
    assert!(factor.im.is_finite());
    assert_eq!(exponent.re, 0.5);
    assert_eq!(exponent.im, -14.134_725);
}
