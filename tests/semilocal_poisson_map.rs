use riemann_ndim_bench::semilocal_poisson::SemilocalPoissonMonoid;
use riemann_ndim_bench::semilocal_trace_contract::{FinitePlaceSet, SemilocalSpaceContract};

#[test]
fn source_space_contract_drives_the_poisson_monoid() {
    let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
    let space = SemilocalSpaceContract::qs_self_dual(places);
    let monoid = SemilocalPoissonMonoid::from_space(&space);

    assert_eq!(monoid.elements_through(15), vec![1, 5, 7, 11, 13]);
    assert_eq!(monoid.places().finite_primes(), &[2, 3]);
}

#[test]
fn compactly_truncated_e_sum_is_finite_and_auditable() {
    let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
    let monoid = SemilocalPoissonMonoid::new(places);
    let result = monoid
        .finite_e_sum(9.0, 15, |m| m as f64)
        .expect("positive semilocal modulus");

    // M_S intersect [1,15] = {1,5,7,11,13}; their sum is 37.
    assert_eq!(result.term_count(), 5);
    assert_eq!(result.max_m(), 15);
    assert_eq!(result.raw_sum(), 37.0);
    assert_eq!(result.value(), 111.0);
}

#[test]
fn adding_a_finite_place_removes_exactly_its_multiples_from_the_prefix() {
    let arch = SemilocalPoissonMonoid::new(FinitePlaceSet::new(vec![]).unwrap());
    let with_two = SemilocalPoissonMonoid::new(FinitePlaceSet::new(vec![2]).unwrap());

    assert_eq!(arch.elements_through(8), vec![1, 2, 3, 4, 5, 6, 7, 8]);
    assert_eq!(with_two.elements_through(8), vec![1, 3, 5, 7]);
}
