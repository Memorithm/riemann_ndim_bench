use riemann_ndim_bench::weil_boundary::MultiplicativeSupport;
use riemann_ndim_bench::weil_support::{
    WeilSupportWindow, convolution_square_support_envelope, prime_is_excluded_by_support,
    source_place_set_below_integer_bound,
};

#[test]
fn compact_h1_support_inside_q_window_stays_inside_convolution_window() {
    let window = WeilSupportWindow::new(10.0).unwrap();
    let h1 = MultiplicativeSupport::new(0.4, 2.4).unwrap();
    assert!(window.contains_h1_support(h1));

    let convolution = convolution_square_support_envelope(h1).unwrap();
    assert!(window.contains_convolution_support(convolution));

    for prime in [11_u64, 13, 17, 101] {
        assert!(
            prime_is_excluded_by_support(prime, convolution),
            "prime {prime} should be excluded by support {convolution:?}"
        );
    }
}

#[test]
fn source_place_set_is_exactly_the_primes_strictly_below_q_for_integer_q() {
    let q10 = source_place_set_below_integer_bound(10).unwrap();
    assert_eq!(q10.finite_primes(), &[2, 3, 5, 7]);

    let q6 = source_place_set_below_integer_bound(6).unwrap();
    assert_eq!(q6.finite_primes(), &[2, 3, 5]);
}

#[test]
fn support_fact_and_conjectural_sufficiency_are_not_encoded_as_the_same_claim() {
    let window = WeilSupportWindow::new(10.0).unwrap();
    assert!(window.source_set_contains_prime(7));
    assert!(!window.source_set_contains_prime(11));

    // These are only membership facts for the source set S(q). There is
    // intentionally no API that returns a "Weil positivity proved" boolean.
}
