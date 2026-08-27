use riemann_ndim_bench::semilocal::{
    ProlateParity, alpha_sequence, archimedean_a2, archimedean_a2_prime,
    build_kprime_closed, build_kprime_unsimplified, crossing_derivatives,
    merged_response_stats, merged_total_abs_derivative,
    sign_corrected_min_diagonal_dominance_margin,
};

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    let error = (actual - expected).abs();
    assert!(
        error <= tolerance,
        "actual={actual:.16e} expected={expected:.16e} error={error:.3e} tolerance={tolerance:.3e}"
    );
}

#[test]
fn alpha_recurrence_matches_first_exact_terms() {
    let alpha = alpha_sequence(7);
    let expected = [
        1.0,
        -0.5,
        3.0 / 8.0,
        -5.0 / 16.0,
        35.0 / 128.0,
        -63.0 / 256.0,
        231.0 / 1024.0,
    ];

    for (actual, expected) in alpha.into_iter().zip(expected) {
        assert_close(actual, expected, 1e-15);
    }
}

#[test]
fn source_a2_derivative_matches_unsimplified_first_order_coefficient() {
    let alpha = alpha_sequence(64);

    for n in 0..63 {
        let r_n = 2.0_f64.sqrt() * (alpha[n + 1] - alpha[n]);
        let unsimplified = 2.0 * archimedean_a2(n) * r_n;
        let closed = archimedean_a2_prime(n, alpha[n]);
        assert_close(closed, unsimplified, 2e-13);
    }
}

#[test]
fn closed_kprime_matches_source_derived_formula() {
    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        let closed = build_kprime_closed(256, parity);
        let source = build_kprime_unsimplified(256, parity);

        assert_eq!(closed.diagonal().len(), source.diagonal().len());
        assert_eq!(closed.off_diagonal().len(), source.off_diagonal().len());

        let max_diag_error = closed
            .diagonal()
            .iter()
            .zip(source.diagonal())
            .map(|(left, right)| (left - right).abs())
            .fold(0.0_f64, f64::max);
        let max_off_error = closed
            .off_diagonal()
            .iter()
            .zip(source.off_diagonal())
            .map(|(left, right)| (left - right).abs())
            .fold(0.0_f64, f64::max);

        assert!(
            max_diag_error < 5e-14,
            "{parity:?}: diagonal error {max_diag_error:e}"
        );
        assert!(
            max_off_error < 5e-14,
            "{parity:?}: off-diagonal error {max_off_error:e}"
        );
    }
}

#[test]
fn sign_corrected_kprime_is_strictly_diagonally_dominant_through_1024() {
    for block_size in [1, 2, 4, 16, 64, 256, 1024] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let margin = sign_corrected_min_diagonal_dominance_margin(block_size, parity);
            assert!(
                margin > 0.0,
                "block_size={block_size} parity={parity:?} margin={margin:e}"
            );
        }
    }
}

#[test]
fn every_tested_crossing_has_the_proved_parity_sign() {
    for block_size in [1, 4, 16, 64] {
        let plus = crossing_derivatives(block_size, ProlateParity::WPlus).unwrap();
        let minus = crossing_derivatives(block_size, ProlateParity::WMinus).unwrap();

        assert!(plus.iter().all(|crossing| crossing.lambda_prime < 0.0));
        assert!(minus.iter().all(|crossing| crossing.lambda_prime > 0.0));
    }
}

#[test]
fn reproduces_independent_large_prime_regression_targets() {
    let targets = [
        (16, 0.2826320800294, 0.2256334290480),
        (24, 0.2611374499942, 0.1874482544496),
        (32, 0.2457556928890, 0.1604809159874),
    ];

    for (m, expected_mean, expected_trimmed) in targets {
        let stats = merged_response_stats(m).unwrap();
        assert_close(m as f64 * stats.mean_abs, expected_mean, 5e-12);
        assert_close(
            m as f64 * stats.trimmed_mean_abs,
            expected_trimmed,
            5e-12,
        );
    }
}

#[test]
fn reproduces_documented_m128_shape_statistics() {
    // m=128 exercises the full merged EVD/statistics path while keeping the
    // default test suite below the expensive m=1024 dense eigensolve.
    let m = 128;
    let stats = merged_response_stats(m).unwrap();
    assert_close(m as f64 * stats.mean_abs, 0.1754882381797, 5e-11);
    assert_close(m as f64 * stats.trimmed_mean_abs, 0.07059034265272, 5e-11);
    assert_close(m as f64 * stats.rms, 0.4130934070885, 5e-11);
    assert_close((m as f64).sqrt() * stats.linf, 0.2471011380024, 5e-11);
}

#[test]
#[ignore = "expensive dense faer EVD; run before promoting high-m asymptotics"]
fn reproduces_high_block_total_response_and_signs_through_1024() {
    let targets = [
        (128, 3.970845543531),
        (256, 4.640481894221),
        (512, 5.359223651882),
        (1024, 6.126883687871),
    ];

    for (m, expected) in targets {
        let actual = merged_total_abs_derivative(m).unwrap();
        assert_close(actual, expected, 5e-9);
    }

    let plus = crossing_derivatives(1024, ProlateParity::WPlus).unwrap();
    let minus = crossing_derivatives(1024, ProlateParity::WMinus).unwrap();
    assert!(plus.iter().all(|crossing| crossing.lambda_prime < 0.0));
    assert!(minus.iter().all(|crossing| crossing.lambda_prime > 0.0));
}
