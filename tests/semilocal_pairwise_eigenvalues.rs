use riemann_ndim_bench::semilocal::{ProlateParity, merged_total_abs_derivative};
use riemann_ndim_bench::semilocal_tridiagonal::{
    PairwiseError, crossing_derivatives_tridiagonal,
    first_order_sqrt_spectrum_tridiagonal_eigenvalues_only, pairwise_first_order_derivatives,
    pairwise_first_order_trace_derivative, pairwise_total_variation_derivative,
    pairwise_total_variation_extrapolation, quadratic_even_power_extrapolate,
};

fn assert_close_scaled(actual: f64, expected: f64, relative: f64, absolute: f64, label: &str) {
    let error = (actual - expected).abs();
    let scale = actual.abs().max(expected.abs());
    let tolerance = absolute.max(relative * scale);
    assert!(
        error <= tolerance,
        "{label}: actual={actual:.16e} expected={expected:.16e} error={error:.3e} tolerance={tolerance:.3e}"
    );
}

#[test]
fn eigenvalues_only_zero_step_matches_full_tridiagonal_spectrum() {
    for block_size in [1_usize, 2, 4, 8, 16, 32, 64, 128] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let exact = crossing_derivatives_tridiagonal(block_size, parity).unwrap();
            let eigenvalues_only =
                first_order_sqrt_spectrum_tridiagonal_eigenvalues_only(block_size, parity, 0.0)
                    .unwrap();

            assert_eq!(exact.len(), eigenvalues_only.len());
            for (index, (reference, candidate)) in exact.iter().zip(&eigenvalues_only).enumerate() {
                assert_close_scaled(
                    *candidate,
                    reference.lambda,
                    2e-13,
                    2e-15,
                    &format!("m={block_size} parity={parity:?} index={index} lambda"),
                );
            }
        }
    }
}

#[test]
fn pairwise_quadratic_extrapolation_matches_exact_first_order_total_variation() {
    for block_size in [8_usize, 16, 32, 64, 128] {
        let exact = merged_total_abs_derivative(block_size).unwrap();
        let diagnostic = pairwise_total_variation_extrapolation(block_size, 5.0e-4).unwrap();

        assert_close_scaled(
            diagnostic.quadratic_h2_h4,
            exact,
            2e-8,
            2e-11,
            &format!("m={block_size} extrapolated total variation"),
        );
    }
}

#[test]
fn pairwise_trace_preserves_parity_signs() {
    for block_size in [8_usize, 32, 128] {
        let plus =
            pairwise_first_order_trace_derivative(block_size, ProlateParity::WPlus, 2.5e-4)
                .unwrap();
        let minus =
            pairwise_first_order_trace_derivative(block_size, ProlateParity::WMinus, 2.5e-4)
                .unwrap();

        assert!(plus < 0.0, "m={block_size}: W+ trace derivative={plus:e}");
        assert!(minus > 0.0, "m={block_size}: W- trace derivative={minus:e}");
    }
}

#[test]
fn empty_pairwise_block_has_zero_variation() {
    assert!(
        pairwise_first_order_derivatives(0, ProlateParity::WPlus, 5.0e-4)
            .unwrap()
            .is_empty()
    );
    assert_eq!(pairwise_total_variation_derivative(0, 5.0e-4).unwrap(), 0.0);
}

#[test]
fn invalid_pairwise_steps_are_rejected() {
    for step in [0.0, -1.0e-4, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            pairwise_total_variation_derivative(8, step),
            Err(PairwiseError::InvalidStep { .. })
        ));
    }
}

#[test]
fn quadratic_even_power_extrapolation_eliminates_h2_and_h4_terms() {
    let exact = 3.25_f64;
    let a2 = -1.75_f64;
    let a4 = 0.625_f64;
    let h = 0.4_f64;

    let model = |step: f64| exact + a2 * step.powi(2) + a4 * step.powi(4);
    let extrapolated = quadratic_even_power_extrapolate(model(h), model(h / 2.0), model(h / 4.0));

    assert_close_scaled(extrapolated, exact, 0.0, 2e-14, "synthetic even-power extrapolation");
}

#[test]
#[ignore = "expensive high-dimensional pairwise checkpoint"]
fn reproduces_documented_m8192_pairwise_checkpoint() {
    let diagnostic = pairwise_total_variation_extrapolation(8192, 5.0e-4).unwrap();

    assert_close_scaled(diagnostic.d_h, 8.722779121464640, 0.0, 5e-9, "m=8192 D(h)");
    assert_close_scaled(
        diagnostic.d_h2,
        8.722587591872799,
        0.0,
        5e-9,
        "m=8192 D(h/2)",
    );
    assert_close_scaled(
        diagnostic.d_h4,
        8.722540424778753,
        0.0,
        5e-9,
        "m=8192 D(h/4)",
    );
    assert_close_scaled(
        diagnostic.quadratic_h2_h4,
        8.722524765996642,
        0.0,
        5e-9,
        "m=8192 Q123",
    );
}

#[test]
#[ignore = "expensive high-dimensional pairwise checkpoint"]
fn reproduces_documented_m16384_pairwise_checkpoint() {
    let diagnostic = pairwise_total_variation_extrapolation(16384, 5.0e-4).unwrap();

    assert_close_scaled(diagnostic.d_h, 9.685679471063766, 0.0, 5e-9, "m=16384 D(h)");
    assert_close_scaled(
        diagnostic.d_h2,
        9.685297949399001,
        0.0,
        5e-9,
        "m=16384 D(h/2)",
    );
    assert_close_scaled(
        diagnostic.d_h4,
        9.685202556488507,
        0.0,
        5e-9,
        "m=16384 D(h/4)",
    );
    assert_close_scaled(
        diagnostic.quadratic_h2_h4,
        9.685170757741071,
        0.0,
        5e-9,
        "m=16384 Q123",
    );
}
