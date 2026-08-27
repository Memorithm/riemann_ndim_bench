use riemann_ndim_bench::semilocal::{ProlateParity, crossing_derivatives};
use riemann_ndim_bench::semilocal_tridiagonal::crossing_derivatives_tridiagonal;

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
fn empty_tridiagonal_block_returns_no_crossings() {
    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        assert!(crossing_derivatives_tridiagonal(0, parity).unwrap().is_empty());
    }
}

#[test]
fn tridiagonal_evd_matches_dense_reference_through_128() {
    for block_size in [1_usize, 2, 4, 8, 16, 32, 64, 128] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let dense = crossing_derivatives(block_size, parity).unwrap();
            let tridiagonal = crossing_derivatives_tridiagonal(block_size, parity).unwrap();

            assert_eq!(dense.len(), tridiagonal.len());

            for (reference, optimized) in dense.iter().zip(&tridiagonal) {
                assert_eq!(reference.parity, optimized.parity);
                assert_eq!(reference.parity_index, optimized.parity_index);

                let prefix = format!(
                    "m={block_size} parity={parity:?} index={}",
                    reference.parity_index
                );

                assert_close_scaled(
                    optimized.mu,
                    reference.mu,
                    2e-13,
                    2e-15,
                    &format!("{prefix} mu"),
                );
                assert_close_scaled(
                    optimized.lambda,
                    reference.lambda,
                    2e-13,
                    2e-15,
                    &format!("{prefix} lambda"),
                );
                assert_close_scaled(
                    optimized.mu_prime,
                    reference.mu_prime,
                    2e-11,
                    2e-14,
                    &format!("{prefix} mu_prime"),
                );
                assert_close_scaled(
                    optimized.lambda_prime,
                    reference.lambda_prime,
                    2e-11,
                    2e-14,
                    &format!("{prefix} lambda_prime"),
                );
            }
        }
    }
}

#[test]
fn tridiagonal_path_preserves_proved_parity_signs() {
    for block_size in [1_usize, 4, 16, 64, 128] {
        let plus = crossing_derivatives_tridiagonal(block_size, ProlateParity::WPlus).unwrap();
        let minus = crossing_derivatives_tridiagonal(block_size, ProlateParity::WMinus).unwrap();

        assert!(plus.iter().all(|crossing| crossing.lambda_prime < 0.0));
        assert!(minus.iter().all(|crossing| crossing.lambda_prime > 0.0));
    }
}
