use riemann_ndim_bench::semilocal::ProlateParity;
use riemann_ndim_bench::semilocal_resolvent::SignCorrectedResolventTraceKernel;
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
fn row_contributions_sum_to_resolvent_trace() {
    for block_size in [0_usize, 1, 2, 8, 64, 256] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let kernel = SignCorrectedResolventTraceKernel::new(block_size, parity);
            for shift in [0.0_f64, 1.0e-8, 1.0e-4, 1.0, 100.0] {
                let rows = kernel.row_contributions(shift).unwrap();
                let row_sum: f64 = rows.iter().sum();
                let trace = kernel.trace(shift).unwrap();

                assert_eq!(rows.len(), block_size);
                assert!(rows.iter().all(|value| value.is_finite()));
                assert_close_scaled(
                    row_sum,
                    trace,
                    5.0e-14,
                    5.0e-14,
                    &format!("m={block_size} parity={parity:?} shift={shift:e} row sum"),
                );
            }
        }
    }
}

#[test]
fn row_resolved_trace_matches_spectral_oracle() {
    for block_size in [8_usize, 32, 128] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let spectral = crossing_derivatives_tridiagonal(block_size, parity).unwrap();
            let kernel = SignCorrectedResolventTraceKernel::new(block_size, parity);
            let sign = parity.sign_correction();

            for shift in [0.0_f64, 1.0e-6, 1.0e-2, 1.0] {
                let expected: f64 = spectral
                    .iter()
                    .map(|crossing| sign * crossing.mu_prime / (crossing.mu + shift))
                    .sum();
                let actual: f64 = kernel.row_contributions(shift).unwrap().into_iter().sum();

                assert_close_scaled(
                    actual,
                    expected,
                    5.0e-11,
                    1.0e-11,
                    &format!("m={block_size} parity={parity:?} shift={shift:e} spectral row sum"),
                );
            }
        }
    }
}

#[test]
fn row_resolution_preserves_finite_section_boundary_rows() {
    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        let kernel = SignCorrectedResolventTraceKernel::new(64, parity);
        let rows = kernel.row_contributions(1.0e-3).unwrap();

        assert_eq!(rows.len(), 64);
        assert!(rows[0].is_finite());
        assert!(rows[63].is_finite());
        assert!(rows.iter().any(|value| value.abs() > 0.0));
    }
}
