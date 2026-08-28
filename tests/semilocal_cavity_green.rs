use faer::linalg::solvers::SelfAdjointEigen;
use faer::Side;

use riemann_ndim_bench::semilocal::{
    ProlateParity, build_k0, build_kprime_closed,
};
use riemann_ndim_bench::semilocal_resolvent::SignCorrectedResolventTraceKernel;

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
fn cavity_green_bands_match_dense_spectral_oracle() {
    for block_size in [1_usize, 2, 4, 16, 64] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let k0 = build_k0(block_size, parity);
            let dense = k0.to_dense();
            let evd = SelfAdjointEigen::new(dense.as_ref(), Side::Lower).unwrap();
            let eigenvectors = evd.U();
            let eigenvalues = evd.S();
            let kernel = SignCorrectedResolventTraceKernel::new(block_size, parity);

            for shift in [0.0_f64, 1.0e-8, 1.0e-4, 1.0] {
                let cavity = kernel.cavity_green_bands(shift).unwrap();

                assert_eq!(cavity.diagonal().len(), block_size);
                assert_eq!(
                    cavity.off_diagonal().len(),
                    block_size.saturating_sub(1)
                );
                assert_eq!(cavity.left_denominators().len(), block_size);
                assert_eq!(cavity.right_denominators().len(), block_size);
                assert!(cavity.left_denominators().iter().all(|value| *value > 0.0));
                assert!(cavity.right_denominators().iter().all(|value| *value > 0.0));

                for i in 0..block_size {
                    let expected: f64 = (0..block_size)
                        .map(|j| {
                            let u = eigenvectors[(i, j)];
                            u * u / (eigenvalues[j] + shift)
                        })
                        .sum();
                    assert_close_scaled(
                        cavity.diagonal()[i],
                        expected,
                        2.0e-10,
                        2.0e-11,
                        &format!(
                            "m={block_size} parity={parity:?} shift={shift:e} diag[{i}]"
                        ),
                    );
                }

                for i in 0..block_size.saturating_sub(1) {
                    let expected: f64 = (0..block_size)
                        .map(|j| {
                            eigenvectors[(i, j)] * eigenvectors[(i + 1, j)]
                                / (eigenvalues[j] + shift)
                        })
                        .sum();
                    assert_close_scaled(
                        cavity.off_diagonal()[i],
                        expected,
                        2.0e-10,
                        2.0e-11,
                        &format!(
                            "m={block_size} parity={parity:?} shift={shift:e} off[{i}]"
                        ),
                    );
                }
            }
        }
    }
}

#[test]
fn cavity_green_reconstructs_weighted_resolvent_trace() {
    for block_size in [1_usize, 2, 8, 64, 256] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let kernel = SignCorrectedResolventTraceKernel::new(block_size, parity);
            let h = build_kprime_closed(block_size, parity);
            let sign = parity.sign_correction();

            for shift in [0.0_f64, 1.0e-6, 1.0e-2, 1.0] {
                let cavity = kernel.cavity_green_bands(shift).unwrap();
                let mut trace = 0.0;

                for i in 0..block_size {
                    trace += sign * h.diagonal()[i] * cavity.diagonal()[i];
                }
                for i in 0..block_size.saturating_sub(1) {
                    trace += 2.0 * sign * h.off_diagonal()[i] * cavity.off_diagonal()[i];
                }

                assert_close_scaled(
                    trace,
                    kernel.trace(shift).unwrap(),
                    2.0e-12,
                    2.0e-12,
                    &format!("m={block_size} parity={parity:?} shift={shift:e} cavity trace"),
                );
            }
        }
    }
}

#[test]
fn empty_cavity_green_has_empty_bands() {
    let kernel = SignCorrectedResolventTraceKernel::new(0, ProlateParity::WPlus);
    let cavity = kernel.cavity_green_bands(0.0).unwrap();

    assert!(cavity.left_denominators().is_empty());
    assert!(cavity.right_denominators().is_empty());
    assert!(cavity.diagonal().is_empty());
    assert!(cavity.off_diagonal().is_empty());
}
