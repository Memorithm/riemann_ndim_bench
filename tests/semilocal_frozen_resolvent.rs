use std::f64::consts::PI;

use riemann_ndim_bench::semilocal::{ProlateParity, build_k0, build_kprime_closed};
use riemann_ndim_bench::semilocal_resolvent::{
    ResolventTraceError, SignCorrectedResolventTraceKernel,
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

fn frozen_symbol_quadrature(row: usize, parity: ProlateParity, shift: f64) -> f64 {
    let block_size = row + 2;
    let k0 = build_k0(block_size, parity);
    let kprime = build_kprime_closed(block_size, parity);
    let sign = parity.sign_correction();

    let b = 0.5 * (k0.off_diagonal()[row - 1] + k0.off_diagonal()[row]);
    let o = 0.5 * sign * (kprime.off_diagonal()[row - 1] + kprime.off_diagonal()[row]);
    let a = k0.diagonal()[row] + shift;
    let d = sign * kprime.diagonal()[row];

    let panels = 32_768_usize;
    let step = (2.0 * PI) / panels as f64;
    let mut sum = 0.0;

    for panel in 0..=panels {
        let theta = panel as f64 * step;
        let integrand = (d - 2.0 * o * theta.cos()) / (a - 2.0 * b * theta.cos()) / (2.0 * PI);
        let weight = if panel == 0 || panel == panels {
            1.0
        } else if panel % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * integrand;
    }

    sum * step / 3.0
}

#[test]
fn closed_frozen_row_density_matches_direct_symbol_quadrature() {
    for row in [4_usize, 32, 128] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let kernel = SignCorrectedResolventTraceKernel::new(row + 2, parity);
            for shift in [0.0_f64, 1.0e-4, 1.0] {
                let closed = kernel.frozen_row_resolvent_density(row, shift).unwrap();
                let quadrature = frozen_symbol_quadrature(row, parity, shift);

                assert_close_scaled(
                    closed,
                    quadrature,
                    2.0e-10,
                    2.0e-12,
                    &format!("row={row} parity={parity:?} shift={shift:e} frozen model"),
                );
            }
        }
    }
}

#[test]
fn frozen_model_requires_an_interior_row() {
    let kernel = SignCorrectedResolventTraceKernel::new(8, ProlateParity::WPlus);

    for row in [0_usize, 7, 8, usize::MAX] {
        assert!(matches!(
            kernel.frozen_row_resolvent_density(row, 1.0e-3),
            Err(ResolventTraceError::InvalidInteriorRow { .. })
        ));
    }
}

#[test]
fn frozen_model_is_positive_in_the_asymptotic_rows_tested() {
    for row in [4_usize, 16, 64, 256] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let kernel = SignCorrectedResolventTraceKernel::new(row + 2, parity);
            for shift in [0.0_f64, 1.0e-6, 1.0e-3, 1.0] {
                let density = kernel.frozen_row_resolvent_density(row, shift).unwrap();
                assert!(
                    density.is_finite() && density > 0.0,
                    "row={row} parity={parity:?} shift={shift:e} density={density:e}"
                );
            }
        }
    }
}
