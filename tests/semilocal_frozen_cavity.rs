use std::f64::consts::PI;

use riemann_ndim_bench::semilocal::{ProlateParity, build_k0, build_kprime_closed};
use riemann_ndim_bench::semilocal_frozen_cavity::frozen_row_cavity_fixed_point;
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

fn frozen_symbol_quadrature(
    block_size: usize,
    row: usize,
    parity: ProlateParity,
    shift: f64,
) -> (f64, f64, f64) {
    let k0 = build_k0(block_size, parity);
    let kprime = build_kprime_closed(block_size, parity);
    let sign = parity.sign_correction();

    let a = k0.diagonal()[row] + shift;
    let b = 0.5 * (k0.off_diagonal()[row - 1] + k0.off_diagonal()[row]);
    let d = sign * kprime.diagonal()[row];
    let o = 0.5 * sign * (kprime.off_diagonal()[row - 1] + kprime.off_diagonal()[row]);

    let panels = 32_768_usize;
    let step = (2.0 * PI) / panels as f64;
    let mut diagonal_sum = 0.0;
    let mut off_diagonal_sum = 0.0;
    let mut weighted_sum = 0.0;

    for panel in 0..=panels {
        let theta = panel as f64 * step;
        let cosine = theta.cos();
        let denominator = a + 2.0 * b * cosine;
        let diagonal_integrand = 1.0 / denominator / (2.0 * PI);
        let off_diagonal_integrand = cosine / denominator / (2.0 * PI);
        let weighted_integrand = (d + 2.0 * o * cosine) / denominator / (2.0 * PI);
        let weight = if panel == 0 || panel == panels {
            1.0
        } else if panel % 2 == 0 {
            2.0
        } else {
            4.0
        };

        diagonal_sum += weight * diagonal_integrand;
        off_diagonal_sum += weight * off_diagonal_integrand;
        weighted_sum += weight * weighted_integrand;
    }

    let factor = step / 3.0;
    (
        diagonal_sum * factor,
        off_diagonal_sum * factor,
        weighted_sum * factor,
    )
}

#[test]
fn frozen_cavity_point_solves_the_schur_fixed_point_and_contracts() {
    for row in [4_usize, 16, 64, 256] {
        let block_size = row + 2;
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            for shift in [0.0_f64, 1.0e-6, 1.0e-3, 1.0] {
                let model = frozen_row_cavity_fixed_point(block_size, row, parity, shift).unwrap();
                let a = model.diagonal_coefficient();
                let b = model.edge_coefficient();
                let q = model.cavity_denominator();

                assert_close_scaled(
                    q,
                    a - b * b / q,
                    2.0e-13,
                    2.0e-14,
                    &format!("row={row} parity={parity:?} shift={shift:e} fixed point"),
                );
                assert!(q.is_finite() && q > 0.0);
                assert!(
                    model.contraction_factor().is_finite()
                        && model.contraction_factor() >= 0.0
                        && model.contraction_factor() < 1.0,
                    "row={row} parity={parity:?} shift={shift:e} contraction={:.16e}",
                    model.contraction_factor()
                );
                assert_close_scaled(
                    model.green_diagonal(),
                    1.0 / (2.0 * q - a),
                    2.0e-12,
                    2.0e-12,
                    &format!("row={row} parity={parity:?} shift={shift:e} Green diagonal"),
                );
                assert_close_scaled(
                    model.green_off_diagonal(),
                    -b * model.green_diagonal() / q,
                    2.0e-13,
                    2.0e-14,
                    &format!("row={row} parity={parity:?} shift={shift:e} Green off diagonal"),
                );
            }
        }
    }
}

#[test]
fn frozen_cavity_green_bands_match_independent_symbol_quadrature() {
    for row in [4_usize, 32, 128] {
        let block_size = row + 2;
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            for shift in [0.0_f64, 1.0e-4, 1.0] {
                let model = frozen_row_cavity_fixed_point(block_size, row, parity, shift).unwrap();
                let (diagonal, off_diagonal, weighted) =
                    frozen_symbol_quadrature(block_size, row, parity, shift);

                assert_close_scaled(
                    model.green_diagonal(),
                    diagonal,
                    2.0e-10,
                    2.0e-12,
                    &format!("row={row} parity={parity:?} shift={shift:e} symbol diagonal"),
                );
                assert_close_scaled(
                    model.green_off_diagonal(),
                    off_diagonal,
                    2.0e-10,
                    2.0e-12,
                    &format!("row={row} parity={parity:?} shift={shift:e} symbol off diagonal"),
                );
                assert_close_scaled(
                    model.weighted_density(),
                    weighted,
                    2.0e-10,
                    2.0e-12,
                    &format!("row={row} parity={parity:?} shift={shift:e} symbol density"),
                );
            }
        }
    }
}

#[test]
fn frozen_cavity_density_matches_the_existing_closed_resolvent_model() {
    for row in [4_usize, 16, 64, 256] {
        let block_size = row + 2;
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let kernel = SignCorrectedResolventTraceKernel::new(block_size, parity);
            for shift in [0.0_f64, 1.0e-6, 1.0e-3, 1.0] {
                let model = frozen_row_cavity_fixed_point(block_size, row, parity, shift).unwrap();
                let closed = kernel.frozen_row_resolvent_density(row, shift).unwrap();

                assert_close_scaled(
                    model.weighted_density(),
                    closed,
                    3.0e-13,
                    3.0e-13,
                    &format!("row={row} parity={parity:?} shift={shift:e} bridge density"),
                );
            }
        }
    }
}

#[test]
fn frozen_cavity_validates_shift_and_interior_row_contracts() {
    assert!(matches!(
        frozen_row_cavity_fixed_point(8, 4, ProlateParity::WPlus, -1.0e-6),
        Err(ResolventTraceError::InvalidShift { .. })
    ));
    assert!(matches!(
        frozen_row_cavity_fixed_point(8, 4, ProlateParity::WPlus, f64::NAN),
        Err(ResolventTraceError::InvalidShift { .. })
    ));

    for row in [0_usize, 7, 8, usize::MAX] {
        assert!(matches!(
            frozen_row_cavity_fixed_point(8, row, ProlateParity::WPlus, 0.0),
            Err(ResolventTraceError::InvalidInteriorRow { .. })
        ));
    }
}
