use std::f64::consts::PI;

use riemann_ndim_bench::semilocal::{ProlateParity, build_k0};
use riemann_ndim_bench::semilocal_zero_shift_cavity::{
    zero_shift_left_cavity_closed_form, zero_shift_left_cavity_denominators,
};
use riemann_ndim_bench::semilocal_zero_shift_response::{
    zero_shift_left_cavity_relative_derivative, zero_shift_left_cavity_response,
    zero_shift_left_cavity_shift_derivatives,
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
fn exact_shift_derivative_matches_differentiated_schur_recurrence() {
    for block_size in [1_usize, 2, 8, 64, 256] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let k0 = build_k0(block_size, parity);
            let left = zero_shift_left_cavity_denominators(block_size, parity);
            let closed = zero_shift_left_cavity_shift_derivatives(block_size, parity);
            let mut recursive = Vec::with_capacity(block_size);

            if block_size > 0 {
                let mut derivative = 1.0;
                recursive.push(derivative);
                for (edge, previous_left) in k0.off_diagonal().iter().zip(left.iter()) {
                    derivative = 1.0 + edge * edge / (previous_left * previous_left) * derivative;
                    recursive.push(derivative);
                }
            }

            for (row, expected) in closed.iter().copied().enumerate() {
                assert_close_scaled(
                    recursive[row],
                    expected,
                    3.0e-12,
                    3.0e-13,
                    &format!("m={block_size} row={row} parity={parity:?} shift derivative"),
                );
            }
        }
    }
}

#[test]
fn relative_first_shift_response_is_universal() {
    let expected = 8.0 * PI / 3.0;
    assert_close_scaled(
        zero_shift_left_cavity_relative_derivative(),
        expected,
        1.0e-15,
        1.0e-15,
        "exported relative derivative",
    );

    for row in [0_usize, 1, 2, 8, 64, 256] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let response = zero_shift_left_cavity_response(row, parity);
            assert_close_scaled(
                response.relative_shift_derivative(),
                expected,
                2.0e-15,
                2.0e-15,
                &format!("row={row} parity={parity:?} relative derivative"),
            );
        }
    }
}

#[test]
fn constant_gap_closes_the_shift_response_induction() {
    let expected_gap = 3.0 / (8.0 * PI);

    for row in [0_usize, 1, 2, 8, 64, 256] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let closed = zero_shift_left_cavity_closed_form(row, parity);
            let gap = closed.denominator() - closed.incoming_schur_correction();
            assert_close_scaled(
                gap,
                expected_gap,
                3.0e-14,
                3.0e-15,
                &format!("row={row} parity={parity:?} induction gap"),
            );
        }
    }
}

#[test]
fn closed_derivative_has_the_direct_degree_formula() {
    for row in [0_usize, 1, 2, 8, 64, 256] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let response = zero_shift_left_cavity_response(row, parity);
            let degree = response.degree() as f64;
            let expected =
                (2.0 * degree + 1.0) * (2.0 * degree + 3.0) / (3.0 * (4.0 * degree + 1.0));
            assert_close_scaled(
                response.shift_derivative(),
                expected,
                3.0e-15,
                3.0e-15,
                &format!("row={row} parity={parity:?} degree formula"),
            );
        }
    }
}
