use riemann_ndim_bench::semilocal::{ProlateParity, build_k0};
use riemann_ndim_bench::semilocal_resolvent::SignCorrectedResolventTraceKernel;
use riemann_ndim_bench::semilocal_zero_shift_cavity::{
    zero_shift_left_cavity_closed_form, zero_shift_left_cavity_denominators,
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
fn exact_closed_form_matches_recursive_zero_shift_left_cavity() {
    for block_size in [1_usize, 2, 8, 64, 256] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let kernel = SignCorrectedResolventTraceKernel::new(block_size, parity);
            let cavity = kernel.cavity_green_bands(0.0).unwrap();
            let closed = zero_shift_left_cavity_denominators(block_size, parity);

            for row in 0..block_size {
                assert_close_scaled(
                    cavity.left_denominators()[row],
                    closed[row],
                    2.0e-12,
                    2.0e-13,
                    &format!("m={block_size} row={row} parity={parity:?}"),
                );
            }
        }
    }
}

#[test]
fn closed_form_closes_the_exact_schur_induction() {
    let block_size = 256_usize;

    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        let k0 = build_k0(block_size, parity);

        for row in 0..block_size {
            let current = zero_shift_left_cavity_closed_form(row, parity);

            assert_close_scaled(
                current.reconstructed_diagonal(),
                k0.diagonal()[row],
                5.0e-14,
                5.0e-15,
                &format!("row={row} parity={parity:?} diagonal reconstruction"),
            );

            if row > 0 {
                let previous = zero_shift_left_cavity_closed_form(row - 1, parity);
                let edge = k0.off_diagonal()[row - 1];
                let schur_correction = edge * edge / previous.denominator();

                assert_close_scaled(
                    schur_correction,
                    current.incoming_schur_correction(),
                    8.0e-14,
                    8.0e-15,
                    &format!("row={row} parity={parity:?} incoming Schur correction"),
                );
            } else {
                assert_eq!(current.incoming_schur_correction(), 0.0);
            }
        }
    }
}

#[test]
fn recursive_left_cavity_prefix_is_independent_of_right_boundary() {
    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        let short = SignCorrectedResolventTraceKernel::new(64, parity)
            .cavity_green_bands(0.0)
            .unwrap();
        let long = SignCorrectedResolventTraceKernel::new(256, parity)
            .cavity_green_bands(0.0)
            .unwrap();

        for row in 0..64 {
            assert_close_scaled(
                short.left_denominators()[row],
                long.left_denominators()[row],
                2.0e-14,
                2.0e-15,
                &format!("row={row} parity={parity:?} prefix invariance"),
            );
        }
    }
}
