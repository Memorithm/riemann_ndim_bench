use riemann_ndim_bench::semilocal::ProlateParity;
use riemann_ndim_bench::semilocal_cavity_transport::cavity_error_transport;
use riemann_ndim_bench::semilocal_resolvent::ResolventTraceError;

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
fn finite_cavity_errors_split_exactly_into_transport_and_drift() {
    for block_size in [16_usize, 64, 256] {
        let rows = [2_usize, block_size / 4, block_size / 2, block_size - 3];
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            for shift in [0.0_f64, 1.0e-8, 1.0e-4, 1.0] {
                for row in rows {
                    let transport =
                        cavity_error_transport(block_size, row, parity, shift).unwrap();

                    for value in [
                        transport.left_error(),
                        transport.left_previous_error(),
                        transport.left_transport_factor(),
                        transport.left_drift(),
                        transport.right_error(),
                        transport.right_next_error(),
                        transport.right_transport_factor(),
                        transport.right_drift(),
                    ] {
                        assert!(value.is_finite());
                    }

                    assert!(transport.left_transport_factor() > 0.0);
                    assert!(transport.right_transport_factor() > 0.0);

                    assert_close_scaled(
                        transport.left_error(),
                        transport.reconstructed_left_error(),
                        5.0e-12,
                        5.0e-13,
                        &format!(
                            "m={block_size} row={row} parity={parity:?} shift={shift:e} left transport"
                        ),
                    );
                    assert_close_scaled(
                        transport.right_error(),
                        transport.reconstructed_right_error(),
                        5.0e-12,
                        5.0e-13,
                        &format!(
                            "m={block_size} row={row} parity={parity:?} shift={shift:e} right transport"
                        ),
                    );
                }
            }
        }
    }
}

#[test]
fn transport_identity_resolves_left_and_right_boundary_influence_separately() {
    let block_size = 128_usize;
    let row = 48_usize;

    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        for shift in [0.0_f64, 1.0e-6, 1.0e-2, 1.0] {
            let transport = cavity_error_transport(block_size, row, parity, shift).unwrap();

            let left_propagated =
                transport.left_transport_factor() * transport.left_previous_error();
            let right_propagated =
                transport.right_transport_factor() * transport.right_next_error();

            assert_close_scaled(
                transport.left_error() - transport.left_drift(),
                left_propagated,
                5.0e-12,
                5.0e-13,
                &format!("parity={parity:?} shift={shift:e} isolated left boundary term"),
            );
            assert_close_scaled(
                transport.right_error() - transport.right_drift(),
                right_propagated,
                5.0e-12,
                5.0e-13,
                &format!("parity={parity:?} shift={shift:e} isolated right boundary term"),
            );
        }
    }
}

#[test]
fn cavity_transport_requires_two_interior_neighbors() {
    let block_size = 16_usize;
    for row in [0_usize, 1, 14, 15, 16, usize::MAX] {
        assert!(matches!(
            cavity_error_transport(block_size, row, ProlateParity::WPlus, 0.0),
            Err(ResolventTraceError::InvalidInteriorRow { .. })
        ));
    }
}

#[test]
fn cavity_transport_preserves_shift_validation() {
    for shift in [-1.0e-6_f64, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            cavity_error_transport(16, 8, ProlateParity::WPlus, shift),
            Err(ResolventTraceError::InvalidShift { .. })
        ));
    }
}
