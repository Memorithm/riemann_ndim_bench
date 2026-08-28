use riemann_ndim_bench::semilocal::ProlateParity;
use riemann_ndim_bench::semilocal_cavity_drift::cavity_drift_factorization;

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
fn drift_splits_exactly_into_edge_and_fixed_point_variation() {
    for block_size in [16_usize, 64, 256] {
        let rows = [
            4_usize,
            (block_size / 4).max(4),
            block_size / 2,
            block_size - 4,
        ];
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            for shift in [0.0_f64, 1.0e-8, 1.0e-4, 1.0] {
                for row in rows {
                    let (transport, factors) =
                        cavity_drift_factorization(block_size, row, parity, shift).unwrap();

                    for value in [
                        factors.left_edge_drift(),
                        factors.left_fixed_point_drift(),
                        factors.right_edge_drift(),
                        factors.right_fixed_point_drift(),
                    ] {
                        assert!(value.is_finite());
                    }

                    assert_close_scaled(
                        transport.left_drift(),
                        factors.reconstructed_left_drift(),
                        8.0e-12,
                        8.0e-13,
                        &format!(
                            "m={block_size} row={row} parity={parity:?} shift={shift:e} left drift"
                        ),
                    );
                    assert_close_scaled(
                        transport.right_drift(),
                        factors.reconstructed_right_drift(),
                        8.0e-12,
                        8.0e-13,
                        &format!(
                            "m={block_size} row={row} parity={parity:?} shift={shift:e} right drift"
                        ),
                    );
                }
            }
        }
    }
}

#[test]
fn transport_multiplier_factors_into_local_contraction_and_corrections() {
    for block_size in [16_usize, 64, 256] {
        let rows = [
            4_usize,
            (block_size / 4).max(4),
            block_size / 2,
            block_size - 4,
        ];
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            for shift in [0.0_f64, 1.0e-8, 1.0e-4, 1.0] {
                for row in rows {
                    let (transport, factors) =
                        cavity_drift_factorization(block_size, row, parity, shift).unwrap();

                    for value in [
                        factors.left_local_contraction(),
                        factors.left_edge_ratio(),
                        factors.left_cavity_ratio(),
                        factors.right_local_contraction(),
                        factors.right_edge_ratio(),
                        factors.right_cavity_ratio(),
                    ] {
                        assert!(value.is_finite() && value > 0.0);
                    }

                    assert!(factors.left_local_contraction() < 1.0);
                    assert!(factors.right_local_contraction() < 1.0);

                    assert_close_scaled(
                        transport.left_transport_factor(),
                        factors.reconstructed_left_transport_factor(),
                        5.0e-13,
                        5.0e-14,
                        &format!(
                            "m={block_size} row={row} parity={parity:?} shift={shift:e} left factor"
                        ),
                    );
                    assert_close_scaled(
                        transport.right_transport_factor(),
                        factors.reconstructed_right_transport_factor(),
                        5.0e-13,
                        5.0e-14,
                        &format!(
                            "m={block_size} row={row} parity={parity:?} shift={shift:e} right factor"
                        ),
                    );
                }
            }
        }
    }
}

#[test]
fn exact_factorization_does_not_assume_uniform_total_contraction() {
    let block_size = 512_usize;
    let row = 256_usize;

    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        let (transport, factors) = cavity_drift_factorization(block_size, row, parity, 0.0).unwrap();

        assert!(factors.left_local_contraction() < 1.0);
        assert!(factors.right_local_contraction() < 1.0);
        assert!(transport.left_transport_factor().is_finite());
        assert!(transport.right_transport_factor().is_finite());
    }
}
