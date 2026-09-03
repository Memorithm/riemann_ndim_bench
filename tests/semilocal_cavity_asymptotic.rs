use riemann_ndim_bench::semilocal::{ProlateParity, build_k0};
use riemann_ndim_bench::semilocal_cavity_asymptotic::{
    frozen_cavity_denominator_second_order, frozen_contraction_asymptotic_diagnostic,
    frozen_contraction_second_order, frozen_diagonal_second_order, frozen_edge_second_order,
    frozen_soft_gap_second_order,
};
use riemann_ndim_bench::semilocal_frozen_cavity::frozen_row_cavity_fixed_point;

#[test]
fn asymptotic_approximants_require_positive_degree() {
    assert_eq!(frozen_diagonal_second_order(0), None);
    assert_eq!(frozen_edge_second_order(0), None);
    assert_eq!(frozen_soft_gap_second_order(0), None);
    assert_eq!(frozen_cavity_denominator_second_order(0), None);
    assert_eq!(frozen_contraction_second_order(0), None);
}

#[test]
fn source_derived_local_series_track_exact_k0_coefficients() {
    for row in [64_usize, 128, 256] {
        let block_size = row + 2;
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let degree = parity.degree(row);
            let degree_f = degree as f64;
            let k0 = build_k0(block_size, parity);
            let diagonal = k0.diagonal()[row];
            let edge = 0.5 * (k0.off_diagonal()[row - 1] + k0.off_diagonal()[row]);
            let gap = diagonal - 2.0 * edge;
            let fixed = frozen_row_cavity_fixed_point(block_size, row, parity, 0.0).unwrap();

            let diagonal_residual =
                diagonal - frozen_diagonal_second_order(degree).unwrap();
            let edge_residual = edge - frozen_edge_second_order(degree).unwrap();
            let gap_residual = gap - frozen_soft_gap_second_order(degree).unwrap();
            let cavity_residual = fixed.cavity_denominator()
                - frozen_cavity_denominator_second_order(degree).unwrap();

            assert!(diagonal_residual.abs() * degree_f.powi(3) < 0.01);
            assert!(edge_residual.abs() * degree_f.powi(3) < 0.03);
            assert!(gap_residual.abs() * degree_f.powi(3) < 0.05);
            assert!(cavity_residual.abs() * degree_f.powi(3) < 0.03);
        }
    }
}

#[test]
fn frozen_contraction_has_the_predicted_cumulative_soft_edge_scaling() {
    for row in [64_usize, 128, 256, 512] {
        let block_size = row + 2;
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let diagnostic =
                frozen_contraction_asymptotic_diagnostic(block_size, row, parity).unwrap();
            let degree = diagnostic.degree() as f64;

            assert!(diagnostic.exact_contraction() > 0.0);
            assert!(diagnostic.exact_contraction() < 1.0);
            assert!((diagnostic.scaled_first_order_gap() - 1.0).abs() < 0.01);
            assert!((diagnostic.scaled_second_order_remainder() - 0.75).abs() < 0.01);
            assert!(diagnostic.residual().abs() * degree.powi(3) < 1.5);
        }
    }
}
