use riemann_ndim_bench::toeplitz::{LogLattice, SymmetricToeplitz};
use std::f64::consts::LN_2;

const EPS: f64 = 1.0e-12;

#[test]
fn paper_scale_has_expected_dimension() {
    let lattice = LogLattice::new(LN_2, 1.0e-3).unwrap();
    assert_eq!(lattice.max_index(), 693);
    assert_eq!(lattice.dimension(), 694);
    assert!((lattice.q() - 1.0e-3_f64.exp()).abs() < EPS);
}

#[test]
fn sampled_kernel_uses_omega_weight() {
    let lattice = LogLattice::new(0.25, 0.1).unwrap();
    let matrix = SymmetricToeplitz::sample_normalized_kernel(lattice, |x| 2.0 + x).unwrap();
    let expected = [0.2, 0.21, 0.22];
    for (&actual, &target) in matrix.first_row().iter().zip(expected.iter()) {
        assert!((actual - target).abs() < EPS);
    }
}

#[test]
fn matrix_free_application_matches_hand_calculation() {
    let matrix = SymmetricToeplitz::from_first_row(vec![2.0, 1.0, 0.5]).unwrap();
    let result = matrix.apply(&[1.0, 2.0, -1.0]).unwrap();
    let expected = [3.5, 4.0, 0.5];
    for (&actual, &target) in result.iter().zip(expected.iter()) {
        assert!((actual - target).abs() < EPS);
    }
}

#[test]
fn toeplitz_form_is_not_positive_by_construction() {
    let matrix = SymmetricToeplitz::from_first_row(vec![0.0, 1.0]).unwrap();
    let positive_direction = matrix.quadratic_form(&[1.0, 1.0]).unwrap();
    let negative_direction = matrix.quadratic_form(&[1.0, -1.0]).unwrap();
    assert!(positive_direction > 0.0);
    assert!(negative_direction < 0.0);
}

#[test]
fn two_by_two_spectrum_matches_closed_form() {
    let matrix = SymmetricToeplitz::from_first_row(vec![2.0, 1.0]).unwrap();
    let eigenvalues = matrix.eigenvalues().unwrap();
    assert_eq!(eigenvalues.len(), 2);
    assert!((eigenvalues[0] - 1.0).abs() < EPS);
    assert!((eigenvalues[1] - 3.0).abs() < EPS);
    assert!((matrix.largest_eigenvalue().unwrap() - 3.0).abs() < EPS);
}

#[test]
fn constant_sampled_kernel_has_rank_one_spectrum() {
    let lattice = LogLattice::new(0.25, 0.1).unwrap();
    let matrix = SymmetricToeplitz::sample_normalized_kernel(lattice, |_| 1.0).unwrap();
    let eigenvalues = matrix.eigenvalues().unwrap();
    let largest = *eigenvalues.last().unwrap();
    assert!((largest - lattice.dimension() as f64 * lattice.omega()).abs() < EPS);
    for &eigenvalue in &eigenvalues[..eigenvalues.len() - 1] {
        assert!(eigenvalue.abs() < 10.0 * EPS);
    }
}
