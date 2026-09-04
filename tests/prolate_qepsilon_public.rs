use std::f64::consts::LN_2;

use riemann_ndim_bench::prolate_qepsilon::{ProlateBasis, ProlateQepsilonKernel};
use riemann_ndim_bench::toeplitz::{LogLattice, SymmetricToeplitz};

#[test]
fn public_prolate_basis_reproduces_source_modes_and_boundary_series() {
    let basis = ProlateBasis::compute(128, 6).unwrap();
    let expected_eigenvalues = [
        0.999971,
        -0.979485,
        0.524086,
        -0.0589766,
        0.00273233,
        -0.0000762914,
    ];
    let expected_contributions = [11.9719, 8.77574, 2.20528, 0.0433983, 0.000125459];

    for (mode, (&actual, &expected)) in basis
        .eigenvalues()
        .iter()
        .zip(expected_eigenvalues.iter())
        .enumerate()
    {
        assert!(
            (actual - expected).abs() < 8.0e-6,
            "mode {mode}: eigenvalue={actual:.12e} expected={expected:.12e}"
        );
    }

    for (mode, &expected) in expected_contributions.iter().enumerate() {
        let actual = basis.epsilon_prime_contribution(mode);
        assert!(
            (actual - expected).abs() < 8.0e-5 * expected.abs().max(1.0),
            "mode {mode}: contribution={actual:.10} expected={expected:.10}"
        );
    }

    assert!(
        (basis.epsilon_prime() - 22.9965).abs() < 2.0e-3,
        "epsilon_prime={:.10}",
        basis.epsilon_prime()
    );
}

#[test]
fn public_qepsilon_kernel_reproduces_archimedean_toeplitz_benchmark() {
    let kernel = ProlateQepsilonKernel::new(64, 11, 24).unwrap();
    assert!((kernel.epsilon_prime() - 22.9965).abs() < 3.0e-4);
    assert_eq!(kernel.q_epsilon(1.0).unwrap(), 0.0);

    let lattice = LogLattice::new(LN_2, 1.0e-3).unwrap();
    let operator = SymmetricToeplitz::sample_normalized_kernel(lattice, |x| {
        kernel.normalized_log_kernel(x).unwrap()
    })
    .unwrap();

    let eigenvalues = operator.eigenvalues().unwrap();
    let largest = eigenvalues[eigenvalues.len() - 1];
    let second = eigenvalues[eigenvalues.len() - 2];

    assert!((largest - 1.05177).abs() < 8.0e-5, "largest={largest:.9}");
    assert!((second - 0.687925).abs() < 8.0e-5, "second={second:.9}");
}
