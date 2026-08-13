#[path = "support/prolate_basis.rs"]
mod prolate_basis;
#[path = "../src/quadrature.rs"]
mod quadrature;

use prolate_basis::Basis;
use quadrature::GaussLegendreUnit;
use riemann_ndim_bench::toeplitz::{LogLattice, SymmetricToeplitz};
use std::f64::consts::LN_2;

fn c_n(basis: &Basis, mode: usize, rho: f64, integration: &GaussLegendreUnit) -> f64 {
    if rho == 1.0 {
        return 0.0;
    }

    let lower = rho.recip();
    let span = 1.0 - lower;
    let integral = integration
        .nodes()
        .iter()
        .zip(integration.weights().iter())
        .map(|(&u, &weight)| {
            let x = lower + span * u;
            let left = x * basis.derivative(mode, x);
            let right_x = rho * x;
            let right = right_x * basis.derivative(mode, right_x);
            weight * left * right
        })
        .sum::<f64>()
        * span;

    rho.sqrt() * integral + rho.powf(-1.5) * basis.derivative(mode, lower) * basis.value(mode, 1.0)
        - rho.powf(1.5) * basis.value(mode, 1.0) * basis.derivative(mode, rho)
}

fn q_epsilon(basis: &Basis, rho: f64, integration: &GaussLegendreUnit) -> f64 {
    basis
        .lambdas
        .iter()
        .enumerate()
        .map(|(mode, &lambda)| {
            let lambda2 = lambda * lambda;
            lambda2 / (1.0 - lambda2) * c_n(basis, mode, rho, integration)
        })
        .sum()
}

#[test]
fn reproduces_published_qepsilon_toeplitz_spectrum() {
    // Equation (100) is evaluated with the first eleven prolate modes.
    // Connes--Consani Appendix F bounds the omitted remainder uniformly on [1, 2].
    let basis = Basis::compute(64, 11);
    let integration = GaussLegendreUnit::new(24).unwrap();
    let epsilon_prime = basis.epsilon_prime();
    assert!((epsilon_prime - 22.9965).abs() < 3.0e-4);

    let lattice = LogLattice::new(LN_2, 1.0e-3).unwrap();
    let operator = SymmetricToeplitz::sample_normalized_kernel(lattice, |x| {
        let rho = x.exp();
        q_epsilon(&basis, rho, &integration) / (2.0 * epsilon_prime)
    })
    .unwrap();

    let eigenvalues = operator.eigenvalues().unwrap();
    let largest = eigenvalues[eigenvalues.len() - 1];
    let second = eigenvalues[eigenvalues.len() - 2];

    // Section 6 reports approximately 1.05177 and 0.687925 for omega=10^-3.
    assert!((largest - 1.05177).abs() < 8.0e-5, "largest={largest:.9}");
    assert!((second - 0.687925).abs() < 8.0e-5, "second={second:.9}");
}
