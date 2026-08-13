#[path = "../src/quadrature.rs"]
mod quadrature;

use faer::{Mat, Side, linalg::solvers::SelfAdjointEigen};
use quadrature::GaussLegendreUnit;
use std::f64::consts::PI;

fn computed_even_eigenvalues(order: usize, count: usize) -> Vec<f64> {
    let q = GaussLegendreUnit::new(order).unwrap();
    let x = q.nodes();
    let w = q.weights();
    let matrix = Mat::from_fn(order, order, |i, j| {
        2.0 * (w[i] * w[j]).sqrt() * (2.0 * PI * x[i] * x[j]).cos()
    });
    let decomposition = SelfAdjointEigen::new(matrix.as_ref(), Side::Lower).unwrap();
    let diagonal = decomposition.S().column_vector();
    let mut values = (0..order).map(|i| diagonal[i]).collect::<Vec<_>>();
    values.sort_by(|left, right| right.abs().total_cmp(&left.abs()));
    values.truncate(count);
    values
}

#[test]
fn reproduces_published_even_prolate_eigenvalues() {
    let actual = computed_even_eigenvalues(128, 6);
    let expected = [
        0.999971,
        -0.979485,
        0.524086,
        -0.0589766,
        0.00273233,
        -0.0000762914,
    ];

    for (index, (&got, &target)) in actual.iter().zip(expected.iter()).enumerate() {
        assert!(
            (got - target).abs() < 8.0e-6,
            "mode {index}: got {got:.12e}, expected {target:.12e}"
        );
    }
}
