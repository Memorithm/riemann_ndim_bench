#[path = "../src/quadrature.rs"]
mod quadrature;

use faer::{Mat, Side, linalg::solvers::SelfAdjointEigen};
use quadrature::GaussLegendreUnit;
use std::f64::consts::PI;

#[test]
fn reproduces_published_boundary_contributions() {
    let order = 128;
    let q = GaussLegendreUnit::new(order).unwrap();
    let x = q.nodes();
    let w = q.weights();
    let c = 2.0 * PI;
    let matrix = Mat::from_fn(order, order, |i, j| {
        2.0 * (w[i] * w[j]).sqrt() * (c * x[i] * x[j]).cos()
    });
    let decomposition = SelfAdjointEigen::new(matrix.as_ref(), Side::Lower).unwrap();
    let eigenvalues = decomposition.S().column_vector();
    let eigenvectors = decomposition.U();
    let mut indices = (0..order).collect::<Vec<_>>();
    indices.sort_by(|&a, &b| eigenvalues[b].abs().total_cmp(&eigenvalues[a].abs()));

    let expected = [11.9719, 8.77574, 2.20528, 0.0433983, 0.000125459];
    let mut total = 0.0;

    for (mode, &target) in expected.iter().enumerate() {
        let k = indices[mode];
        let lambda = eigenvalues[k];
        let boundary = (0..order)
            .map(|i| w[i].sqrt() * eigenvectors[(i, k)] * (c * x[i]).cos())
            .sum::<f64>()
            * 2.0
            / lambda;
        let lambda2 = lambda * lambda;
        let contribution = lambda2 / (1.0 - lambda2) * boundary * boundary;
        total += contribution;
        assert!(
            (contribution - target).abs() < 8.0e-5 * target.abs().max(1.0),
            "mode {mode}: got {contribution:.10}, expected {target:.10}"
        );
    }

    let sixth = indices[5];
    let lambda = eigenvalues[sixth];
    let boundary = (0..order)
        .map(|i| w[i].sqrt() * eigenvectors[(i, sixth)] * (c * x[i]).cos())
        .sum::<f64>()
        * 2.0
        / lambda;
    let lambda2 = lambda * lambda;
    total += lambda2 / (1.0 - lambda2) * boundary * boundary;
    assert!((total - 22.9965).abs() < 2.0e-3, "total={total:.10}");
}
