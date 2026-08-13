use crate::quadrature::GaussLegendreUnit;
use faer::{Mat, Side, linalg::solvers::SelfAdjointEigen};
use std::f64::consts::PI;

pub struct Basis {
    pub nodes: Vec<f64>,
    pub weights: Vec<f64>,
    pub lambdas: Vec<f64>,
    pub samples: Vec<Vec<f64>>,
}

impl Basis {
    pub fn compute(order: usize, count: usize) -> Self {
        let q = GaussLegendreUnit::new(order).unwrap();
        let nodes = q.nodes().to_vec();
        let weights = q.weights().to_vec();
        let c = 2.0 * PI;
        let matrix = Mat::from_fn(order, order, |i, j| {
            2.0 * (weights[i] * weights[j]).sqrt() * (c * nodes[i] * nodes[j]).cos()
        });
        let decomposition = SelfAdjointEigen::new(matrix.as_ref(), Side::Lower).unwrap();
        let diagonal = decomposition.S().column_vector();
        let vectors = decomposition.U();
        let mut indices = (0..order).collect::<Vec<_>>();
        indices.sort_by(|&a, &b| diagonal[b].abs().total_cmp(&diagonal[a].abs()));
        indices.truncate(count);
        let lambdas = indices.iter().map(|&k| diagonal[k]).collect();
        let samples = indices
            .iter()
            .map(|&k| (0..order).map(|i| vectors[(i, k)] / weights[i].sqrt()).collect())
            .collect();
        Self {
            nodes,
            weights,
            lambdas,
            samples,
        }
    }

    pub fn value(&self, mode: usize, y: f64) -> f64 {
        let c = 2.0 * PI;
        let sum = self
            .nodes
            .iter()
            .zip(self.weights.iter())
            .zip(self.samples[mode].iter())
            .map(|((&x, &w), &f)| w * f * (c * x * y).cos())
            .sum::<f64>();
        2.0 * sum / self.lambdas[mode]
    }

    pub fn derivative(&self, mode: usize, y: f64) -> f64 {
        let c = 2.0 * PI;
        let sum = self
            .nodes
            .iter()
            .zip(self.weights.iter())
            .zip(self.samples[mode].iter())
            .map(|((&x, &w), &f)| w * f * x * (c * x * y).sin())
            .sum::<f64>();
        -2.0 * c * sum / self.lambdas[mode]
    }

    pub fn epsilon_prime(&self) -> f64 {
        (0..self.lambdas.len())
            .map(|mode| {
                let lambda2 = self.lambdas[mode] * self.lambdas[mode];
                let boundary = self.value(mode, 1.0);
                lambda2 / (1.0 - lambda2) * boundary * boundary
            })
            .sum()
    }
}
