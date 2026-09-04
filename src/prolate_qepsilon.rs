//! Source-locked archimedean prolate / `Q epsilon` construction.
//!
//! This module promotes the Connes--Consani archimedean source regression out
//! of integration-test scaffolding so that later Riemann-specific Weil-bridge
//! work can reuse the same mathematical object. The formulas follow
//! `docs/Q_EPSILON_SPEC.md` and are not fitted to the published spectrum.
//!
//! Reproducing this kernel validates the archimedean source calculation. It
//! does not identify any finite eigenvalue with a zeta zero and does not prove
//! the Riemann hypothesis.

use std::f64::consts::PI;
use std::fmt;

use faer::linalg::solvers::{EvdError, SelfAdjointEigen};
use faer::{Mat, Side};

use crate::quadrature::{GaussLegendreUnit, QuadratureError};

/// Error returned while constructing or evaluating the source-locked kernel.
#[derive(Debug)]
pub enum ProlateQepsilonError {
    /// The requested number of prolate modes must lie in `1..=order`.
    InvalidModeCount { order: usize, count: usize },
    /// `Q epsilon` is evaluated from the source formula on `rho >= 1`.
    InvalidRho { rho: f64 },
    /// Gauss--Legendre quadrature construction failed.
    Quadrature(QuadratureError),
    /// The symmetric prolate eigensolver failed.
    Eigensolver(EvdError),
}

impl fmt::Display for ProlateQepsilonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidModeCount { order, count } => write!(
                f,
                "prolate mode count must satisfy 1 <= count <= order: order={order}, count={count}"
            ),
            Self::InvalidRho { rho } => {
                write!(
                    f,
                    "Q epsilon source evaluation requires finite rho >= 1: {rho}"
                )
            }
            Self::Quadrature(error) => write!(f, "quadrature construction failed: {error:?}"),
            Self::Eigensolver(error) => write!(f, "prolate eigensolver failed: {error:?}"),
        }
    }
}

impl std::error::Error for ProlateQepsilonError {}

impl From<QuadratureError> for ProlateQepsilonError {
    fn from(value: QuadratureError) -> Self {
        Self::Quadrature(value)
    }
}

impl From<EvdError> for ProlateQepsilonError {
    fn from(value: EvdError) -> Self {
        Self::Eigensolver(value)
    }
}

/// Even prolate basis used by the archimedean `epsilon` / `Q epsilon` series.
///
/// The discretized integral operator is
///
/// `2 * integral_0^1 cos(2 pi x y) f(y) dy`,
///
/// evaluated with Gauss--Legendre quadrature. Eigenmodes are ordered by
/// decreasing absolute eigenvalue, matching the source regression convention.
#[derive(Clone, Debug)]
pub struct ProlateBasis {
    nodes: Vec<f64>,
    weights: Vec<f64>,
    eigenvalues: Vec<f64>,
    samples: Vec<Vec<f64>>,
}

impl ProlateBasis {
    /// Construct the first `count` even prolate modes using `order` quadrature
    /// nodes on `[0,1]`.
    pub fn compute(order: usize, count: usize) -> Result<Self, ProlateQepsilonError> {
        if count == 0 || count > order {
            return Err(ProlateQepsilonError::InvalidModeCount { order, count });
        }

        let quadrature = GaussLegendreUnit::new(order)?;
        let nodes = quadrature.nodes().to_vec();
        let weights = quadrature.weights().to_vec();
        let c = 2.0 * PI;
        let matrix = Mat::from_fn(order, order, |i, j| {
            2.0 * (weights[i] * weights[j]).sqrt() * (c * nodes[i] * nodes[j]).cos()
        });
        let decomposition = SelfAdjointEigen::new(matrix.as_ref(), Side::Lower)?;
        let diagonal = decomposition.S().column_vector();
        let vectors = decomposition.U();

        let mut indices = (0..order).collect::<Vec<_>>();
        indices.sort_by(|&left, &right| diagonal[right].abs().total_cmp(&diagonal[left].abs()));
        indices.truncate(count);

        let eigenvalues = indices.iter().map(|&index| diagonal[index]).collect();
        let samples = indices
            .iter()
            .map(|&column| {
                (0..order)
                    .map(|row| vectors[(row, column)] / weights[row].sqrt())
                    .collect()
            })
            .collect();

        Ok(Self {
            nodes,
            weights,
            eigenvalues,
            samples,
        })
    }

    /// Number of retained source modes.
    #[inline]
    pub fn mode_count(&self) -> usize {
        self.eigenvalues.len()
    }

    /// Prolate integral-operator eigenvalues in decreasing absolute value.
    #[inline]
    pub fn eigenvalues(&self) -> &[f64] {
        &self.eigenvalues
    }

    /// Quadrature nodes used to construct the basis.
    #[inline]
    pub fn nodes(&self) -> &[f64] {
        &self.nodes
    }

    /// Quadrature weights used to construct the basis.
    #[inline]
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// Evaluate the source continuation of one retained prolate mode.
    pub fn value(&self, mode: usize, y: f64) -> f64 {
        let c = 2.0 * PI;
        let sum = self
            .nodes
            .iter()
            .zip(self.weights.iter())
            .zip(self.samples[mode].iter())
            .map(|((&x, &weight), &sample)| weight * sample * (c * x * y).cos())
            .sum::<f64>();
        2.0 * sum / self.eigenvalues[mode]
    }

    /// Evaluate the derivative of the source continuation with respect to `y`.
    pub fn derivative(&self, mode: usize, y: f64) -> f64 {
        let c = 2.0 * PI;
        let sum = self
            .nodes
            .iter()
            .zip(self.weights.iter())
            .zip(self.samples[mode].iter())
            .map(|((&x, &weight), &sample)| weight * sample * x * (c * x * y).sin())
            .sum::<f64>();
        -2.0 * c * sum / self.eigenvalues[mode]
    }

    /// Published boundary-series contribution of one retained mode to
    /// `epsilon'(1+)`.
    pub fn epsilon_prime_contribution(&self, mode: usize) -> f64 {
        let lambda = self.eigenvalues[mode];
        let lambda2 = lambda * lambda;
        let boundary = self.value(mode, 1.0);
        lambda2 / (1.0 - lambda2) * boundary * boundary
    }

    /// Truncated source series for `epsilon'(1+)` over the retained modes.
    pub fn epsilon_prime(&self) -> f64 {
        (0..self.mode_count())
            .map(|mode| self.epsilon_prime_contribution(mode))
            .sum()
    }
}

/// Reusable source-locked `Q epsilon` evaluator.
#[derive(Clone, Debug)]
pub struct ProlateQepsilonKernel {
    basis: ProlateBasis,
    integration: GaussLegendreUnit,
    epsilon_prime: f64,
}

impl ProlateQepsilonKernel {
    /// Construct the source kernel from independent prolate and integration
    /// quadrature orders.
    pub fn new(
        prolate_order: usize,
        mode_count: usize,
        integration_order: usize,
    ) -> Result<Self, ProlateQepsilonError> {
        let basis = ProlateBasis::compute(prolate_order, mode_count)?;
        let integration = GaussLegendreUnit::new(integration_order)?;
        let epsilon_prime = basis.epsilon_prime();
        Ok(Self {
            basis,
            integration,
            epsilon_prime,
        })
    }

    /// Retained prolate basis.
    #[inline]
    pub fn basis(&self) -> &ProlateBasis {
        &self.basis
    }

    /// Truncated source value of `epsilon'(1+)` used for normalization.
    #[inline]
    pub fn epsilon_prime(&self) -> f64 {
        self.epsilon_prime
    }

    /// Evaluate the source coefficient `C_n(rho)` from equation (99).
    pub fn c_n(&self, mode: usize, rho: f64) -> Result<f64, ProlateQepsilonError> {
        checked_rho(rho)?;
        if rho == 1.0 {
            return Ok(0.0);
        }

        let lower = rho.recip();
        let span = 1.0 - lower;
        let integral = self
            .integration
            .nodes()
            .iter()
            .zip(self.integration.weights().iter())
            .map(|(&unit_x, &weight)| {
                let x = lower + span * unit_x;
                let left = x * self.basis.derivative(mode, x);
                let right_x = rho * x;
                let right = right_x * self.basis.derivative(mode, right_x);
                weight * left * right
            })
            .sum::<f64>()
            * span;

        Ok(rho.sqrt() * integral
            + rho.powf(-1.5) * self.basis.derivative(mode, lower) * self.basis.value(mode, 1.0)
            - rho.powf(1.5) * self.basis.value(mode, 1.0) * self.basis.derivative(mode, rho))
    }

    /// Evaluate the retained source series `Q epsilon(rho)` for `rho >= 1`.
    pub fn q_epsilon(&self, rho: f64) -> Result<f64, ProlateQepsilonError> {
        checked_rho(rho)?;
        let mut total = 0.0;
        for (mode, &lambda) in self.basis.eigenvalues().iter().enumerate() {
            let lambda2 = lambda * lambda;
            total += lambda2 / (1.0 - lambda2) * self.c_n(mode, rho)?;
        }
        Ok(total)
    }

    /// Evaluate the normalized kernel
    /// `chi(log rho) = Q epsilon(rho) / (2 epsilon'(1+))`.
    pub fn normalized_q_epsilon(&self, rho: f64) -> Result<f64, ProlateQepsilonError> {
        Ok(self.q_epsilon(rho)? / (2.0 * self.epsilon_prime))
    }

    /// Evaluate the even logarithmic kernel used by the Toeplitz lattice.
    ///
    /// `log_distance` is mapped to `rho = exp(|log_distance|)` using the source
    /// symmetry around `rho=1`.
    pub fn normalized_log_kernel(&self, log_distance: f64) -> Result<f64, ProlateQepsilonError> {
        if !log_distance.is_finite() {
            return Err(ProlateQepsilonError::InvalidRho { rho: log_distance });
        }
        self.normalized_q_epsilon(log_distance.abs().exp())
    }
}

fn checked_rho(rho: f64) -> Result<(), ProlateQepsilonError> {
    if rho.is_finite() && rho >= 1.0 {
        Ok(())
    } else {
        Err(ProlateQepsilonError::InvalidRho { rho })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_mode_counts_and_rho() {
        assert!(matches!(
            ProlateBasis::compute(8, 0),
            Err(ProlateQepsilonError::InvalidModeCount { .. })
        ));
        assert!(matches!(
            ProlateBasis::compute(8, 9),
            Err(ProlateQepsilonError::InvalidModeCount { .. })
        ));

        let kernel = ProlateQepsilonKernel::new(16, 4, 8).unwrap();
        assert!(matches!(
            kernel.q_epsilon(0.5),
            Err(ProlateQepsilonError::InvalidRho { .. })
        ));
        assert!(matches!(
            kernel.q_epsilon(f64::NAN),
            Err(ProlateQepsilonError::InvalidRho { .. })
        ));
    }

    #[test]
    fn q_epsilon_vanishes_at_the_source_boundary() {
        let kernel = ProlateQepsilonKernel::new(32, 6, 12).unwrap();
        assert_eq!(kernel.q_epsilon(1.0).unwrap(), 0.0);
        assert_eq!(kernel.normalized_log_kernel(0.0).unwrap(), 0.0);
    }
}
