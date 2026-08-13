use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuadratureError {
    OrderTooSmall,
    DidNotConverge,
}

#[derive(Debug, Clone)]
pub struct GaussLegendreUnit {
    nodes: Vec<f64>,
    weights: Vec<f64>,
}

impl GaussLegendreUnit {
    pub fn new(order: usize) -> Result<Self, QuadratureError> {
        if order < 2 {
            return Err(QuadratureError::OrderTooSmall);
        }

        let mut nodes = vec![0.0; order];
        let mut weights = vec![0.0; order];
        let half = order.div_ceil(2);

        for i in 0..half {
            let mut z = (PI * (i as f64 + 0.75) / (order as f64 + 0.5)).cos();
            let mut converged = false;

            for _ in 0..64 {
                let (pn, pnm1) = legendre_pair(order, z);
                let derivative = order as f64 * (z * pn - pnm1) / (z * z - 1.0);
                let next = z - pn / derivative;
                if (next - z).abs() <= 4.0 * f64::EPSILON * next.abs().max(1.0) {
                    z = next;
                    converged = true;
                    break;
                }
                z = next;
            }

            if !converged {
                return Err(QuadratureError::DidNotConverge);
            }

            let (pn, pnm1) = legendre_pair(order, z);
            let derivative = order as f64 * (z * pn - pnm1) / (z * z - 1.0);
            let full_weight = 2.0 / ((1.0 - z * z) * derivative * derivative);
            let weight = 0.5 * full_weight;

            let left = i;
            let right = order - 1 - i;
            nodes[left] = 0.5 * (1.0 - z);
            weights[left] = weight;
            nodes[right] = 0.5 * (1.0 + z);
            weights[right] = weight;
        }

        Ok(Self { nodes, weights })
    }

    pub fn nodes(&self) -> &[f64] {
        &self.nodes
    }

    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    pub fn order(&self) -> usize {
        self.nodes.len()
    }

    pub fn integrate(&self, mut f: impl FnMut(f64) -> f64) -> f64 {
        self.nodes
            .iter()
            .zip(self.weights.iter())
            .map(|(&x, &w)| w * f(x))
            .sum()
    }
}

fn legendre_pair(order: usize, x: f64) -> (f64, f64) {
    if order == 0 {
        return (1.0, 1.0);
    }

    let mut pnm1 = 1.0;
    let mut pn = x;
    if order == 1 {
        return (pn, pnm1);
    }

    for n in 2..=order {
        let next = ((2 * n - 1) as f64 * x * pn - (n - 1) as f64 * pnm1) / n as f64;
        pnm1 = pn;
        pn = next;
    }
    (pn, pnm1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrates_polynomials() {
        let quadrature = GaussLegendreUnit::new(16).unwrap();
        for degree in [0_i32, 1, 2, 5, 9] {
            let numerical = quadrature.integrate(|x| x.powi(degree));
            let exact = 1.0 / (degree as f64 + 1.0);
            assert!((numerical - exact).abs() < 2.0e-14);
        }
    }
}
