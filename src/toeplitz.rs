use faer::{Mat, Side, linalg::solvers::SelfAdjointEigen};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LogLattice {
    interval_length: f64,
    omega: f64,
    q: f64,
    max_index: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LatticeError {
    NonFiniteParameter,
    NonPositiveInterval,
    NonPositiveStep,
}

impl LogLattice {
    /// Discretization used in Connes--Consani §6.1:
    /// omega = log(q), j omega in [0,a], N = floor(a / omega).
    pub fn new(interval_length: f64, omega: f64) -> Result<Self, LatticeError> {
        if !interval_length.is_finite() || !omega.is_finite() {
            return Err(LatticeError::NonFiniteParameter);
        }
        if interval_length <= 0.0 {
            return Err(LatticeError::NonPositiveInterval);
        }
        if omega <= 0.0 {
            return Err(LatticeError::NonPositiveStep);
        }

        let max_index = (interval_length / omega).floor() as usize;
        Ok(Self {
            interval_length,
            omega,
            q: omega.exp(),
            max_index,
        })
    }

    pub fn interval_length(self) -> f64 {
        self.interval_length
    }

    pub fn omega(self) -> f64 {
        self.omega
    }

    pub fn q(self) -> f64 {
        self.q
    }

    pub fn max_index(self) -> usize {
        self.max_index
    }

    /// Number of lattice points j = 0,...,N.
    pub fn dimension(self) -> usize {
        self.max_index + 1
    }

    pub fn x(self, index: usize) -> Option<f64> {
        (index <= self.max_index).then_some(index as f64 * self.omega)
    }

    pub fn rho(self, index: usize) -> Option<f64> {
        self.x(index).map(f64::exp)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SymmetricToeplitz {
    first_row: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToeplitzError {
    Empty,
    NonFiniteCoefficient,
    DimensionMismatch,
    DecompositionFailed,
}

impl SymmetricToeplitz {
    /// Builds a real self-adjoint Toeplitz operator from its first row.
    /// No positivity is imposed or assumed.
    pub fn from_first_row(first_row: Vec<f64>) -> Result<Self, ToeplitzError> {
        if first_row.is_empty() {
            return Err(ToeplitzError::Empty);
        }
        if first_row.iter().any(|value| !value.is_finite()) {
            return Err(ToeplitzError::NonFiniteCoefficient);
        }
        Ok(Self { first_row })
    }

    /// Builds the normalized finite operator appearing in Connes--Consani
    /// equations (105)--(106): T_q has entries omega * chi(|i-j| omega),
    /// where chi(x) = (Q epsilon)(exp(x)) / (2 epsilon'(1+)).
    ///
    /// This function only performs the discretization. The caller must provide
    /// `chi`; the research bench does not synthesize or infer that kernel.
    pub fn sample_normalized_kernel(
        lattice: LogLattice,
        mut chi: impl FnMut(f64) -> f64,
    ) -> Result<Self, ToeplitzError> {
        let mut first_row = Vec::with_capacity(lattice.dimension());
        for k in 0..lattice.dimension() {
            let x = k as f64 * lattice.omega();
            first_row.push(lattice.omega() * chi(x));
        }
        Self::from_first_row(first_row)
    }

    pub fn dimension(&self) -> usize {
        self.first_row.len()
    }

    pub fn first_row(&self) -> &[f64] {
        &self.first_row
    }

    pub fn entry(&self, row: usize, col: usize) -> Option<f64> {
        if row >= self.dimension() || col >= self.dimension() {
            return None;
        }
        Some(self.first_row[row.abs_diff(col)])
    }

    /// Matrix-free O(N^2) application. It stores only N coefficients rather
    /// than an N x N matrix.
    pub fn apply(&self, input: &[f64]) -> Result<Vec<f64>, ToeplitzError> {
        if input.len() != self.dimension() {
            return Err(ToeplitzError::DimensionMismatch);
        }

        let n = self.dimension();
        let mut output = vec![0.0; n];
        for (i, out) in output.iter_mut().enumerate() {
            let mut sum = 0.0;
            for (j, &value) in input.iter().enumerate() {
                sum += self.first_row[i.abs_diff(j)] * value;
            }
            *out = sum;
        }
        Ok(output)
    }

    pub fn quadratic_form(&self, vector: &[f64]) -> Result<f64, ToeplitzError> {
        let image = self.apply(vector)?;
        Ok(vector
            .iter()
            .zip(image.iter())
            .map(|(&left, &right)| left * right)
            .sum())
    }

    /// Dense materialization is deliberately isolated for validation and exact
    /// finite-dimensional eigendecomposition. Production kernels can continue
    /// to use the matrix-free `apply` path.
    pub fn dense(&self) -> Mat<f64> {
        Mat::from_fn(self.dimension(), self.dimension(), |i, j| {
            self.first_row[i.abs_diff(j)]
        })
    }

    /// Returns the real eigenvalues in nondecreasing order.
    pub fn eigenvalues(&self) -> Result<Vec<f64>, ToeplitzError> {
        let matrix = self.dense();
        let decomposition = SelfAdjointEigen::new(matrix.as_ref(), Side::Lower)
            .map_err(|_| ToeplitzError::DecompositionFailed)?;
        let diagonal = decomposition.S().column_vector();
        let mut values = (0..self.dimension())
            .map(|index| diagonal[index])
            .collect::<Vec<_>>();
        values.sort_by(f64::total_cmp);
        Ok(values)
    }

    pub fn largest_eigenvalue(&self) -> Result<f64, ToeplitzError> {
        self.eigenvalues()?
            .last()
            .copied()
            .ok_or(ToeplitzError::Empty)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let lattice = LogLattice::new(0.3, 0.1).unwrap();
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
        let expected = [3.5, 4.0, 1.5];
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
        let lattice = LogLattice::new(0.3, 0.1).unwrap();
        let matrix = SymmetricToeplitz::sample_normalized_kernel(lattice, |_| 1.0).unwrap();
        let eigenvalues = matrix.eigenvalues().unwrap();
        let largest = *eigenvalues.last().unwrap();
        assert!((largest - lattice.dimension() as f64 * lattice.omega()).abs() < EPS);
        for &eigenvalue in &eigenvalues[..eigenvalues.len() - 1] {
            assert!(eigenvalue.abs() < 10.0 * EPS);
        }
    }
}
