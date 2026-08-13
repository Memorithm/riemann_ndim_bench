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
    ParameterNotFinite,
    IntervalNotPositive,
    StepNotPositive,
}

impl LogLattice {
    pub fn new(interval_length: f64, omega: f64) -> Result<Self, LatticeError> {
        if !interval_length.is_finite() || !omega.is_finite() {
            return Err(LatticeError::ParameterNotFinite);
        }
        if interval_length <= 0.0 {
            return Err(LatticeError::IntervalNotPositive);
        }
        if omega <= 0.0 {
            return Err(LatticeError::StepNotPositive);
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
    pub fn from_first_row(first_row: Vec<f64>) -> Result<Self, ToeplitzError> {
        if first_row.is_empty() {
            return Err(ToeplitzError::Empty);
        }
        if first_row.iter().any(|value| !value.is_finite()) {
            return Err(ToeplitzError::NonFiniteCoefficient);
        }
        Ok(Self { first_row })
    }

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

    pub fn apply(&self, input: &[f64]) -> Result<Vec<f64>, ToeplitzError> {
        if input.len() != self.dimension() {
            return Err(ToeplitzError::DimensionMismatch);
        }

        let mut output = vec![0.0; self.dimension()];
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

    pub fn dense(&self) -> Mat<f64> {
        Mat::from_fn(self.dimension(), self.dimension(), |i, j| {
            self.first_row[i.abs_diff(j)]
        })
    }

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
