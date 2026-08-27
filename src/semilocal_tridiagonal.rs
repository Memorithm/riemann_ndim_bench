//! Direct tridiagonal eigensolver paths for the finite semilocal prolate blocks.
//!
//! `semilocal::crossing_derivatives` intentionally remains the dense reference
//! implementation. This module provides two independently testable paths:
//!
//! - a full tridiagonal EVD with eigenvectors for exact Rayleigh derivatives;
//! - an eigenvalues-only centered finite-difference path for high dimensions.
//!
//! The latter operates on the first-order linearized matrix
//! `K(0) + step * K'(0)` and forms spectral differences mode by mode before
//! summation, avoiding the cancellation of subtracting two large trace totals.

use std::fmt;

use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::linalg::evd::{
    ComputeEigenvectors, SelfAdjointEvdParams, self_adjoint_evd_scratch,
    tridiagonal_self_adjoint_evd,
};
use faer::{Auto, Col, Mat, MatRef, Par, Spec};

use crate::semilocal::{
    CrossingDerivative, ProlateParity, SemilocalError, SymmetricTridiagonal, build_k0,
    build_kprime_closed,
};

fn quadratic_form_column(
    matrix: &SymmetricTridiagonal,
    eigenvectors: MatRef<'_, f64>,
    column: usize,
) -> f64 {
    let mut value = 0.0;

    for i in 0..matrix.len() {
        let ui = eigenvectors[(i, column)];
        value += matrix.diagonal()[i] * ui * ui;
    }

    for i in 0..matrix.off_diagonal().len() {
        value += 2.0
            * matrix.off_diagonal()[i]
            * eigenvectors[(i, column)]
            * eigenvectors[(i + 1, column)];
    }

    value
}

fn padded_subdiagonal(
    block_size: usize,
    matrix: &SymmetricTridiagonal,
    derivative: Option<(&SymmetricTridiagonal, f64)>,
) -> Col<f64> {
    Col::from_fn(block_size, |i| {
        if i + 1 >= block_size {
            0.0
        } else if let Some((matrix_prime, step)) = derivative {
            matrix.off_diagonal()[i] + step * matrix_prime.off_diagonal()[i]
        } else {
            matrix.off_diagonal()[i]
        }
    })
}

/// Compute every finite crossing derivative by feeding `K(0)` directly to
/// faer's tridiagonal self-adjoint eigensolver.
///
/// The dense [`crate::semilocal::crossing_derivatives`] path is retained as an
/// independent reference. This function computes full eigenvectors because
/// the first derivative requires the Rayleigh quotient `u_j^T K'(0) u_j`.
pub fn crossing_derivatives_tridiagonal(
    block_size: usize,
    parity: ProlateParity,
) -> Result<Vec<CrossingDerivative>, SemilocalError> {
    if block_size == 0 {
        return Ok(Vec::new());
    }

    let k0 = build_k0(block_size, parity);
    let kprime = build_kprime_closed(block_size, parity);

    let diagonal = Col::from_fn(block_size, |i| k0.diagonal()[i]);
    let subdiagonal = padded_subdiagonal(block_size, &k0, None);
    let mut eigenvalues = Col::<f64>::zeros(block_size);
    let mut eigenvectors = Mat::<f64>::zeros(block_size, block_size);

    let params = Spec::new(<SelfAdjointEvdParams as Auto<f64>>::auto());
    // The full-eigenvector path intentionally uses faer's public dense scratch
    // helper as a conservative workspace bound. Eigenvector storage is O(m^2)
    // regardless; the separate eigenvalues-only path below uses O(m) workspace.
    let mut memory = MemBuffer::new(self_adjoint_evd_scratch::<f64>(
        block_size,
        ComputeEigenvectors::Yes,
        Par::Seq,
        params,
    ));

    tridiagonal_self_adjoint_evd(
        diagonal.as_diagonal(),
        subdiagonal.as_diagonal(),
        eigenvalues.as_diagonal_mut(),
        Some(eigenvectors.as_mut()),
        Par::Seq,
        MemStack::new(&mut memory),
        params,
    )?;

    let mut out = Vec::with_capacity(block_size);

    for j in 0..block_size {
        let mu = eigenvalues[j];
        if mu <= 0.0 || !mu.is_finite() {
            return Err(SemilocalError::NonPositiveEigenvalue {
                index: j,
                value: mu,
            });
        }

        let lambda = mu.sqrt();
        let mu_prime = quadratic_form_column(&kprime, eigenvectors.as_ref(), j);
        let lambda_prime = mu_prime / (2.0 * lambda);

        out.push(CrossingDerivative {
            parity,
            parity_index: j,
            mu,
            lambda,
            mu_prime,
            lambda_prime,
        });
    }

    Ok(out)
}

/// Return the sorted square-root spectrum of the first-order linearization
/// `K(0) + step K'(0)` using no eigenvectors.
///
/// For real `f64`, faer's `u=None` tridiagonal path needs only two temporary
/// real vectors. The explicit [`StackReq`] below therefore keeps the auxiliary
/// workspace O(m); output storage is O(m) as well.
///
/// `step` must be finite and small enough that the finite matrix remains
/// positive definite.
pub fn first_order_sqrt_spectrum_tridiagonal_eigenvalues_only(
    block_size: usize,
    parity: ProlateParity,
    step: f64,
) -> Result<Vec<f64>, SemilocalError> {
    assert!(step.is_finite(), "finite-difference step must be finite");

    if block_size == 0 {
        return Ok(Vec::new());
    }

    let k0 = build_k0(block_size, parity);
    let kprime = build_kprime_closed(block_size, parity);
    let diagonal = Col::from_fn(block_size, |i| {
        k0.diagonal()[i] + step * kprime.diagonal()[i]
    });
    let subdiagonal = padded_subdiagonal(block_size, &k0, Some((&kprime, step)));
    let mut eigenvalues = Col::<f64>::zeros(block_size);

    let params = Spec::new(<SelfAdjointEvdParams as Auto<f64>>::auto());
    // Inspection of faer 0.24.4's real-valued `u=None` implementation shows
    // that it allocates exactly two simultaneous n-by-1 real temporaries before
    // invoking the in-place tridiagonal QR algorithm. This public dyn-stack
    // requirement mirrors that storage without allocating an n-by-n buffer.
    let workspace = StackReq::new::<f64>(block_size).array(2);
    let mut memory = MemBuffer::new(workspace);

    tridiagonal_self_adjoint_evd(
        diagonal.as_diagonal(),
        subdiagonal.as_diagonal(),
        eigenvalues.as_diagonal_mut(),
        None,
        Par::Seq,
        MemStack::new(&mut memory),
        params,
    )?;

    let mut spectrum = Vec::with_capacity(block_size);
    for j in 0..block_size {
        let mu = eigenvalues[j];
        if mu <= 0.0 || !mu.is_finite() {
            return Err(SemilocalError::NonPositiveEigenvalue {
                index: j,
                value: mu,
            });
        }
        spectrum.push(mu.sqrt());
    }

    Ok(spectrum)
}

/// Error specific to centered pairwise finite differences.
#[derive(Debug)]
pub enum PairwiseError {
    /// Underlying semilocal eigensolver or positivity error.
    Semilocal(SemilocalError),
    /// Centered differences require a finite strictly positive step.
    InvalidStep { value: f64 },
}

impl fmt::Display for PairwiseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Semilocal(error) => write!(f, "{error}"),
            Self::InvalidStep { value } => {
                write!(f, "finite-difference step must be positive and finite: {value:e}")
            }
        }
    }
}

impl std::error::Error for PairwiseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Semilocal(error) => Some(error),
            Self::InvalidStep { .. } => None,
        }
    }
}

impl From<SemilocalError> for PairwiseError {
    fn from(value: SemilocalError) -> Self {
        Self::Semilocal(value)
    }
}

fn checked_pairwise_step(step: f64) -> Result<(), PairwiseError> {
    if step.is_finite() && step > 0.0 {
        Ok(())
    } else {
        Err(PairwiseError::InvalidStep { value: step })
    }
}

/// Centered mode-by-mode derivatives
/// `(lambda_j(+h) - lambda_j(-h)) / (2h)` for one parity sector.
pub fn pairwise_first_order_derivatives(
    block_size: usize,
    parity: ProlateParity,
    step: f64,
) -> Result<Vec<f64>, PairwiseError> {
    checked_pairwise_step(step)?;

    let forward = first_order_sqrt_spectrum_tridiagonal_eigenvalues_only(
        block_size, parity, step,
    )?;
    let backward = first_order_sqrt_spectrum_tridiagonal_eigenvalues_only(
        block_size, parity, -step,
    )?;

    Ok(forward
        .into_iter()
        .zip(backward)
        .map(|(plus, minus)| (plus - minus) / (2.0 * step))
        .collect())
}

/// Centered derivative of `Tr sqrt(K)` for one parity sector, obtained by
/// summing already-paired mode derivatives rather than subtracting two traces.
pub fn pairwise_first_order_trace_derivative(
    block_size: usize,
    parity: ProlateParity,
    step: f64,
) -> Result<f64, PairwiseError> {
    Ok(pairwise_first_order_derivatives(block_size, parity, step)?
        .into_iter()
        .sum())
}

/// Total first-order variation on `W+ union W-`, with absolute values taken
/// mode by mode before summation.
pub fn pairwise_total_variation_derivative(
    block_size: usize,
    step: f64,
) -> Result<f64, PairwiseError> {
    checked_pairwise_step(step)?;

    let mut total = 0.0;
    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        total += pairwise_first_order_derivatives(block_size, parity, step)?
            .into_iter()
            .map(f64::abs)
            .sum::<f64>();
    }
    Ok(total)
}

/// Three-level centered extrapolation diagnostics for a hierarchy
/// `h, h/2, h/4`, assuming the leading errors are `O(h^2)` and `O(h^4)`.
#[derive(Clone, Copy, Debug)]
pub struct PairwiseExtrapolation {
    pub d_h: f64,
    pub d_h2: f64,
    pub d_h4: f64,
    pub raw_difference_ratio: f64,
    pub richardson_h_h2: f64,
    pub richardson_h2_h4: f64,
    pub quadratic_h2_h4: f64,
}

/// Eliminate exact `h^2` and `h^4` terms from three values at `h, h/2, h/4`.
#[inline]
pub fn quadratic_even_power_extrapolate(d_h: f64, d_h2: f64, d_h4: f64) -> f64 {
    (d_h - 20.0 * d_h2 + 64.0 * d_h4) / 45.0
}

/// Compute the pairwise total variation at `h, h/2, h/4` and return the
/// order-2/order-4 diagnostics used by the Phase-4 high-dimensional checks.
pub fn pairwise_total_variation_extrapolation(
    block_size: usize,
    step: f64,
) -> Result<PairwiseExtrapolation, PairwiseError> {
    checked_pairwise_step(step)?;

    let d_h = pairwise_total_variation_derivative(block_size, step)?;
    let d_h2 = pairwise_total_variation_derivative(block_size, step / 2.0)?;
    let d_h4 = pairwise_total_variation_derivative(block_size, step / 4.0)?;
    let denominator = d_h2 - d_h4;
    let raw_difference_ratio = if denominator == 0.0 {
        f64::NAN
    } else {
        (d_h - d_h2) / denominator
    };
    let richardson_h_h2 = (4.0 * d_h2 - d_h) / 3.0;
    let richardson_h2_h4 = (4.0 * d_h4 - d_h2) / 3.0;
    let quadratic_h2_h4 = quadratic_even_power_extrapolate(d_h, d_h2, d_h4);

    Ok(PairwiseExtrapolation {
        d_h,
        d_h2,
        d_h4,
        raw_difference_ratio,
        richardson_h_h2,
        richardson_h2_h4,
        quadratic_h2_h4,
    })
}
