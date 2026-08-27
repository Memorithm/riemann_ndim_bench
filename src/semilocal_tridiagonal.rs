//! Direct tridiagonal eigensolver path for the finite semilocal prolate blocks.
//!
//! `semilocal::crossing_derivatives` intentionally remains the dense reference
//! implementation.  This module feeds the already-tridiagonal `K(0)` directly
//! to faer 0.24's low-level self-adjoint tridiagonal EVD and computes the same
//! Rayleigh derivatives against the exact `K'(0)` tridiagonal.  Keeping both
//! paths available makes the optimization independently regression-testable.

use faer::dyn_stack::{MemBuffer, MemStack};
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

/// Compute every finite crossing derivative by feeding `K(0)` directly to
/// faer's tridiagonal self-adjoint eigensolver.
///
/// The dense [`crate::semilocal::crossing_derivatives`] path is retained as an
/// independent reference.  This function computes full eigenvectors because
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
    // faer's low-level routine expects a diagonal view with the same dimension
    // as the main diagonal and reads entries 0..n-2.  The final padded entry is
    // therefore explicitly zero.
    let subdiagonal = Col::from_fn(block_size, |i| {
        if i + 1 < block_size {
            k0.off_diagonal()[i]
        } else {
            0.0
        }
    });
    let mut eigenvalues = Col::<f64>::zeros(block_size);
    let mut eigenvectors = Mat::<f64>::zeros(block_size, block_size);

    let params = Spec::new(<SelfAdjointEvdParams as Auto<f64>>::auto());
    // faer does not currently expose a dedicated public scratch-size helper
    // for `tridiagonal_self_adjoint_evd`.  The documented dense self-adjoint
    // scratch requirement is a conservative superset and keeps this path on a
    // public, versioned API.  A later eigenvalues-only path can use a tighter
    // workspace because it does not need eigenvectors.
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
