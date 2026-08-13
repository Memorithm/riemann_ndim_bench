mod numerics;
mod toeplitz;

use num_complex::Complex64;
use numerics::{xi_symmetry_residual, zeta_checked};
use riemann_ndim_bench::{PiRadialGeometry, SpectralPoint};
use std::f64::consts::LN_2;
use toeplitz::{LogLattice, SymmetricToeplitz};

fn main() {
    println!("Riemann N-dimensional bench - phase 2");
    println!("Geometry, independent zeta/xi checks, and sourced Toeplitz discretization.\n");

    println!("sigma      mirror      raw_radius       normalized       log_normalized");
    println!("--------------------------------------------------------------------------");

    for sigma in [0.1_f64, 0.3, 0.49, 0.5, 0.51, 0.7, 0.9] {
        let point = SpectralPoint::new(sigma, 0.0);
        let mirror = point.critical_line_reflection();
        let raw = PiRadialGeometry::raw_radius(sigma);
        let normalized = PiRadialGeometry::normalized_radius(sigma);
        let log_normalized = PiRadialGeometry::log_normalized_radius(sigma);

        println!(
            "{sigma:<10.3} {mirror_sigma:<10.3} {raw:<16.10} {normalized:<16.10} {log_normalized:<+16.10}",
            mirror_sigma = mirror.sigma,
        );
    }

    println!(
        "\ncritical raw radius = {:.15}",
        PiRadialGeometry::critical_radius()
    );

    println!("\nIndependent numerical checks:");
    for s in [
        Complex64::new(2.0, 0.0),
        Complex64::new(0.5, 14.0),
        Complex64::new(0.37, 9.25),
    ] {
        let estimate = zeta_checked(s).expect("zeta evaluation failed");
        println!(
            "s={:.3}{:+.3}i  zeta={:.12}{:+.12}i  delta={:.3e}  N={}  B={}",
            s.re,
            s.im,
            estimate.value.re,
            estimate.value.im,
            estimate.cross_resolution_delta,
            estimate.fine_n,
            estimate.bernoulli_terms,
        );
    }

    let probe = Complex64::new(0.37, 9.25);
    let residual = xi_symmetry_residual(probe).expect("xi evaluation failed");
    println!("xi(s)=xi(1-s) relative residual at probe: {residual:.3e}");

    println!("\nN-dimensional Toeplitz infrastructure:");
    let lattice = LogLattice::new(LN_2, 1.0e-3).expect("valid logarithmic lattice");
    println!(
        "interval=ln(2), omega={:.3e}, q={:.12}, max_index={}, dimension={}",
        lattice.omega(),
        lattice.q(),
        lattice.max_index(),
        lattice.dimension()
    );

    let synthetic = SymmetricToeplitz::from_first_row(vec![2.0, 1.0])
        .expect("valid synthetic Toeplitz operator");
    let eigenvalues = synthetic
        .eigenvalues()
        .expect("self-adjoint eigendecomposition failed");
    println!("synthetic 2x2 spectrum = {eigenvalues:?}");
    println!("Q-epsilon kernel is intentionally not synthesized in phase 2.");
}
