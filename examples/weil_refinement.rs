use std::env;
use std::error::Error;

use riemann_ndim_bench::semilocal_compact_archimedean::PositiveRational;
use riemann_ndim_bench::weil_refinement::{
    WeilQuadratureLevel, audit_finite_weil_quadrature_refinement,
};
use riemann_ndim_bench::weil_support_sweep::WeilSupportWindow;

fn parse_usize(index: usize, default: usize) -> Result<usize, Box<dyn Error>> {
    match env::args().nth(index) {
        Some(value) => Ok(value.parse()?),
        None => Ok(default),
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let dimension = parse_usize(1, 4)?;
    let lower = PositiveRational::new(1, 2)?;
    let upper = PositiveRational::new(7, 2)?;
    let window = WeilSupportWindow::new(lower, upper);
    let levels = [
        WeilQuadratureLevel::new(48, 48, 64, 64),
        WeilQuadratureLevel::new(72, 72, 96, 96),
        WeilQuadratureLevel::new(96, 96, 128, 128),
    ];

    let audit = audit_finite_weil_quadrature_refinement(window, dimension, &levels)?;
    println!("dimension={dimension}");
    println!("support_lower=1/2");
    println!("support_upper=7/2");
    println!(
        "correlation_order,archimedean_order,boundary_order,gram_order,raw_lambda_min,generalized_lambda_min,gram_condition,max_boundary_residual,max_pairing_asymmetry,max_whitened_asymmetry"
    );
    for sample in audit.samples() {
        let level = sample.level();
        println!(
            "{},{},{},{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
            level.correlation_order(),
            level.archimedean_order(),
            level.boundary_order(),
            level.gram_order(),
            sample.raw_minimum_eigenvalue(),
            sample.generalized_minimum_eigenvalue(),
            sample.gram_condition_number(),
            sample.max_boundary_residual(),
            sample.max_pairing_asymmetry(),
            sample.max_whitened_asymmetry(),
        );
    }

    let (raw_min, raw_max) = audit.raw_observed_interval();
    let (generalized_min, generalized_max) = audit.generalized_observed_interval();
    println!("raw_observed_min={raw_min:.15e}");
    println!("raw_observed_max={raw_max:.15e}");
    println!("raw_observed_span={:.15e}", audit.raw_observed_span());
    println!("generalized_observed_min={generalized_min:.15e}");
    println!("generalized_observed_max={generalized_max:.15e}");
    println!(
        "generalized_observed_span={:.15e}",
        audit.generalized_observed_span()
    );
    if let Some(delta) = audit.last_raw_delta() {
        println!("last_raw_delta={delta:.15e}");
    }
    if let Some(delta) = audit.last_generalized_delta() {
        println!("last_generalized_delta={delta:.15e}");
    }

    Ok(())
}
