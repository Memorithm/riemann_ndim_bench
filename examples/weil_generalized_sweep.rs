use std::env;
use std::error::Error;

use riemann_ndim_bench::semilocal_compact_archimedean::{CompactArchimedeanBump, PositiveRational};
use riemann_ndim_bench::weil_generalized_spectrum::audit_finite_weil_generalized_spectrum;

fn parse_usize(index: usize, default: usize) -> Result<usize, Box<dyn Error>> {
    match env::args().nth(index) {
        Some(value) => Ok(value.parse()?),
        None => Ok(default),
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let max_dimension = parse_usize(1, 6)?;
    let correlation_order = parse_usize(2, 96)?;
    let archimedean_order = parse_usize(3, 96)?;
    let boundary_order = parse_usize(4, 128)?;
    let gram_order = parse_usize(5, 128)?;

    let bump =
        CompactArchimedeanBump::new(PositiveRational::new(1, 2)?, PositiveRational::new(7, 2)?)?;
    let audit = audit_finite_weil_generalized_spectrum(
        bump,
        max_dimension,
        correlation_order,
        archimedean_order,
        boundary_order,
        gram_order,
    )?;
    let rows = audit.principal_sweep()?;

    println!("support_lower=1/2");
    println!("support_upper=7/2");
    println!("max_dimension={max_dimension}");
    println!("correlation_order={correlation_order}");
    println!("archimedean_order={archimedean_order}");
    println!("boundary_order={boundary_order}");
    println!("gram_order={gram_order}");
    println!(
        "max_boundary_residual={:.15e}",
        audit.pairing().max_boundary_residual()
    );
    println!(
        "max_pairing_asymmetry={:.15e}",
        audit.pairing().max_raw_pairing_asymmetry()
    );
    println!(
        "max_whitened_asymmetry={:.15e}",
        audit.max_whitened_asymmetry()
    );
    println!(
        "dimension,raw_lambda_min,generalized_lambda_min,gram_lambda_min,gram_lambda_max,gram_condition"
    );
    for row in rows {
        println!(
            "{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
            row.dimension(),
            row.raw_minimum_eigenvalue(),
            row.generalized_minimum_eigenvalue(),
            row.gram_minimum_eigenvalue(),
            row.gram_maximum_eigenvalue(),
            row.gram_condition_number(),
        );
    }

    Ok(())
}
