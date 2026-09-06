use std::env;
use std::error::Error;

use riemann_ndim_bench::semilocal_compact_archimedean::PositiveRational;
use riemann_ndim_bench::weil_support_sweep::{WeilSupportWindow, audit_finite_weil_support_sweep};

fn parse_usize(index: usize, default: usize) -> Result<usize, Box<dyn Error>> {
    match env::args().nth(index) {
        Some(value) => Ok(value.parse()?),
        None => Ok(default),
    }
}

fn rational(numerator: u64, denominator: u64) -> Result<PositiveRational, Box<dyn Error>> {
    Ok(PositiveRational::new(numerator, denominator)?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let max_dimension = parse_usize(1, 4)?;
    let correlation_order = parse_usize(2, 96)?;
    let archimedean_order = parse_usize(3, 96)?;
    let boundary_order = parse_usize(4, 128)?;
    let gram_order = parse_usize(5, 128)?;

    // Manufactured width sweep around the same arithmetic midpoint rho=2.
    // These windows are numerical audit parameters, not privileged supports.
    let windows = [
        WeilSupportWindow::new(rational(3, 4)?, rational(13, 4)?),
        WeilSupportWindow::new(rational(1, 2)?, rational(7, 2)?),
        WeilSupportWindow::new(rational(1, 4)?, rational(15, 4)?),
    ];

    let audits = audit_finite_weil_support_sweep(
        &windows,
        max_dimension,
        correlation_order,
        archimedean_order,
        boundary_order,
        gram_order,
    )?;

    println!("max_dimension={max_dimension}");
    println!("correlation_order={correlation_order}");
    println!("archimedean_order={archimedean_order}");
    println!("boundary_order={boundary_order}");
    println!("gram_order={gram_order}");
    println!(
        "lower_num,lower_den,upper_num,upper_den,dimension,raw_lambda_min,generalized_lambda_min,gram_lambda_min,gram_lambda_max,gram_condition,max_boundary_residual,max_pairing_asymmetry,max_whitened_asymmetry"
    );

    for audit in audits {
        let lower = audit.window().lower();
        let upper = audit.window().upper();
        for row in audit.rows() {
            println!(
                "{},{},{},{},{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
                lower.numerator(),
                lower.denominator(),
                upper.numerator(),
                upper.denominator(),
                row.dimension(),
                row.raw_minimum_eigenvalue(),
                row.generalized_minimum_eigenvalue(),
                row.gram_minimum_eigenvalue(),
                row.gram_maximum_eigenvalue(),
                row.gram_condition_number(),
                audit.max_boundary_residual(),
                audit.max_pairing_asymmetry(),
                audit.max_whitened_asymmetry(),
            );
        }
    }

    Ok(())
}
