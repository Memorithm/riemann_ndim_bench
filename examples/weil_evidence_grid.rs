use std::env;
use std::error::Error;

use riemann_ndim_bench::semilocal_compact_archimedean::PositiveRational;
use riemann_ndim_bench::weil_evidence_grid::audit_finite_weil_evidence_grid;
use riemann_ndim_bench::weil_refinement::WeilQuadratureLevel;
use riemann_ndim_bench::weil_support_sweep::WeilSupportWindow;

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
    let windows = [
        WeilSupportWindow::new(rational(3, 4)?, rational(13, 4)?),
        WeilSupportWindow::new(rational(1, 2)?, rational(7, 2)?),
        WeilSupportWindow::new(rational(1, 4)?, rational(15, 4)?),
    ];
    let levels = [
        WeilQuadratureLevel::new(48, 48, 64, 64),
        WeilQuadratureLevel::new(72, 72, 96, 96),
        WeilQuadratureLevel::new(96, 96, 128, 128),
    ];

    let grid = audit_finite_weil_evidence_grid(&windows, max_dimension, &levels)?;
    println!("max_dimension={max_dimension}");
    println!(
        "lower_num,lower_den,upper_num,upper_den,dimension,level_index,correlation_order,archimedean_order,boundary_order,gram_order,raw_lambda_min,generalized_lambda_min,gram_condition,max_boundary_residual,max_pairing_asymmetry,max_whitened_asymmetry,raw_observed_span,generalized_observed_span,last_raw_delta,last_generalized_delta,max_gram_condition"
    );

    for cell in grid.cells() {
        let lower = cell.window().lower();
        let upper = cell.window().upper();
        let last_raw_delta = cell.last_raw_delta().unwrap_or(0.0);
        let last_generalized_delta = cell.last_generalized_delta().unwrap_or(0.0);
        let max_gram_condition = cell.maximum_gram_condition_number();
        for (level_index, sample) in cell.samples().iter().enumerate() {
            let level = sample.level();
            println!(
                "{},{},{},{},{},{},{},{},{},{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
                lower.numerator(),
                lower.denominator(),
                upper.numerator(),
                upper.denominator(),
                cell.dimension(),
                level_index,
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
                cell.raw_observed_span(),
                cell.generalized_observed_span(),
                last_raw_delta,
                last_generalized_delta,
                max_gram_condition,
            );
        }
    }

    Ok(())
}
