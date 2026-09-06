use riemann_ndim_bench::semilocal_compact_archimedean::PositiveRational;
use riemann_ndim_bench::weil_bernstein_evidence_grid::audit_finite_weil_bernstein_evidence_grid;
use riemann_ndim_bench::weil_bernstein_subspace::BernsteinIndexSubspace;
use riemann_ndim_bench::weil_refinement::WeilQuadratureLevel;
use riemann_ndim_bench::weil_support_sweep::WeilSupportWindow;

fn rational(
    numerator: u64,
    denominator: u64,
) -> Result<PositiveRational, Box<dyn std::error::Error>> {
    Ok(PositiveRational::new(numerator, denominator)?)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
    let subspace = BernsteinIndexSubspace::new(6, vec![0, 2, 4, 6])?;
    let grid = audit_finite_weil_bernstein_evidence_grid(&windows, &subspace, &levels)?;

    println!(
        "window,level,lower_num,lower_den,upper_num,upper_den,correlation_order,archimedean_order,boundary_order,gram_order,bernstein_raw_min,bernstein_generalized_min,leading_generalized_min,family_delta,bernstein_gram_condition,leading_gram_condition,boundary_residual,pairing_asymmetry,whitened_asymmetry"
    );
    for (window_index, cell) in grid.cells().iter().enumerate() {
        let window = cell.window();
        for (level_index, sample) in cell.samples().iter().copied().enumerate() {
            let level = sample.level();
            println!(
                "{window_index},{level_index},{},{},{},{},{},{},{},{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
                window.lower().numerator(),
                window.lower().denominator(),
                window.upper().numerator(),
                window.upper().denominator(),
                level.correlation_order(),
                level.archimedean_order(),
                level.boundary_order(),
                level.gram_order(),
                sample.bernstein_raw_minimum_eigenvalue(),
                sample.bernstein_generalized_minimum_eigenvalue(),
                sample.leading_legendre_generalized_minimum_eigenvalue(),
                sample.generalized_family_delta(),
                sample.bernstein_gram_condition_number(),
                sample.leading_legendre_gram_condition_number(),
                sample.max_boundary_residual(),
                sample.max_pairing_asymmetry(),
                sample.max_whitened_asymmetry(),
            );
        }
        eprintln!(
            "window={window_index}, bernstein_span={:.15e}, leading_span={:.15e}, family_delta_span={:.15e}, last_bernstein_delta={:?}, last_leading_delta={:?}, last_family_delta_change={:?}",
            cell.bernstein_observed_span(),
            cell.leading_observed_span(),
            cell.family_delta_observed_span(),
            cell.last_bernstein_delta(),
            cell.last_leading_delta(),
            cell.last_family_delta_change(),
        );
    }

    Ok(())
}
