use riemann_ndim_bench::semilocal_compact_archimedean::{CompactArchimedeanBump, PositiveRational};
use riemann_ndim_bench::weil_sine_pairing_refinement::{
    SinePairingRefinementLevel, audit_finite_weil_sine_pairing_refinement,
};
use riemann_ndim_bench::weil_sine_truncation::SineModeSet;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bump =
        CompactArchimedeanBump::new(PositiveRational::new(1, 2)?, PositiveRational::new(7, 2)?)?;
    let modes = SineModeSet::new(vec![1, 2, 3])?;
    let parent_dimensions = [4_usize, 6, 8];
    let levels = [
        SinePairingRefinementLevel::new(48, 32, 32, 48),
        SinePairingRefinementLevel::new(72, 48, 48, 72),
        SinePairingRefinementLevel::new(96, 64, 64, 96),
    ];
    let grid = audit_finite_weil_sine_pairing_refinement(
        bump,
        &modes,
        &parent_dimensions,
        &levels,
    )?;

    println!(
        "parent_dimension,coefficient_order,correlation_order,archimedean_order,boundary_order,pairing_amplitude,direct_projected_residual,normalized_residual,parent_matrix_projection_residual,pole_residual,archimedean_residual,prime_residual,direct_asymmetry,projected_asymmetry,direct_boundary,projected_boundary"
    );
    for cell in grid.cells() {
        for sample in cell.samples().iter().copied() {
            let level = sample.level();
            println!(
                "{},{},{},{},{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
                cell.parent_dimension(),
                level.coefficient_quadrature_order(),
                level.correlation_order(),
                level.archimedean_order(),
                level.boundary_order(),
                sample.max_pairing_amplitude(),
                sample.max_direct_projected_pairing_residual(),
                sample.max_normalized_pairing_residual(),
                sample.max_parent_matrix_projection_residual(),
                sample.max_pole_term_residual(),
                sample.max_archimedean_term_residual(),
                sample.max_prime_total_residual(),
                sample.max_direct_pairing_asymmetry(),
                sample.max_projected_pairing_asymmetry(),
                sample.max_direct_boundary_residual(),
                sample.max_projected_boundary_residual(),
            );
        }
        eprintln!(
            "parent_dimension={}, pairing_residual_span={:.15e}, normalized_residual_span={:.15e}, last_pairing_delta={:?}, last_normalized_delta={:?}, max_matrix_projection_residual={:.15e}, max_asymmetry={:.15e}, max_boundary={:.15e}",
            cell.parent_dimension(),
            cell.pairing_residual_observed_span(),
            cell.normalized_residual_observed_span(),
            cell.last_pairing_residual_delta(),
            cell.last_normalized_residual_delta(),
            cell.max_parent_matrix_projection_residual(),
            cell.max_pairing_asymmetry(),
            cell.max_boundary_residual(),
        );
    }
    Ok(())
}
