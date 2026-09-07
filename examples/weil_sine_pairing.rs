use riemann_ndim_bench::semilocal_compact_archimedean::{
    CompactArchimedeanBump, PositiveRational,
};
use riemann_ndim_bench::weil_compact_pairing::CompactWeilPairingConfig;
use riemann_ndim_bench::weil_sine_pairing::{
    DirectSinePairingAuditConfig, audit_finite_weil_direct_sine_pairing,
};
use riemann_ndim_bench::weil_sine_truncation::SineModeSet;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bump =
        CompactArchimedeanBump::new(PositiveRational::new(1, 2)?, PositiveRational::new(7, 2)?)?;
    let modes = SineModeSet::new(vec![1, 2, 3])?;
    let parent_dimensions = [4_usize, 6, 8, 10];
    let config = DirectSinePairingAuditConfig::new(
        96,
        CompactWeilPairingConfig::new(64, 64, 96),
    );
    let audit =
        audit_finite_weil_direct_sine_pairing(bump, &modes, &parent_dimensions, config)?;

    println!(
        "parent_dimension,pairing_amplitude,direct_projected_residual,normalized_residual,parent_matrix_projection_residual,pole_residual,archimedean_residual,prime_residual,direct_asymmetry,projected_asymmetry,direct_boundary,projected_boundary"
    );
    for sample in audit.samples().iter().copied() {
        println!(
            "{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
            sample.parent_dimension(),
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
        "last_pairing_residual_delta={:?}",
        audit.last_pairing_residual_delta(),
    );
    Ok(())
}
