use riemann_ndim_bench::semilocal_compact_archimedean::{CompactArchimedeanBump, PositiveRational};
use riemann_ndim_bench::weil_compact_pairing::CompactWeilPairingConfig;
use riemann_ndim_bench::weil_sine_continuity::{
    SineContinuityProbeConfig, audit_finite_weil_sine_continuity,
};
use riemann_ndim_bench::weil_sine_truncation::SineModeSet;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bump =
        CompactArchimedeanBump::new(PositiveRational::new(1, 2)?, PositiveRational::new(7, 2)?)?;
    let modes = SineModeSet::new(vec![1, 2, 3])?;
    let parent_dimensions = [4_usize, 6, 8, 10];
    let config = SineContinuityProbeConfig::new(
        96,
        96,
        CompactWeilPairingConfig::new(64, 64, 96),
    );
    let audit = audit_finite_weil_sine_continuity(bump, &modes, &parent_dimensions, config)?;

    println!(
        "parent_dimension,max_l2_error,max_relative_l2_error,max_pairing_residual,max_observed_l2_continuity_quotient,max_pole_residual,max_archimedean_residual,max_prime_residual,max_pairing_asymmetry"
    );
    for sample in audit.samples() {
        println!(
            "{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
            sample.parent_dimension(),
            sample.max_l2_error(),
            sample.max_relative_l2_error(),
            sample.max_pairing_residual(),
            sample.max_observed_l2_continuity_quotient(),
            sample.max_pole_residual(),
            sample.max_archimedean_residual(),
            sample.max_prime_residual(),
            sample.max_pairing_asymmetry(),
        );
    }

    eprintln!("last_l2_error_delta={:?}", audit.last_l2_error_delta());
    eprintln!(
        "last_pairing_residual_delta={:?}",
        audit.last_pairing_residual_delta()
    );
    eprintln!(
        "last_continuity_quotient_delta={:?}",
        audit.last_continuity_quotient_delta()
    );
    Ok(())
}
