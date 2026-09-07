use riemann_ndim_bench::semilocal_compact_archimedean::{CompactArchimedeanBump, PositiveRational};
use riemann_ndim_bench::weil_sine_direct_q::{
    DirectSineQAuditConfig, audit_finite_weil_direct_sine_q,
};
use riemann_ndim_bench::weil_sine_truncation::SineModeSet;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bump =
        CompactArchimedeanBump::new(PositiveRational::new(1, 2)?, PositiveRational::new(7, 2)?)?;
    let modes = SineModeSet::new(vec![1, 2, 3])?;
    let parent_dimensions = [4_usize, 6, 8, 10];
    let config = DirectSineQAuditConfig::new(96, 128, 512);
    let audit = audit_finite_weil_direct_sine_q(bump, &modes, &parent_dimensions, config)?;

    println!(
        "parent_dimension,generator_residual,direct_q_amplitude,projected_q_amplitude,q_abs_residual,q_normalized_residual,direct_boundary_residual,projected_boundary_residual"
    );
    for sample in audit.samples().iter().copied() {
        println!(
            "{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
            sample.parent_dimension(),
            sample.max_generator_residual(),
            sample.max_direct_q_amplitude(),
            sample.max_projected_q_amplitude(),
            sample.max_q_absolute_residual(),
            sample.normalized_q_residual(),
            sample.max_direct_boundary_residual(),
            sample.max_projected_boundary_residual(),
        );
    }

    eprintln!(
        "last_q_residual_delta={:?}, last_normalized_q_residual_delta={:?}",
        audit.last_q_residual_delta(),
        audit.last_normalized_q_residual_delta(),
    );
    Ok(())
}
