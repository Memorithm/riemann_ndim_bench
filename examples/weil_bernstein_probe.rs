use riemann_ndim_bench::semilocal_compact_archimedean::{CompactArchimedeanBump, PositiveRational};
use riemann_ndim_bench::weil_bernstein_probe::audit_finite_weil_bernstein_probes;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bump =
        CompactArchimedeanBump::new(PositiveRational::new(1, 2)?, PositiveRational::new(7, 2)?)?;
    let audit = audit_finite_weil_bernstein_probes(
        bump,
        6,
        &[0, 2, 4, 6],
        72,
        72,
        96,
        96,
    )?;

    println!(
        "degree,index,raw_quadratic,gram_norm_squared,generalized_rayleigh,boundary_residual,reconstruction_residual,coefficient_l1"
    );
    for probe in audit.probes() {
        println!(
            "{},{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
            audit.degree(),
            probe.index(),
            probe.raw_quadratic_value(),
            probe.gram_norm_squared(),
            probe.generalized_rayleigh_quotient(),
            probe.max_boundary_residual(),
            probe.max_reconstruction_residual(),
            probe.coefficient_l1_norm(),
        );
    }

    eprintln!(
        "parent_dimension={}, parent_gram_condition={:.15e}, parent_pairing_asymmetry={:.15e}, parent_whitened_asymmetry={:.15e}",
        audit.parent_dimension(),
        audit.parent_gram_condition_number(),
        audit.parent_max_raw_pairing_asymmetry(),
        audit.parent_max_whitened_asymmetry(),
    );
    Ok(())
}
