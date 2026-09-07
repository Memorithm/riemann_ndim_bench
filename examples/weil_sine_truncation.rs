use riemann_ndim_bench::semilocal_compact_archimedean::{CompactArchimedeanBump, PositiveRational};
use riemann_ndim_bench::weil_refinement::WeilQuadratureLevel;
use riemann_ndim_bench::weil_sine_truncation::{SineModeSet, audit_finite_weil_sine_truncation};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bump =
        CompactArchimedeanBump::new(PositiveRational::new(1, 2)?, PositiveRational::new(7, 2)?)?;
    let modes = SineModeSet::new(vec![1, 2, 3])?;
    let parent_dimensions = [4_usize, 6, 8, 10];
    let level = WeilQuadratureLevel::new(64, 64, 96, 96);
    let audit = audit_finite_weil_sine_truncation(bump, &modes, &parent_dimensions, 96, level)?;

    println!(
        "parent_dimension,reconstruction_residual,coefficient_l1,raw_min,generalized_min,leading_generalized_min,family_delta,gram_condition,leading_gram_condition,boundary_residual,pairing_asymmetry,whitened_asymmetry"
    );
    for sample in audit.samples().iter().copied() {
        println!(
            "{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
            sample.parent_dimension(),
            sample.max_reconstruction_residual(),
            sample.max_coefficient_l1_norm(),
            sample.raw_minimum_eigenvalue(),
            sample.generalized_minimum_eigenvalue(),
            sample.leading_legendre_generalized_minimum_eigenvalue(),
            sample.generalized_family_delta(),
            sample.gram_condition_number(),
            sample.leading_legendre_gram_condition_number(),
            sample.max_boundary_residual(),
            sample.max_pairing_asymmetry(),
            sample.max_whitened_asymmetry(),
        );
    }

    eprintln!(
        "generalized_span={:.15e}, family_delta_span={:.15e}, last_generalized_delta={:?}, last_reconstruction_delta={:?}, last_family_delta_change={:?}",
        audit.generalized_observed_span(),
        audit.family_delta_observed_span(),
        audit.last_generalized_delta(),
        audit.last_reconstruction_delta(),
        audit.last_family_delta_change(),
    );
    Ok(())
}
