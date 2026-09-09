use riemann_ndim_bench::semilocal_compact_archimedean::{
    CompactArchimedeanBump, PositiveRational,
};
use riemann_ndim_bench::weil_sine_taylor_h1::{
    bound_sine_taylor_pairing_perturbation, derive_sine_taylor_h1_envelope,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bump = CompactArchimedeanBump::new(
        PositiveRational::new(1, 2)?,
        PositiveRational::new(7, 2)?,
    )?;

    println!(
        "mode,degree,error_l2_bound,derivative_error_l2_bound,target_l2_bound,target_derivative_l2_bound"
    );
    for degree in [4_usize, 6, 8, 10, 12] {
        let audit = derive_sine_taylor_h1_envelope(bump, 1, degree)?;
        let h1 = audit.h1_bounds();
        println!(
            "{},{},{:.15e},{:.15e},{:.15e},{:.15e}",
            audit.mode(),
            audit.degree(),
            h1.error_l2(),
            h1.derivative_error_l2(),
            h1.target_l2(),
            h1.target_log_derivative_l2(),
        );
    }

    let pairing = bound_sine_taylor_pairing_perturbation(bump, 1, 10, 2, 12)?;
    println!(
        "pairing_bound_mode1_degree10_mode2_degree12={:.15e}",
        pairing.pairing().total_pairing_error_bound()
    );
    Ok(())
}
