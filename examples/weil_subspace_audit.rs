use riemann_ndim_bench::semilocal_compact_archimedean::{
    CompactArchimedeanBump, PositiveRational,
};
use riemann_ndim_bench::weil_subspace::{
    LegendreDegreeSubspace, audit_finite_weil_legendre_subspaces,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bump = CompactArchimedeanBump::new(
        PositiveRational::new(1, 2)?,
        PositiveRational::new(7, 2)?,
    )?;

    let subspaces = [
        LegendreDegreeSubspace::new(vec![0, 1, 2, 3])?,
        LegendreDegreeSubspace::new(vec![0, 2, 4, 6])?,
    ];
    let audit = audit_finite_weil_legendre_subspaces(
        bump,
        &subspaces,
        72,
        72,
        96,
        96,
    )?;

    println!(
        "subspace,degrees,dimension,raw_min,generalized_min,gram_condition,boundary_residual,whitened_asymmetry"
    );
    for (index, spectrum) in audit.spectra().iter().enumerate() {
        let degrees = spectrum
            .degrees()
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(":");
        println!(
            "{index},{degrees},{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
            spectrum.dimension(),
            spectrum.minimum_raw_eigenvalue(),
            spectrum.minimum_generalized_eigenvalue(),
            spectrum.gram_condition_number(),
            spectrum.max_boundary_residual(),
            spectrum.max_whitened_asymmetry(),
        );
    }

    eprintln!(
        "parent_dimension={}, parent_max_raw_pairing_asymmetry={:.15e}",
        audit.computed_legendre_dimension(),
        audit.full_matrix_max_raw_pairing_asymmetry(),
    );
    Ok(())
}
