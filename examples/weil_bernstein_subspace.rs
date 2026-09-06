use riemann_ndim_bench::semilocal_compact_archimedean::{CompactArchimedeanBump, PositiveRational};
use riemann_ndim_bench::weil_bernstein_subspace::{
    BernsteinIndexSubspace, audit_finite_weil_bernstein_subspace,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bump =
        CompactArchimedeanBump::new(PositiveRational::new(1, 2)?, PositiveRational::new(7, 2)?)?;
    let selected = BernsteinIndexSubspace::new(6, vec![0, 2, 4, 6])?;
    let audit = audit_finite_weil_bernstein_subspace(bump, &selected, 72, 72, 96, 96)?;

    println!(
        "family,degree,indices,dimension,raw_min,generalized_min,gram_condition,boundary_residual,whitened_asymmetry"
    );
    let indices = audit
        .indices()
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(":");
    println!(
        "bernstein,{},{indices},{},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
        audit.degree(),
        audit.dimension(),
        audit.minimum_raw_eigenvalue(),
        audit.minimum_generalized_eigenvalue(),
        audit.gram_condition_number(),
        audit.max_boundary_residual(),
        audit.max_whitened_asymmetry(),
    );
    println!(
        "leading_legendre,{},0:1:2:3,4,,{:.15e},{:.15e},,",
        audit.degree(),
        audit.leading_legendre_generalized_minimum(),
        audit.leading_legendre_gram_condition_number(),
    );
    eprintln!(
        "parent_dimension={}, parent_max_raw_pairing_asymmetry={:.15e}",
        audit.parent_dimension(),
        audit.parent_max_raw_pairing_asymmetry(),
    );
    Ok(())
}
