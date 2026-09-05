use std::f64::consts::PI;

use riemann_ndim_bench::semilocal_trace_contract::{
    BasicCharacterNormalization, FinitePlaceSet, SemilocalCutoff, SemilocalSpaceContract,
    SemilocalTraceContract, symmetric_test_value,
};

#[test]
fn source_contract_keeps_place_character_and_cutoff_data_together() {
    let places = FinitePlaceSet::new(vec![2, 3, 5]).unwrap();
    let space = SemilocalSpaceContract::qs_self_dual(places);
    let cutoff = SemilocalCutoff::new((2.0 * PI).exp()).unwrap();
    let contract = SemilocalTraceContract::new(space, cutoff);

    assert_eq!(contract.space().places().finite_primes(), &[2, 3, 5]);
    assert_eq!(contract.space().places().place_count(), 4);
    assert_eq!(
        contract.space().character_normalization(),
        BasicCharacterNormalization::QsSelfDual
    );
    assert!((contract.cutoff().quantized_band_endpoint() - 2.0).abs() < 4.0e-15);
}

#[test]
fn theorem_2_5_scalar_normalizations_are_source_locked() {
    let cutoff = SemilocalCutoff::new(8.0).unwrap();
    let expected = 2.0 * 1.25 * 8.0_f64.ln();
    assert_eq!(cutoff.theorem_2_5_leading_term(1.25), expected);
    assert_eq!(symmetric_test_value(9.0, 2.0).unwrap(), 6.0);
}
