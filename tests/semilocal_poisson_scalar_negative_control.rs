use std::f64::consts::PI;

use riemann_ndim_bench::semilocal_poisson::SemilocalPoissonMonoid;
use riemann_ndim_bench::semilocal_trace_contract::FinitePlaceSet;

const COEFFICIENTS: [f64; 3] = [1.0 / 3.0, -4.0 / 3.0, 1.0];
const GAUSSIAN_RATES: [f64; 3] = [1.0, 4.0, 9.0];

fn manufactured_schwartz(x: f64) -> f64 {
    COEFFICIENTS
        .into_iter()
        .zip(GAUSSIAN_RATES)
        .map(|(coefficient, rate)| coefficient * (-PI * rate * x * x).exp())
        .sum()
}

fn manufactured_fourier(x: f64) -> f64 {
    COEFFICIENTS
        .into_iter()
        .zip(GAUSSIAN_RATES)
        .map(|(coefficient, rate)| coefficient / rate.sqrt() * (-PI * x * x / rate).exp())
        .sum()
}

fn gaussian_tail_bound(scale: f64, transformed: bool, max_m: u64) -> f64 {
    let start = max_m as f64;
    COEFFICIENTS
        .into_iter()
        .zip(GAUSSIAN_RATES)
        .map(|(coefficient, rate)| {
            let amplitude = if transformed {
                coefficient.abs() / rate.sqrt()
            } else {
                coefficient.abs()
            };
            let exponent_rate = if transformed {
                PI * scale * scale / rate
            } else {
                PI * rate * scale * scale
            };
            amplitude * (-exponent_rate * start * start).exp() / (2.0 * exponent_rate * start)
        })
        .sum()
}

#[test]
fn manufactured_fixture_satisfies_the_archimedean_poisson_boundary_conditions() {
    assert!(manufactured_schwartz(0.0).abs() < 2.0e-16);
    assert!(manufactured_fourier(0.0).abs() < 2.0e-16);
}

#[test]
fn deleting_finite_place_multiples_on_the_real_axis_is_not_semilocal_poisson() {
    // Source equation (4.6) is an identity on the genuine semilocal quotient,
    // with p-adic coordinates and a self-dual basic character.  This test
    // deliberately evaluates the tempting but incorrect scalar surrogate:
    // keep only the real coordinate and replace the integer sum by M_S.
    let places = FinitePlaceSet::new(vec![2]).unwrap();
    let monoid = SemilocalPoissonMonoid::new(places);
    let x = 0.75_f64;
    let reciprocal = 1.0 / x;
    let max_m = 31_u64;

    let transformed_side = monoid
        .finite_e_sum(x, max_m, |m| manufactured_fourier(m as f64 * x))
        .unwrap();
    let reciprocal_side = monoid
        .finite_e_sum(reciprocal, max_m, |m| {
            manufactured_schwartz(m as f64 * reciprocal)
        })
        .unwrap();

    let residual = (transformed_side.value() - reciprocal_side.value()).abs();
    let omitted_bound = x.sqrt() * gaussian_tail_bound(x, true, max_m)
        + reciprocal.sqrt() * gaussian_tail_bound(reciprocal, false, max_m);
    let declared_roundoff = 2.0e-14;

    assert!(omitted_bound < 1.0e-40);
    assert!(
        residual > 1_000_000.0 * (omitted_bound + declared_roundoff),
        "naive scalar reduction unexpectedly looks Poisson-invariant: residual={residual:.16e}, tail={omitted_bound:.3e}"
    );
    assert!(residual > 4.0e-2);
}
