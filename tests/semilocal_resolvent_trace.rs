use std::f64::consts::PI;

use riemann_ndim_bench::semilocal::{ProlateParity, build_kprime_closed};
use riemann_ndim_bench::semilocal_resolvent::{
    ResolventTraceError, SignCorrectedResolventTraceKernel,
    sign_corrected_resolvent_trace_tridiagonal,
};
use riemann_ndim_bench::semilocal_tridiagonal::crossing_derivatives_tridiagonal;

fn assert_close_scaled(actual: f64, expected: f64, relative: f64, absolute: f64, label: &str) {
    let error = (actual - expected).abs();
    let scale = actual.abs().max(expected.abs());
    let tolerance = absolute.max(relative * scale);
    assert!(
        error <= tolerance,
        "{label}: actual={actual:.16e} expected={expected:.16e} error={error:.3e} tolerance={tolerance:.3e}"
    );
}

fn square_root_trace_from_resolvent(
    block_size: usize,
    parity: ProlateParity,
    panels: usize,
) -> f64 {
    assert!(panels > 0 && panels.is_multiple_of(2));

    let kernel = SignCorrectedResolventTraceKernel::new(block_size, parity);
    let sign = parity.sign_correction();
    let h_trace = sign
        * build_kprime_closed(block_size, parity)
            .diagonal()
            .iter()
            .sum::<f64>();

    // With x = tan(theta),
    // 1/2 Tr(K^-1/2 H) = (1/pi) int_0^infinity Tr[(K+x^2 I)^-1 H] dx
    // becomes a smooth finite-interval integral. At theta=pi/2 the integrand
    // tends to Tr(H)/pi because the resolvent is asymptotic to x^-2 I.
    let step = (0.5 * PI) / panels as f64;
    let mut weighted_sum = 0.0;

    for panel in 0..=panels {
        let integrand = if panel == panels {
            h_trace / PI
        } else {
            let theta = panel as f64 * step;
            let x = theta.tan();
            kernel.trace(x * x).unwrap() * (1.0 + x * x) / PI
        };
        let weight = if panel == 0 || panel == panels {
            1.0
        } else if panel % 2 == 0 {
            2.0
        } else {
            4.0
        };
        weighted_sum += weight * integrand;
    }

    weighted_sum * step / 3.0
}

#[test]
fn selected_inverse_resolvent_trace_matches_spectral_oracle() {
    for block_size in [1_usize, 2, 4, 8, 32, 128] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let spectral = crossing_derivatives_tridiagonal(block_size, parity).unwrap();
            let kernel = SignCorrectedResolventTraceKernel::new(block_size, parity);
            let sign = parity.sign_correction();

            for shift in [0.0_f64, 1.0e-9, 1.0e-6, 1.0e-3, 1.0, 100.0] {
                let expected: f64 = spectral
                    .iter()
                    .map(|crossing| sign * crossing.mu_prime / (crossing.mu + shift))
                    .sum();
                let actual = kernel.trace(shift).unwrap();

                assert_close_scaled(
                    actual,
                    expected,
                    5.0e-11,
                    1.0e-11,
                    &format!("m={block_size} parity={parity:?} shift={shift:e}"),
                );
            }
        }
    }
}

#[test]
fn resolvent_integral_reproduces_exact_square_root_trace_derivative() {
    for block_size in [8_usize, 32, 128] {
        for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
            let sign = parity.sign_correction();
            let exact: f64 = 0.5
                * crossing_derivatives_tridiagonal(block_size, parity)
                    .unwrap()
                    .into_iter()
                    .map(|crossing| sign * crossing.mu_prime / crossing.lambda)
                    .sum::<f64>();
            let integrated = square_root_trace_from_resolvent(block_size, parity, 512);

            assert_close_scaled(
                integrated,
                exact,
                2.0e-10,
                2.0e-10,
                &format!("m={block_size} parity={parity:?} square-root resolvent integral"),
            );
        }
    }
}

#[test]
fn sign_corrected_resolvent_trace_is_positive_and_decreases_with_shift() {
    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        let kernel = SignCorrectedResolventTraceKernel::new(64, parity);
        let mut previous = f64::INFINITY;
        for shift in [0.0_f64, 1.0e-6, 1.0e-3, 1.0, 100.0] {
            let trace = kernel.trace(shift).unwrap();
            assert!(trace.is_finite() && trace > 0.0);
            assert!(
                trace < previous,
                "parity={parity:?} shift={shift:e}: trace={trace:e} previous={previous:e}"
            );
            previous = trace;
        }
    }
}

#[test]
fn empty_resolvent_block_has_zero_trace() {
    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        assert_eq!(
            sign_corrected_resolvent_trace_tridiagonal(0, parity, 0.0).unwrap(),
            0.0
        );
    }
}

#[test]
fn invalid_resolvent_shifts_are_rejected() {
    let kernel = SignCorrectedResolventTraceKernel::new(8, ProlateParity::WPlus);
    for shift in [-1.0e-6_f64, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            kernel.trace(shift),
            Err(ResolventTraceError::InvalidShift { .. })
        ));
    }
}
