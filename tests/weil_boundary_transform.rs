use std::f64::consts::LN_2;

use riemann_ndim_bench::weil_boundary::{
    MultiplicativeSupport, critical_boundary_moments, q_on_support,
};

fn smooth_log_bump_value_and_second(rho: f64) -> (f64, f64) {
    let x = rho.ln();
    let u = x / LN_2;
    if u.abs() >= 1.0 {
        return (0.0, 0.0);
    }

    let v = 1.0 - u * u;
    let value = (-1.0 / v).exp();
    let first_log_exponent_derivative = -2.0 * u / (v * v);
    let second_log_exponent_derivative = -2.0 / (v * v) - 8.0 * u * u / (v * v * v);
    let second_u_derivative = (second_log_exponent_derivative
        + first_log_exponent_derivative * first_log_exponent_derivative)
        * value;
    let second_log_derivative = second_u_derivative / (LN_2 * LN_2);
    (value, second_log_derivative)
}

#[test]
fn q_preserves_declared_compact_support_exactly() {
    let support = MultiplicativeSupport::new(0.5, 2.0).unwrap();

    for rho in [0.125_f64, 0.49, 2.01, 8.0] {
        let q_value = q_on_support(support, rho, smooth_log_bump_value_and_second).unwrap();
        assert_eq!(q_value, 0.0, "rho={rho}");
    }
}

#[test]
fn q_image_satisfies_the_two_weil_boundary_moments() {
    let support = MultiplicativeSupport::new(0.5, 2.0).unwrap();
    let moments = critical_boundary_moments(support, 128, |rho| {
        q_on_support(support, rho, smooth_log_bump_value_and_second).unwrap()
    })
    .unwrap();

    assert!(
        moments.plus_half.abs() < 2.0e-11,
        "plus_half={:.16e}",
        moments.plus_half
    );
    assert!(
        moments.minus_half.abs() < 2.0e-11,
        "minus_half={:.16e}",
        moments.minus_half
    );
}

#[test]
fn symmetric_support_matches_the_source_factor_two_window() {
    let support = MultiplicativeSupport::new(2.0_f64.powf(-0.5), 2.0_f64.powf(0.5)).unwrap();
    let width = support.log_upper() - support.log_lower();
    assert!((width - LN_2).abs() < 8.0 * f64::EPSILON);
}
