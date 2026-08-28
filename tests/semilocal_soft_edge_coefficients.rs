use std::f64::consts::PI;

use riemann_ndim_bench::semilocal::{ProlateParity, build_k0, build_kprime_closed};

#[derive(Clone, Copy, Debug)]
struct SoftEdgeRow {
    b_bar: f64,
    potential: f64,
    weight: f64,
    trace_prefactor: f64,
    infrared_scale: f64,
}

fn row_diagnostics(row: usize, parity: ProlateParity) -> SoftEdgeRow {
    assert!(row > 0);

    let block_size = row + 2;
    let k0 = build_k0(block_size, parity);
    let kprime = build_kprime_closed(block_size, parity);

    let backward = k0.off_diagonal()[row - 1];
    let forward = k0.off_diagonal()[row];
    let b_bar = 0.5 * (backward + forward);
    let potential = k0.diagonal()[row] - backward - forward;

    let sign = parity.sign_correction();
    let h_diagonal = sign * kprime.diagonal()[row];
    let h_backward = sign * kprime.off_diagonal()[row - 1];
    let h_forward = sign * kprime.off_diagonal()[row];
    let weight = h_diagonal - h_backward - h_forward;

    assert!(b_bar > 0.0);
    assert!(potential > 0.0);
    assert!(weight > 0.0);

    let row_f = row as f64;
    let trace_prefactor = PI * row_f * weight / b_bar.sqrt();
    let infrared_scale = 4.0 * row_f * (potential / b_bar).sqrt();

    SoftEdgeRow {
        b_bar,
        potential,
        weight,
        trace_prefactor,
        infrared_scale,
    }
}

fn parity_constant(parity: ProlateParity) -> f64 {
    match parity {
        ProlateParity::WPlus => 1.0,
        ProlateParity::WMinus => 5.0,
    }
}

fn assert_near(actual: f64, expected: f64, tolerance: f64, label: &str) {
    let error = (actual - expected).abs();
    assert!(
        error <= tolerance,
        "{label}: actual={actual:.16e} expected={expected:.16e} error={error:.3e} tolerance={tolerance:.3e}"
    );
}

#[test]
fn validates_documented_soft_edge_row_asymptotics() {
    // The Phase-4 notes deliberately evaluate V_i at i=4096 because it is a
    // cancellation of three O(i) coefficients and begins to lose f64 digits
    // at substantially larger rows. The other quantities remain stable at
    // i=16384 and are checked there.
    let potential_row = 4096_usize;
    let stable_row = 16384_usize;

    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        let a = parity_constant(parity);

        let potential_sample = row_diagnostics(potential_row, parity);
        let i = potential_row as f64;
        let scaled_potential = 8.0 * i * (64.0 * PI * i * potential_sample.potential - 1.0);
        assert_near(
            scaled_potential,
            -a,
            3.0e-3,
            &format!("{parity:?} potential first correction"),
        );

        let sample = row_diagnostics(stable_row, parity);
        let i = stable_row as f64;

        // b_bar(i) = i/(4 pi) * [1 + a/(8 i) + O(i^-2)].
        let scaled_b = 8.0 * i * (4.0 * PI * sample.b_bar / i - 1.0);
        assert_near(
            scaled_b,
            a,
            5.0e-5,
            &format!("{parity:?} symmetric edge coefficient"),
        );

        // W_i = 1/[2 pi^(3/2) sqrt(i)] * [1 - a/(16 i) + O(i^-2)].
        let scaled_weight = 16.0 * i * (2.0 * PI.powf(1.5) * i.sqrt() * sample.weight - 1.0);
        assert_near(
            scaled_weight,
            -a,
            2.0e-4,
            &format!("{parity:?} sign-corrected perturbation weight"),
        );

        // C_i = pi i W_i/sqrt(b_bar(i)) = 1 - a/(8 i) + O(i^-2).
        let scaled_prefactor = 8.0 * i * (sample.trace_prefactor - 1.0);
        assert_near(
            scaled_prefactor,
            -a,
            3.0e-4,
            &format!("{parity:?} local trace prefactor"),
        );
    }
}

#[test]
fn validates_intrinsic_soft_edge_infrared_scale() {
    // Combining the documented b_bar and V expansions gives
    // 4 i sqrt(V_i/b_bar(i)) = 1 - a/(8 i) + O(i^-2).
    // This is the intrinsic theta_IR ~ 1/(4 i) regulator used by the formal
    // local trace calculation; no fitted cutoff is introduced here.
    let row = 4096_usize;
    let i = row as f64;

    for parity in [ProlateParity::WPlus, ProlateParity::WMinus] {
        let a = parity_constant(parity);
        let sample = row_diagnostics(row, parity);
        let first_correction = i * (sample.infrared_scale - 1.0);

        assert_near(
            first_correction,
            -a / 8.0,
            2.0e-4,
            &format!("{parity:?} infrared first correction"),
        );
    }
}
