use riemann_ndim_bench::semilocal::{CrossingDerivative, ProlateParity};
use riemann_ndim_bench::semilocal_tridiagonal::crossing_derivatives_tridiagonal;
use std::error::Error;
use std::io;
use std::time::{Duration, Instant};

#[cfg(target_os = "linux")]
use std::fs;

const BLOCK_SIZE: usize = 1024;
const EXPECTED_TOTAL_ABS_DERIVATIVE: f64 = 6.126_883_687_871;
const EXPECTED_M_MEAN_ABS: f64 = 0.095_732_557_622_99;
const EXPECTED_M_TRIMMED_MEAN_ABS: f64 = 0.024_908_301_053_19;
const EXPECTED_M_RMS: f64 = 0.410_603_983_216_0;
const EXPECTED_SQRT_M_LINF: f64 = 0.246_988_190_152_8;

#[derive(Debug)]
struct ParitySummary {
    min_mu: f64,
    max_mu: f64,
    min_lambda: f64,
    max_lambda: f64,
    min_lambda_prime: f64,
    max_lambda_prime: f64,
    negative_derivatives: usize,
    zero_derivatives: usize,
    positive_derivatives: usize,
    non_finite_values: usize,
}

#[derive(Debug)]
struct ResponseSummary {
    mean_abs: f64,
    trimmed_mean_abs: f64,
    rms: f64,
    linf: f64,
    total_abs_derivative: f64,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("semilocal high-m diagnostics failed: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let total_start = Instant::now();
    let rss_before_kib = peak_rss_kib();

    let plus_start = Instant::now();
    let plus = crossing_derivatives_tridiagonal(BLOCK_SIZE, ProlateParity::WPlus)?;
    let plus_elapsed = plus_start.elapsed();
    let rss_after_plus_kib = peak_rss_kib();

    let minus_start = Instant::now();
    let minus = crossing_derivatives_tridiagonal(BLOCK_SIZE, ProlateParity::WMinus)?;
    let minus_elapsed = minus_start.elapsed();
    let peak_rss_kib = peak_rss_kib();
    let total_elapsed = total_start.elapsed();

    let plus_summary = summarize_parity(&plus);
    let minus_summary = summarize_parity(&minus);
    let response = summarize_response(&plus, &minus, BLOCK_SIZE);

    ensure_finite_and_sign_correct(ProlateParity::WPlus, &plus_summary)?;
    ensure_finite_and_sign_correct(ProlateParity::WMinus, &minus_summary)?;

    assert_close(
        response.total_abs_derivative,
        EXPECTED_TOTAL_ABS_DERIVATIVE,
        5e-8,
        "m=1024 total_abs_derivative",
    )?;
    assert_close(
        BLOCK_SIZE as f64 * response.mean_abs,
        EXPECTED_M_MEAN_ABS,
        5e-8,
        "m=1024 m*mean_abs",
    )?;
    assert_close(
        BLOCK_SIZE as f64 * response.trimmed_mean_abs,
        EXPECTED_M_TRIMMED_MEAN_ABS,
        5e-8,
        "m=1024 m*trimmed_mean_abs",
    )?;
    assert_close(
        BLOCK_SIZE as f64 * response.rms,
        EXPECTED_M_RMS,
        5e-8,
        "m=1024 m*rms",
    )?;
    assert_close(
        (BLOCK_SIZE as f64).sqrt() * response.linf,
        EXPECTED_SQRT_M_LINF,
        5e-8,
        "m=1024 sqrt(m)*linf",
    )?;

    println!("semilocal_exact_q0_diagnostics");
    println!("block_size={BLOCK_SIZE}");
    println!("solver=faer_tridiagonal_self_adjoint_evd");
    println!("reference_solver=faer::linalg::solvers::SelfAdjointEigen");
    println!("independent_of_stieltjes_quadrature=true");
    println!("wplus_elapsed_seconds={:.6}", seconds(plus_elapsed));
    println!("wminus_elapsed_seconds={:.6}", seconds(minus_elapsed));
    println!("total_elapsed_seconds={:.6}", seconds(total_elapsed));
    print_optional_u64("rss_before_kib", rss_before_kib);
    print_optional_u64("rss_after_wplus_peak_kib", rss_after_plus_kib);
    print_optional_u64("process_peak_rss_kib", peak_rss_kib);
    print_parity("wplus", &plus_summary);
    print_parity("wminus", &minus_summary);
    println!(
        "merged_total_abs_derivative={:.15e}",
        response.total_abs_derivative
    );
    println!(
        "merged_m_mean_abs={:.15e}",
        BLOCK_SIZE as f64 * response.mean_abs
    );
    println!(
        "merged_m_trimmed_mean_abs={:.15e}",
        BLOCK_SIZE as f64 * response.trimmed_mean_abs
    );
    println!("merged_m_rms={:.15e}", BLOCK_SIZE as f64 * response.rms);
    println!(
        "merged_sqrt_m_linf={:.15e}",
        (BLOCK_SIZE as f64).sqrt() * response.linf
    );
    println!("finite_compression_only=true");
    println!("zeta_zero_identification=false");
    println!("rh_implication=false");

    Ok(())
}

fn summarize_parity(crossings: &[CrossingDerivative]) -> ParitySummary {
    let mut summary = ParitySummary {
        min_mu: f64::INFINITY,
        max_mu: f64::NEG_INFINITY,
        min_lambda: f64::INFINITY,
        max_lambda: f64::NEG_INFINITY,
        min_lambda_prime: f64::INFINITY,
        max_lambda_prime: f64::NEG_INFINITY,
        negative_derivatives: 0,
        zero_derivatives: 0,
        positive_derivatives: 0,
        non_finite_values: 0,
    };

    for crossing in crossings {
        if !crossing.mu.is_finite()
            || !crossing.lambda.is_finite()
            || !crossing.mu_prime.is_finite()
            || !crossing.lambda_prime.is_finite()
        {
            summary.non_finite_values += 1;
        }

        summary.min_mu = summary.min_mu.min(crossing.mu);
        summary.max_mu = summary.max_mu.max(crossing.mu);
        summary.min_lambda = summary.min_lambda.min(crossing.lambda);
        summary.max_lambda = summary.max_lambda.max(crossing.lambda);
        summary.min_lambda_prime = summary.min_lambda_prime.min(crossing.lambda_prime);
        summary.max_lambda_prime = summary.max_lambda_prime.max(crossing.lambda_prime);

        if crossing.lambda_prime < 0.0 {
            summary.negative_derivatives += 1;
        } else if crossing.lambda_prime > 0.0 {
            summary.positive_derivatives += 1;
        } else {
            summary.zero_derivatives += 1;
        }
    }

    summary
}

fn summarize_response(
    plus: &[CrossingDerivative],
    minus: &[CrossingDerivative],
    block_size: usize,
) -> ResponseSummary {
    let mut merged = plus.iter().chain(minus).collect::<Vec<_>>();
    merged.sort_by(|left, right| left.lambda.total_cmp(&right.lambda));

    let normalization = (block_size as f64).sqrt();
    let values = merged
        .iter()
        .map(|crossing| crossing.lambda_prime / normalization)
        .collect::<Vec<_>>();
    let count = values.len() as f64;
    let mean_abs = values.iter().map(|value| value.abs()).sum::<f64>() / count;
    let rms = (values.iter().map(|value| value * value).sum::<f64>() / count).sqrt();
    let linf = values
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let trim = values.len() / 8;
    let trimmed = &values[trim..values.len() - trim];
    let trimmed_mean_abs =
        trimmed.iter().map(|value| value.abs()).sum::<f64>() / trimmed.len() as f64;
    let total_abs_derivative = merged
        .iter()
        .map(|crossing| crossing.lambda_prime.abs())
        .sum();

    ResponseSummary {
        mean_abs,
        trimmed_mean_abs,
        rms,
        linf,
        total_abs_derivative,
    }
}

fn ensure_finite_and_sign_correct(
    parity: ProlateParity,
    summary: &ParitySummary,
) -> Result<(), Box<dyn Error>> {
    if summary.non_finite_values != 0 {
        return Err(io::Error::other(format!(
            "{parity:?}: observed {} non-finite spectral/derivative values",
            summary.non_finite_values
        ))
        .into());
    }

    let sign_ok = match parity {
        ProlateParity::WPlus => {
            summary.negative_derivatives == BLOCK_SIZE
                && summary.zero_derivatives == 0
                && summary.positive_derivatives == 0
        }
        ProlateParity::WMinus => {
            summary.positive_derivatives == BLOCK_SIZE
                && summary.zero_derivatives == 0
                && summary.negative_derivatives == 0
        }
    };

    if !sign_ok {
        return Err(io::Error::other(format!(
            "{parity:?}: derivative sign counts are negative={} zero={} positive={}",
            summary.negative_derivatives, summary.zero_derivatives, summary.positive_derivatives
        ))
        .into());
    }

    if !(summary.min_mu > 0.0 && summary.min_lambda > 0.0) {
        return Err(io::Error::other(format!(
            "{parity:?}: non-positive denominator diagnostics: min_mu={} min_lambda={}",
            summary.min_mu, summary.min_lambda
        ))
        .into());
    }

    Ok(())
}

fn assert_close(
    actual: f64,
    expected: f64,
    tolerance: f64,
    label: &str,
) -> Result<(), Box<dyn Error>> {
    let error = (actual - expected).abs();
    if error > tolerance {
        return Err(io::Error::other(format!(
            "{label}: actual={actual:.16e} expected={expected:.16e} error={error:.3e} tolerance={tolerance:.3e}"
        ))
        .into());
    }
    Ok(())
}

fn print_parity(label: &str, summary: &ParitySummary) {
    println!("{label}_min_mu={:.15e}", summary.min_mu);
    println!("{label}_max_mu={:.15e}", summary.max_mu);
    println!(
        "{label}_relative_min_mu={:.15e}",
        summary.min_mu / summary.max_mu
    );
    println!("{label}_min_lambda={:.15e}", summary.min_lambda);
    println!("{label}_max_lambda={:.15e}", summary.max_lambda);
    println!(
        "{label}_min_rayleigh_denominator={:.15e}",
        2.0 * summary.min_lambda
    );
    println!("{label}_min_lambda_prime={:.15e}", summary.min_lambda_prime);
    println!("{label}_max_lambda_prime={:.15e}", summary.max_lambda_prime);
    println!(
        "{label}_negative_derivatives={}",
        summary.negative_derivatives
    );
    println!("{label}_zero_derivatives={}", summary.zero_derivatives);
    println!(
        "{label}_positive_derivatives={}",
        summary.positive_derivatives
    );
    println!("{label}_non_finite_values={}", summary.non_finite_values);
}

fn seconds(duration: Duration) -> f64 {
    duration.as_secs_f64()
}

#[cfg(target_os = "linux")]
fn peak_rss_kib() -> Option<u64> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    status.lines().find_map(|line| {
        let rest = line.strip_prefix("VmHWM:")?;
        rest.split_whitespace().next()?.parse::<u64>().ok()
    })
}

#[cfg(not(target_os = "linux"))]
fn peak_rss_kib() -> Option<u64> {
    None
}

fn print_optional_u64(label: &str, value: Option<u64>) {
    match value {
        Some(value) => println!("{label}={value}"),
        None => println!("{label}=unavailable"),
    }
}
