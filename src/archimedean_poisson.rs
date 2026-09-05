//! Certified manufactured oracle for the archimedean Poisson identity used by
//! the Riemann-specific semilocal bridge.
//!
//! Connes--Consani, *The Scaling Hamiltonian*, equation (4.3), uses the
//! archimedean Poisson relation
//!
//! `E(F f)(x) = E(f)(x^-1)`, `x > 0`,
//!
//! where, for an even function satisfying `f(0)=Ff(0)=0`,
//!
//! `E(f)(x) = x^(1/2) sum_{n>=1} f(n x)`.
//!
//! This module does not attempt a generic Fourier library. It provides one
//! analytic Gaussian-combination fixture whose Fourier transform and both
//! boundary values are known in closed form. Finite E-sums carry explicit
//! Gaussian tail bounds, so the regression tests the source Poisson mechanism
//! without invoking zeta or a fitted numerical target.

use std::f64::consts::PI;
use std::fmt;

#[derive(Clone, Copy, Debug)]
struct GaussianComponent {
    coefficient: f64,
    scale: f64,
}

const SOURCE_FIXTURE: [GaussianComponent; 3] = [
    GaussianComponent {
        coefficient: 1.0 / 3.0,
        scale: 1.0,
    },
    GaussianComponent {
        coefficient: -4.0 / 3.0,
        scale: 4.0,
    },
    GaussianComponent {
        coefficient: 1.0,
        scale: 9.0,
    },
];

/// Error returned by the certified archimedean Poisson oracle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ArchimedeanPoissonError {
    InvalidScalePoint { x: f64 },
    EmptyPrefix,
}

impl fmt::Display for ArchimedeanPoissonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidScalePoint { x } => {
                write!(
                    f,
                    "Poisson scale x must be finite and strictly positive: {x}"
                )
            }
            Self::EmptyPrefix => write!(f, "Poisson E-sum prefix must contain at least one term"),
        }
    }
}

impl std::error::Error for ArchimedeanPoissonError {}

/// Exact manufactured source fixture
///
/// `f(x) = (1/3)e^(-pi x^2) - (4/3)e^(-4 pi x^2) + e^(-9 pi x^2)`.
///
/// Its coefficients were chosen so that both `f(0)` and `Ff(0)` vanish under
/// the Fourier convention `Ff(xi)=integral f(x)e^(-2 pi i x xi) dx`.
pub fn source_poisson_fixture_value(x: f64) -> f64 {
    mixture_value(&SOURCE_FIXTURE, x)
}

/// Closed Fourier transform of [`source_poisson_fixture_value`].
///
/// Since
///
/// `F[e^(-pi a x^2)](xi) = a^(-1/2)e^(-pi xi^2/a)`,
///
/// the transformed fixture is obtained term by term with no numerical Fourier
/// quadrature.
pub fn source_poisson_fixture_fourier_value(xi: f64) -> f64 {
    SOURCE_FIXTURE
        .iter()
        .map(|component| {
            let transformed = GaussianComponent {
                coefficient: component.coefficient / component.scale.sqrt(),
                scale: component.scale.recip(),
            };
            gaussian_component_value(transformed, xi)
        })
        .sum()
}

/// Boundary value `f(0)`, exposed to keep the source hypothesis auditable.
pub fn source_poisson_fixture_at_zero() -> f64 {
    SOURCE_FIXTURE
        .iter()
        .map(|component| component.coefficient)
        .sum()
}

/// Boundary value `Ff(0)`, exposed to keep the source hypothesis auditable.
pub fn source_poisson_fixture_fourier_at_zero() -> f64 {
    SOURCE_FIXTURE
        .iter()
        .map(|component| component.coefficient / component.scale.sqrt())
        .sum()
}

/// Auditable finite approximation of one archimedean `E` sum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CertifiedESum {
    value: f64,
    max_n: u64,
    absolute_tail_bound: f64,
}

impl CertifiedESum {
    #[inline]
    pub fn value(self) -> f64 {
        self.value
    }

    #[inline]
    pub fn max_n(self) -> u64 {
        self.max_n
    }

    #[inline]
    pub fn absolute_tail_bound(self) -> f64 {
        self.absolute_tail_bound
    }
}

/// Certified comparison of the two sides of equation (4.3).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PoissonIdentityComparison {
    x: f64,
    e_fourier_at_x: CertifiedESum,
    e_original_at_inverse_x: CertifiedESum,
}

impl PoissonIdentityComparison {
    #[inline]
    pub fn x(self) -> f64 {
        self.x
    }

    #[inline]
    pub fn left(self) -> CertifiedESum {
        self.e_fourier_at_x
    }

    #[inline]
    pub fn right(self) -> CertifiedESum {
        self.e_original_at_inverse_x
    }

    #[inline]
    pub fn residual(self) -> f64 {
        self.e_fourier_at_x.value - self.e_original_at_inverse_x.value
    }

    /// Sum of the rigorous truncation bounds on both sides.
    #[inline]
    pub fn combined_tail_bound(self) -> f64 {
        self.e_fourier_at_x.absolute_tail_bound + self.e_original_at_inverse_x.absolute_tail_bound
    }
}

/// Compare the two sides of the source Poisson identity with certified Gaussian
/// truncation bounds.
pub fn compare_source_poisson_identity(
    x: f64,
    max_n: u64,
) -> Result<PoissonIdentityComparison, ArchimedeanPoissonError> {
    checked_x(x)?;
    if max_n == 0 {
        return Err(ArchimedeanPoissonError::EmptyPrefix);
    }

    let transformed = transformed_components();
    let e_fourier_at_x = certified_e_sum(&transformed, x, max_n)?;
    let e_original_at_inverse_x = certified_e_sum(&SOURCE_FIXTURE, x.recip(), max_n)?;

    Ok(PoissonIdentityComparison {
        x,
        e_fourier_at_x,
        e_original_at_inverse_x,
    })
}

fn transformed_components() -> [GaussianComponent; 3] {
    SOURCE_FIXTURE.map(|component| GaussianComponent {
        coefficient: component.coefficient / component.scale.sqrt(),
        scale: component.scale.recip(),
    })
}

fn certified_e_sum(
    components: &[GaussianComponent],
    x: f64,
    max_n: u64,
) -> Result<CertifiedESum, ArchimedeanPoissonError> {
    checked_x(x)?;
    if max_n == 0 {
        return Err(ArchimedeanPoissonError::EmptyPrefix);
    }

    let mut raw_sum = 0.0;
    for n in 1..=max_n {
        raw_sum += mixture_value(components, n as f64 * x);
    }

    let half_density = x.sqrt();
    let absolute_tail_bound = half_density
        * components
            .iter()
            .map(|&component| gaussian_series_tail_bound(component, x, max_n))
            .sum::<f64>();

    Ok(CertifiedESum {
        value: half_density * raw_sum,
        max_n,
        absolute_tail_bound,
    })
}

fn gaussian_series_tail_bound(component: GaussianComponent, x: f64, max_n: u64) -> f64 {
    let c = PI * component.scale * x * x;
    let n = max_n as f64;
    component.coefficient.abs() * (-c * n * n).exp() / (2.0 * c * n)
}

fn mixture_value(components: &[GaussianComponent], x: f64) -> f64 {
    components
        .iter()
        .map(|&component| gaussian_component_value(component, x))
        .sum()
}

fn gaussian_component_value(component: GaussianComponent, x: f64) -> f64 {
    component.coefficient * (-PI * component.scale * x * x).exp()
}

fn checked_x(x: f64) -> Result<(), ArchimedeanPoissonError> {
    if x.is_finite() && x > 0.0 {
        Ok(())
    } else {
        Err(ArchimedeanPoissonError::InvalidScalePoint { x })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixture_satisfies_both_source_boundary_conditions() {
        assert!(source_poisson_fixture_at_zero().abs() < 2.0e-16);
        assert!(source_poisson_fixture_fourier_at_zero().abs() < 2.0e-16);
    }

    #[test]
    fn invalid_scale_and_empty_prefix_are_rejected() {
        assert!(matches!(
            compare_source_poisson_identity(0.0, 16),
            Err(ArchimedeanPoissonError::InvalidScalePoint { .. })
        ));
        assert_eq!(
            compare_source_poisson_identity(1.0, 0).unwrap_err(),
            ArchimedeanPoissonError::EmptyPrefix
        );
    }
}
