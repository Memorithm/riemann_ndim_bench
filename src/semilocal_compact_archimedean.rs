//! Exact support bookkeeping for a manufactured compactly supported
//! archimedean factor used with the finite semilocal `E` bridge.
//!
//! The support endpoints and the positive archimedean coordinate are stored as
//! exact positive rationals. Integer `m` support bounds are therefore derived
//! without floating-point threshold decisions. Floating point is used only to
//! evaluate a smooth bump after exact support membership has been decided.

use std::fmt;

use crate::semilocal_bruhat_e::{
    BruhatEBridgeError, FiniteBruhatEBridgeAudit, compare_finite_bruhat_e_bridge,
};
use crate::semilocal_factorizable_poisson::LocalBallSpec;
use crate::semilocal_qs::{QsArithmeticError, QsRational};
use crate::semilocal_trace_contract::FinitePlaceSet;
use crate::weil_boundary::{MultiplicativeSupport, WeilBoundaryError};

/// Reduced strictly positive rational number `numerator / denominator`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PositiveRational {
    numerator: u64,
    denominator: u64,
}

impl PositiveRational {
    pub fn new(numerator: u64, denominator: u64) -> Result<Self, CompactArchimedeanError> {
        if numerator == 0 {
            return Err(CompactArchimedeanError::ZeroNumerator);
        }
        if denominator == 0 {
            return Err(CompactArchimedeanError::ZeroDenominator);
        }
        let divisor = gcd(numerator, denominator);
        Ok(Self {
            numerator: numerator / divisor,
            denominator: denominator / divisor,
        })
    }

    #[inline]
    pub fn numerator(self) -> u64 {
        self.numerator
    }

    #[inline]
    pub fn denominator(self) -> u64 {
        self.denominator
    }

    #[inline]
    pub fn as_f64(self) -> f64 {
        self.numerator as f64 / self.denominator as f64
    }

    /// Exact floor and remainder data for `self / scale`.
    fn quotient_by(self, scale: Self) -> RationalIntegerQuotient {
        let numerator = u128::from(self.numerator) * u128::from(scale.denominator);
        let denominator = u128::from(self.denominator) * u128::from(scale.numerator);
        RationalIntegerQuotient {
            floor: numerator / denominator,
            has_fraction: numerator % denominator != 0,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct RationalIntegerQuotient {
    floor: u128,
    has_fraction: bool,
}

impl RationalIntegerQuotient {
    /// Smallest integer strictly larger than the represented positive rational.
    fn first_integer_strictly_above(self) -> u128 {
        self.floor + 1
    }

    /// Largest integer strictly smaller than the represented positive rational.
    fn last_integer_strictly_below(self) -> Option<u128> {
        if self.has_fraction {
            Some(self.floor)
        } else {
            self.floor.checked_sub(1)
        }
    }
}

/// Standard smooth bump with an exact rational open support envelope.
///
/// Numerically, inside `(lower, upper)` the profile is
///
/// `exp(-1 / (t (1-t)))`, `t=(rho-lower)/(upper-lower)`,
///
/// and it is exactly zero outside the open interval. Exact rational arithmetic
/// decides whether a sampled point `rho = m * x_infinity` is inside the
/// support; the numerical formula is never used to decide truncation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactArchimedeanBump {
    lower: PositiveRational,
    upper: PositiveRational,
    support: MultiplicativeSupport,
}

impl CompactArchimedeanBump {
    pub fn new(
        lower: PositiveRational,
        upper: PositiveRational,
    ) -> Result<Self, CompactArchimedeanError> {
        let left = u128::from(lower.numerator) * u128::from(upper.denominator);
        let right = u128::from(upper.numerator) * u128::from(lower.denominator);
        if left >= right {
            return Err(CompactArchimedeanError::InvalidSupportOrder { lower, upper });
        }

        let support = MultiplicativeSupport::new(lower.as_f64(), upper.as_f64())?;
        Ok(Self {
            lower,
            upper,
            support,
        })
    }

    #[inline]
    pub fn lower(self) -> PositiveRational {
        self.lower
    }

    #[inline]
    pub fn upper(self) -> PositiveRational {
        self.upper
    }

    #[inline]
    pub fn support(self) -> MultiplicativeSupport {
        self.support
    }

    /// Exact inclusive integer range of `m >= 1` for which
    /// `lower < m*x_infinity < upper` can hold.
    pub fn active_m_bounds(self, x_infinity: PositiveRational) -> Option<(u64, u64)> {
        let lower_q = self.lower.quotient_by(x_infinity);
        let upper_q = self.upper.quotient_by(x_infinity);

        let min_u128 = lower_q.first_integer_strictly_above().max(1);
        let max_u128 = upper_q.last_integer_strictly_below()?;
        if min_u128 > max_u128 || min_u128 > u128::from(u64::MAX) {
            return None;
        }

        let min_m = min_u128 as u64;
        let max_m = max_u128.min(u128::from(u64::MAX)) as u64;
        (min_m <= max_m).then_some((min_m, max_m))
    }

    #[inline]
    pub fn max_m(self, x_infinity: PositiveRational) -> u64 {
        self.active_m_bounds(x_infinity)
            .map(|(_, max_m)| max_m)
            .unwrap_or(0)
    }

    /// Exact support-membership decision for the sampled point
    /// `rho = m*x_infinity`.
    pub fn contains_scaled_m(self, m: u64, x_infinity: PositiveRational) -> bool {
        self.active_m_bounds(x_infinity)
            .is_some_and(|(min_m, max_m)| m >= min_m && m <= max_m)
    }

    /// Numerical bump value after exact support membership has been decided.
    pub fn value_at_scaled_m(self, m: u64, x_infinity: PositiveRational) -> f64 {
        if !self.contains_scaled_m(m, x_infinity) {
            return 0.0;
        }

        let rho = m as f64 * x_infinity.as_f64();
        let lower = self.support.lower();
        let upper = self.support.upper();
        let width = upper - lower;
        let t = ((rho - lower) / width).clamp(f64::MIN_POSITIVE, 1.0 - f64::EPSILON);
        (-1.0 / (t * (1.0 - t))).exp()
    }
}

/// Result of the compact-support wrapper around the finite Bruhat `E` bridge.
#[derive(Clone, Debug, PartialEq)]
pub struct CompactBruhatEBridgeAudit {
    input_sample_count: usize,
    active_sample_count: usize,
    x_infinity: PositiveRational,
    active_m_bounds: Option<(u64, u64)>,
    inner: FiniteBruhatEBridgeAudit,
}

impl CompactBruhatEBridgeAudit {
    #[inline]
    pub fn input_sample_count(&self) -> usize {
        self.input_sample_count
    }

    #[inline]
    pub fn active_sample_count(&self) -> usize {
        self.active_sample_count
    }

    #[inline]
    pub fn x_infinity(&self) -> PositiveRational {
        self.x_infinity
    }

    #[inline]
    pub fn active_m_bounds(&self) -> Option<(u64, u64)> {
        self.active_m_bounds
    }

    #[inline]
    pub fn max_m(&self) -> u64 {
        self.active_m_bounds.map(|(_, max_m)| max_m).unwrap_or(0)
    }

    #[inline]
    pub fn inner(&self) -> &FiniteBruhatEBridgeAudit {
        &self.inner
    }
}

#[derive(Debug)]
pub enum CompactArchimedeanError {
    ZeroNumerator,
    ZeroDenominator,
    InvalidSupportOrder {
        lower: PositiveRational,
        upper: PositiveRational,
    },
    Support(WeilBoundaryError),
    Qs(QsArithmeticError),
    Bridge(BruhatEBridgeError),
}

impl fmt::Display for CompactArchimedeanError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroNumerator => write!(f, "positive rational numerator must be non-zero"),
            Self::ZeroDenominator => write!(f, "positive rational denominator must be non-zero"),
            Self::InvalidSupportOrder { lower, upper } => write!(
                f,
                "compact support must satisfy lower < upper: {}/{} !< {}/{}",
                lower.numerator,
                lower.denominator,
                upper.numerator,
                upper.denominator
            ),
            Self::Support(error) => write!(f, "compact support is not numerically representable: {error}"),
            Self::Qs(error) => write!(f, "Q_S arithmetic failed while applying support: {error}"),
            Self::Bridge(error) => write!(f, "finite Bruhat E bridge failed: {error}"),
        }
    }
}

impl std::error::Error for CompactArchimedeanError {}

impl From<WeilBoundaryError> for CompactArchimedeanError {
    fn from(value: WeilBoundaryError) -> Self {
        Self::Support(value)
    }
}

impl From<QsArithmeticError> for CompactArchimedeanError {
    fn from(value: QsArithmeticError) -> Self {
        Self::Qs(value)
    }
}

impl From<BruhatEBridgeError> for CompactArchimedeanError {
    fn from(value: BruhatEBridgeError) -> Self {
        Self::Bridge(value)
    }
}

/// Run the finite Bruhat `E` bridge with a smooth compact archimedean fixture.
///
/// `max_m` is derived exactly from the rational upper support endpoint and the
/// rational positive archimedean coordinate. Samples whose representative lies
/// outside the exact open support are removed before the finite bridge because
/// their archimedean contribution is certified to be zero.
pub fn compare_compact_bruhat_e_bridge(
    samples: &[QsRational],
    places: &FinitePlaceSet,
    original_balls: &[LocalBallSpec],
    modulus: f64,
    x_infinity: PositiveRational,
    bump: CompactArchimedeanBump,
) -> Result<CompactBruhatEBridgeAudit, CompactArchimedeanError> {
    let active_m_bounds = bump.active_m_bounds(x_infinity);
    let max_m = active_m_bounds.map(|(_, max_m)| max_m).unwrap_or(0);

    let mut active_samples = Vec::new();
    for &sample in samples {
        let decomposition = sample.unit_monoid_decomposition(places)?;
        if bump.contains_scaled_m(decomposition.monoid_element(), x_infinity) {
            active_samples.push(sample);
        }
    }

    let inner = compare_finite_bruhat_e_bridge(
        &active_samples,
        places,
        original_balls,
        modulus,
        max_m,
        |m| bump.value_at_scaled_m(m, x_infinity),
    )?;

    Ok(CompactBruhatEBridgeAudit {
        input_sample_count: samples.len(),
        active_sample_count: active_samples.len(),
        x_infinity,
        active_m_bounds,
        inner,
    })
}

fn gcd(mut left: u64, mut right: u64) -> u64 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_rational_support_derives_integer_m_bounds_without_float_thresholds() {
        let lower = PositiveRational::new(1, 2).unwrap();
        let upper = PositiveRational::new(7, 2).unwrap();
        let x = PositiveRational::new(1, 2).unwrap();
        let bump = CompactArchimedeanBump::new(lower, upper).unwrap();

        // 1/2 < m/2 < 7/2  <=>  1 < m < 7.
        assert_eq!(bump.active_m_bounds(x), Some((2, 6)));
        assert_eq!(bump.max_m(x), 6);
        assert!(!bump.contains_scaled_m(1, x));
        assert!(bump.contains_scaled_m(2, x));
        assert!(bump.contains_scaled_m(6, x));
        assert!(!bump.contains_scaled_m(7, x));
    }

    #[test]
    fn non_integer_upper_ratio_uses_exact_floor_for_last_active_m() {
        let lower = PositiveRational::new(1, 3).unwrap();
        let upper = PositiveRational::new(7, 3).unwrap();
        let x = PositiveRational::new(2, 5).unwrap();
        let bump = CompactArchimedeanBump::new(lower, upper).unwrap();

        // upper/x = (7/3)/(2/5) = 35/6, so m=5 is the last integer below it.
        assert_eq!(bump.max_m(x), 5);
    }

    #[test]
    fn compact_support_filters_boundary_samples_and_drives_finite_e_bound() {
        let places = FinitePlaceSet::new(vec![2]).unwrap();
        let balls = [LocalBallSpec::new(2, 0)];
        let samples = [
            QsRational::new(1, 1, &places).unwrap(),
            QsRational::new(3, 1, &places).unwrap(),
            QsRational::new(6, 1, &places).unwrap(),
            QsRational::new(3, 2, &places).unwrap(),
            QsRational::new(5, 1, &places).unwrap(),
            QsRational::new(10, 1, &places).unwrap(),
            QsRational::new(7, 1, &places).unwrap(),
        ];
        let x = PositiveRational::new(1, 2).unwrap();
        let bump = CompactArchimedeanBump::new(
            PositiveRational::new(1, 2).unwrap(),
            PositiveRational::new(7, 2).unwrap(),
        )
        .unwrap();

        let audit = compare_compact_bruhat_e_bridge(
            &samples,
            &places,
            &balls,
            4.0,
            x,
            bump,
        )
        .unwrap();

        // m=1 and m=7 lie exactly on the open support boundary and are removed.
        assert_eq!(audit.input_sample_count(), 7);
        assert_eq!(audit.active_sample_count(), 5);
        assert_eq!(audit.active_m_bounds(), Some((2, 6)));
        assert_eq!(audit.max_m(), 6);
        assert_eq!(audit.inner().finite_e_sum().max_m(), 6);
        assert_eq!(audit.inner().groups().len(), 2);
        assert_eq!(audit.inner().groups()[0].monoid_representative(), 3);
        assert_eq!(audit.inner().groups()[1].monoid_representative(), 5);
        assert!(audit.inner().e_bridge_residual().abs() <= 2.0e-15);
    }

    #[test]
    fn rational_inputs_are_reduced_canonically() {
        assert_eq!(
            PositiveRational::new(6, 8).unwrap(),
            PositiveRational::new(3, 4).unwrap()
        );
    }
}
