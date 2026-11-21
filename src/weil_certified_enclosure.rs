//! Rigorous, conservative integer enclosures for the ACTUAL compact
//! Riemann--Weil finite pairing of smooth Legendre-weighted Qg functions.
//!
//! The intervals are proved for the ideal C-infinity bump, NOT inferred from
//! binary64 quadrature or empirical error. They are intentionally very wide:
//! the exact sign checker will normally return Undetermined. No RH claim.

use std::fmt;

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_rational_gate::{
    ExactWeilSignAudit, ExactWeilSignError, WeilRationalInterval,
    audit_exact_weil_interval_sign,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CertifiedWeilEnvelopeError {
    EmptyBasis,
    DimensionOverflow,
    ArithmeticOverflow { stage: &'static str },
    ExactGate(ExactWeilSignError),
}

impl fmt::Display for CertifiedWeilEnvelopeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for CertifiedWeilEnvelopeError {}

impl From<ExactWeilSignError> for CertifiedWeilEnvelopeError {
    fn from(value: ExactWeilSignError) -> Self {
        Self::ExactGate(value)
    }
}

/// Analytic bounds for the ideal h_j=Qg_j in x=log(rho).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CertifiedWeilBasisBounds {
    sup_value: i128,
    sup_log_derivative: i128,
}

impl CertifiedWeilBasisBounds {
    pub const fn sup_value(self) -> i128 {
        self.sup_value
    }

    pub const fn sup_log_derivative(self) -> i128 {
        self.sup_log_derivative
    }
}

/// Exact rational intervals, all with common denominator one. The pairing
/// entries enclose the full ideal Weil functional, including archimedean and
/// all eligible prime terms. The zero-shift correlations have their own bounds.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CertifiedCompactWeilMatrixEnvelope {
    dimension: usize,
    log_span_upper: i128,
    ratio_upper: i128,
    max_strict_prime_argument: i128,
    possible_prime_power_count_upper: i128,
    arch_coefficient_abs_upper: i128,
    basis_bounds: Vec<CertifiedWeilBasisBounds>,
    theta_zero: Vec<WeilRationalInterval>,
    pairing: Vec<WeilRationalInterval>,
}

impl CertifiedCompactWeilMatrixEnvelope {
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    pub const fn scale(&self) -> i128 {
        1
    }

    pub const fn log_span_upper(&self) -> i128 {
        self.log_span_upper
    }

    pub const fn ratio_upper(&self) -> i128 {
        self.ratio_upper
    }

    pub const fn max_strict_prime_argument(&self) -> i128 {
        self.max_strict_prime_argument
    }

    pub const fn possible_prime_power_count_upper(&self) -> i128 {
        self.possible_prime_power_count_upper
    }

    pub const fn arch_coefficient_abs_upper(&self) -> i128 {
        self.arch_coefficient_abs_upper
    }

    pub fn basis_bounds(&self) -> &[CertifiedWeilBasisBounds] {
        &self.basis_bounds
    }

    pub fn theta_zero_intervals(&self) -> &[WeilRationalInterval] {
        &self.theta_zero
    }

    pub fn pairing_intervals(&self) -> &[WeilRationalInterval] {
        &self.pairing
    }

    /// Pass only independently established full-entry intervals to the
    /// checked-integer gate. A very wide interval cannot establish a sign.
    pub fn audit_exact_sign(
        &self,
        witness: &[i128],
        head_dimension: Option<usize>,
    ) -> Result<ExactWeilSignAudit, CertifiedWeilEnvelopeError> {
        Ok(audit_exact_weil_interval_sign(
            self.dimension,
            self.scale(),
            &self.pairing,
            witness,
            head_dimension,
        )?)
    }
}

fn overflow(stage: &'static str) -> CertifiedWeilEnvelopeError {
    CertifiedWeilEnvelopeError::ArithmeticOverflow { stage }
}

fn add(a: i128, b: i128, stage: &'static str) -> Result<i128, CertifiedWeilEnvelopeError> {
    a.checked_add(b).ok_or_else(|| overflow(stage))
}

fn sub(a: i128, b: i128, stage: &'static str) -> Result<i128, CertifiedWeilEnvelopeError> {
    a.checked_sub(b).ok_or_else(|| overflow(stage))
}

fn mul(a: i128, b: i128, stage: &'static str) -> Result<i128, CertifiedWeilEnvelopeError> {
    a.checked_mul(b).ok_or_else(|| overflow(stage))
}

fn ceil_ratio(
    numerator: u128,
    denominator: u128,
    stage: &'static str,
) -> Result<i128, CertifiedWeilEnvelopeError> {
    let quotient = numerator / denominator;
    let ceil = quotient
        .checked_add(u128::from(!numerator.is_multiple_of(denominator)))
        .ok_or_else(|| overflow(stage))?;
    i128::try_from(ceil).map_err(|_| overflow(stage))
}

/// Exact integer sup bounds for the ideal compact Qg_j profile.
///
/// Set d=t(1-t), y=1/d, B(t)=exp(-y), and a(t)=1-2t. The formulae
/// B_t=B*a*y^2,
/// B_tt=B*(-2*y^2-2*a^2*y^3+a^2*y^4),
/// B_ttt=B*(12*a*y^3+6*a^3*y^4-6*a*y^4-6*a^3*y^5+a^3*y^6)
/// and exp(-y)*y^k<=k! prove |B|<=1, |B_t|<=2, |B_tt|<=40,
/// |B_ttt|<=1800 on the entire smooth zero extension.
///
/// Since |P_j(z)|<=1 for z in [-1,1], iterated Markov inequalities give
/// |d_t^k P_j(2t-1)| <= (2*j^2)^k for k=1,2,3.
/// Put n=j^2. For g=B*P_j,
/// G1=2+2n, G2=40+8n+4n^2, G3=1800+240n+24n^2+8n^3
/// bound its first three t derivatives; |g|<=1.
/// With x=log rho and A=rho/(b-a)<=ceil(b/(b-a)), one has
/// h=-A*g_t-A^2*g_tt+g/4,
/// h_x=-(A*g_t+3A^2*g_tt+A^3*g_ttt)+A*g_t/4.
/// The safe integer bounds are H0=1+A*G1+A^2*G2 and
/// H1=2*A*G1+3*A^2*G2+A^3*G3. No f64 computation is involved.
pub fn certify_compact_weil_basis_bounds(
    bump: CompactArchimedeanBump,
    degree: usize,
) -> Result<CertifiedWeilBasisBounds, CertifiedWeilEnvelopeError> {
    let b_cross = u128::from(bump.upper().numerator())
        * u128::from(bump.lower().denominator());
    let a_cross = u128::from(bump.lower().numerator())
        * u128::from(bump.upper().denominator());
    let width = b_cross.checked_sub(a_cross).ok_or_else(|| overflow("width"))?;
    let a = ceil_ratio(b_cross, width, "coordinate scale")?;
    let degree = i128::try_from(degree).map_err(|_| overflow("degree"))?;
    let n = mul(degree, degree, "degree squared")?;
    let n2 = mul(n, n, "degree fourth")?;
    let n3 = mul(n2, n, "degree sixth")?;
    let g1 = add(2, mul(2, n, "G1 term")?, "G1")?;
    let g2 = add(
        add(40, mul(8, n, "G2 linear")?, "G2 subtotal")?,
        mul(4, n2, "G2 quadratic")?,
        "G2",
    )?;
    let g3 = add(
        add(
            add(1800, mul(240, n, "G3 linear")?, "G3 subtotal")?,
            mul(24, n2, "G3 quadratic")?,
            "G3 subtotal",
        )?,
        mul(8, n3, "G3 cubic")?,
        "G3",
    )?;
    let a2 = mul(a, a, "scale squared")?;
    let a3 = mul(a2, a, "scale cubed")?;
    let ag1 = mul(a, g1, "A G1")?;
    let a2g2 = mul(a2, g2, "A2 G2")?;
    let sup_value = add(add(1, ag1, "H0 subtotal")?, a2g2, "H0")?;
    let sup_log_derivative = add(
        add(
            mul(2, ag1, "H1 first")?,
            mul(3, a2g2, "H1 second")?,
            "H1 subtotal",
        )?,
        mul(a3, g3, "H1 third")?,
        "H1",
    )?;
    Ok(CertifiedWeilBasisBounds {
        sup_value,
        sup_log_derivative,
    })
}

/// Build proved (but conservative) FULL Weil matrix intervals.
///
/// Let R=b/a, T=ceil(R-1)>=log(R), C=ceil(R)>=R, and
/// K=17+ceil((R+1)/(R-1)). The exact Qg pole is zero.
/// With H0_i, H1_i from certify_compact_weil_basis_bounds:
///   |theta_ij(s)| <= T*H0_i*H0_j
///   |d_s(theta_ij(s)+theta_ij(-s))|
///       <= T*(H1_i*H0_j+H0_i*H1_j).
/// Both shift-derivative bounds follow by differentiating the smooth
/// compact convolution and integrating by parts. They apply for ALL s.
///
/// The actual source archimedean formula is
/// W_R=theta(0)*(gamma+log(4*pi*(R-1)/(R+1)))+int_0^log(R) I(s) ds,
/// I(s)=[(exp(s/2)-1)*2*theta(0)
///       +exp(s/2)*(theta(s)+theta(-s)-2*theta(0))]/(2*sinh s).
/// Since 0<gamma<1, 3<pi<4, |log z|<=max(z,1/z), the coefficient
/// is bounded by K. Also sinh(s)>=s, exp(s/2)<=R<=C and
/// exp(s/2)-1 <= C*s/2. These cancel the apparent singularity at
/// s=0, yielding the safe archimedean absolute bound
///   K*T*H0_i*H0_j
///   +C*T^2*(H0_i*H0_j+H1_i*H0_j+H0_i*H1_j).
///
/// Each prime-power contribution at 2<=n<R has absolute magnitude
/// <=2*T^2*H0_i*H0_j, because log(p)/sqrt(n)<=log(R)<=T.
/// At most N=max(0,ceil(R)-2) integers qualify. Therefore the complete
/// exact pairing lies in [-B,+B] with
/// B=arch_bound +2*N*T^2*H0_i*H0_j.
/// All arithmetic is exact and checked; overflows reject the certificate.
pub fn certify_compact_weil_matrix_enclosure(
    bump: CompactArchimedeanBump,
    dimension: usize,
) -> Result<CertifiedCompactWeilMatrixEnvelope, CertifiedWeilEnvelopeError> {
    if dimension == 0 {
        return Err(CertifiedWeilEnvelopeError::EmptyBasis);
    }
    let entry_count = dimension
        .checked_mul(dimension)
        .ok_or(CertifiedWeilEnvelopeError::DimensionOverflow)?;
    let b_cross = u128::from(bump.upper().numerator())
        * u128::from(bump.lower().denominator());
    let a_cross = u128::from(bump.lower().numerator())
        * u128::from(bump.upper().denominator());
    let difference = b_cross
        .checked_sub(a_cross)
        .ok_or_else(|| overflow("ratio difference"))?;
    let denominator = u128::from(bump.upper().denominator())
        * u128::from(bump.lower().numerator());
    let t = ceil_ratio(difference, denominator, "log span ceiling")?;
    let c = add(t, 1, "ratio ceiling")?;
    let sum = b_cross
        .checked_add(a_cross)
        .ok_or_else(|| overflow("coefficient ratio sum"))?;
    let k = add(
        17,
        ceil_ratio(sum, difference, "coefficient ceiling")?,
        "coefficient bound",
    )?;
    let strict = if b_cross.is_multiple_of(denominator) {
        b_cross / denominator - 1
    } else {
        b_cross / denominator
    };
    let max_strict = i128::try_from(strict).map_err(|_| overflow("prime cutoff"))?;
    let count = sub(max_strict, 1, "candidate prime count")?.max(0);

    let mut basis_bounds = Vec::new();
    basis_bounds
        .try_reserve_exact(dimension)
        .map_err(|_| CertifiedWeilEnvelopeError::DimensionOverflow)?;
    for degree in 0..dimension {
        basis_bounds.push(certify_compact_weil_basis_bounds(bump, degree)?);
    }
    let mut theta_zero = Vec::new();
    let mut pairing = Vec::new();
    theta_zero
        .try_reserve_exact(entry_count)
        .map_err(|_| CertifiedWeilEnvelopeError::DimensionOverflow)?;
    pairing
        .try_reserve_exact(entry_count)
        .map_err(|_| CertifiedWeilEnvelopeError::DimensionOverflow)?;
    let t2 = mul(t, t, "log span squared")?;
    let ct2 = mul(c, t2, "arch integral factor")?;
    let kt = mul(k, t, "arch coefficient factor")?;
    let prime_factor = mul(mul(2, count, "prime count")?, t2, "prime factor")?;

    for i in 0..dimension {
        for j in 0..dimension {
            let ui = basis_bounds[i];
            let uj = basis_bounds[j];
            let uv = mul(ui.sup_value(), uj.sup_value(), "value product")?;
            let ud = add(
                mul(ui.sup_log_derivative(), uj.sup_value(), "left derivative")?,
                mul(ui.sup_value(), uj.sup_log_derivative(), "right derivative")?,
                "derivative cross",
            )?;
            let correlation = mul(t, uv, "zero-shift correlation")?;
            theta_zero.push(if i == j {
                WeilRationalInterval::new(0, correlation)?
            } else {
                WeilRationalInterval::new(-correlation, correlation)?
            });
            let arch = add(
                mul(kt, uv, "arch coefficient")?,
                mul(ct2, add(uv, ud, "arch derivative")?, "arch integral")?,
                "arch total",
            )?;
            let bound = add(arch, mul(prime_factor, uv, "prime total")?, "pairing total")?;
            pairing.push(WeilRationalInterval::new(-bound, bound)?);
        }
    }
    Ok(CertifiedCompactWeilMatrixEnvelope {
        dimension,
        log_span_upper: t,
        ratio_upper: c,
        max_strict_prime_argument: max_strict,
        possible_prime_power_count_upper: count,
        arch_coefficient_abs_upper: k,
        basis_bounds,
        theta_zero,
        pairing,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::PositiveRational;
    use crate::weil_quadratic_matrix::audit_finite_weil_quadratic_matrix;
    use crate::weil_rational_gate::ExactFiniteWeilSign;

    fn rat(n: u64, d: u64) -> PositiveRational {
        PositiveRational::new(n, d).unwrap()
    }

    fn bump(a: (u64, u64), b: (u64, u64)) -> CompactArchimedeanBump {
        CompactArchimedeanBump::new(rat(a.0, a.1), rat(b.0, b.1)).unwrap()
    }

    #[test]
    fn ratio_two_eliminates_the_prime_sector_and_feeds_exact_gate() {
        let result = certify_compact_weil_matrix_enclosure(bump((1, 1), (2, 1)), 2).unwrap();
        assert_eq!(result.dimension(), 2);
        assert_eq!(result.scale(), 1);
        assert_eq!(result.log_span_upper(), 1);
        assert_eq!(result.ratio_upper(), 2);
        assert_eq!(result.arch_coefficient_abs_upper(), 20);
        assert_eq!(result.max_strict_prime_argument(), 1);
        assert_eq!(result.possible_prime_power_count_upper(), 0);
        assert_eq!(result.pairing_intervals().len(), 4);
        assert_eq!(result.pairing_intervals()[1], result.pairing_intervals()[2]);
        assert_eq!(result.theta_zero_intervals()[1], result.theta_zero_intervals()[2]);
        for i in 0..2 {
            let diagonal = result.theta_zero_intervals()[i * 2 + i];
            assert_eq!(diagonal.lower(), 0);
            assert!(diagonal.upper() > 0);
        }
        assert_eq!(
            result.audit_exact_sign(&[1, 0], Some(1)).unwrap().verdict(),
            ExactFiniteWeilSign::Undetermined
        );
    }

    #[test]
    fn exact_support_cutoff_handles_integer_and_noninteger_ratios() {
        let below = certify_compact_weil_matrix_enclosure(bump((1, 1), (3, 2)), 1).unwrap();
        assert_eq!(below.max_strict_prime_argument(), 1);
        assert_eq!(below.possible_prime_power_count_upper(), 0);
        let fractional = certify_compact_weil_matrix_enclosure(bump((1, 1), (5, 2)), 1).unwrap();
        assert_eq!(fractional.max_strict_prime_argument(), 2);
        assert_eq!(fractional.possible_prime_power_count_upper(), 1);
        let integer = certify_compact_weil_matrix_enclosure(bump((1, 1), (3, 1)), 1).unwrap();
        assert_eq!(integer.max_strict_prime_argument(), 2);
        assert_eq!(integer.possible_prime_power_count_upper(), 1);
    }

    #[test]
    fn analytical_intervals_contain_existing_numerical_observations() {
        // Numerical agreement is ONLY a regression check, not the proof.
        for support in [bump((1, 1), (2, 1)), bump((1, 2), (7, 2))] {
            let bounds = certify_compact_weil_matrix_enclosure(support, 2).unwrap();
            let observed = audit_finite_weil_quadratic_matrix(support, 2, 32, 32, 48).unwrap();
            for i in 0..2 {
                for j in 0..2 {
                    let interval = bounds.pairing_intervals()[i * 2 + j];
                    let entry = observed.entry(i, j).unwrap();
                    assert!(entry >= interval.lower() as f64);
                    assert!(entry <= interval.upper() as f64);
                }
            }
        }
    }

    #[test]
    fn degree_increase_preserves_non_decreasing_envelopes() {
        let support = bump((1, 1), (2, 1));
        let mut prev = certify_compact_weil_basis_bounds(support, 0).unwrap();
        for degree in 1..=5 {
            let next = certify_compact_weil_basis_bounds(support, degree).unwrap();
            assert!(next.sup_value() >= prev.sup_value());
            assert!(next.sup_log_derivative() >= prev.sup_log_derivative());
            prev = next;
        }
    }

    #[test]
    fn excessive_degree_or_narrow_support_fails_closed() {
        let support = bump((1, 1), (2, 1));
        assert!(matches!(
            certify_compact_weil_basis_bounds(support, usize::MAX),
            Err(CertifiedWeilEnvelopeError::ArithmeticOverflow { .. })
        ));
        let narrow = bump((1, 1), (1_000_000_000_001, 1_000_000_000_000));
        assert!(matches!(
            certify_compact_weil_matrix_enclosure(narrow, 2),
            Err(CertifiedWeilEnvelopeError::ArithmeticOverflow { .. })
        ));
        assert!(matches!(
            certify_compact_weil_matrix_enclosure(support, 0),
            Err(CertifiedWeilEnvelopeError::EmptyBasis)
        ));
    }
}
