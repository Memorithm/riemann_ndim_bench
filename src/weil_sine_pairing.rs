//! Full mixed-pairing cross-check for direct sine-enriched Weil functions.
//!
//! For each sine mode, the direct route evaluates `Q(bump sin(m*pi*t))`
//! analytically and sends that compact function through the generic mixed Weil
//! evaluator. The projected route uses the finite shifted-Legendre expansion
//! from `weil_sine_truncation`. Its total pairing is also checked against the
//! existing parent Legendre matrix through `C^T A C`.
//!
//! Agreement is finite numerical evidence about these two implementations. It
//! is not a density theorem, complete-space Weil positivity, Conjecture 4.1, or
//! a proof of RH.

use std::fmt;

use crate::semilocal_compact_archimedean::CompactArchimedeanBump;
use crate::weil_boundary::{WeilBoundaryError, WeilBoundaryMoments};
use crate::weil_compact_pairing::{
    CompactWeilEvaluand, CompactWeilPairingConfig, CompactWeilPairingError,
    FiniteCompactWeilPairingAudit, audit_compact_weil_pairing,
};
use crate::weil_quadratic_matrix::{
    CompactWeilBasisFunction, FiniteWeilMatrixError, FiniteWeilQuadraticMatrixAudit,
    audit_finite_weil_quadratic_matrix,
};
use crate::weil_sine_direct_q::{CompactSineWeilFunction, FiniteWeilDirectSineError};
use crate::weil_sine_truncation::{
    FiniteWeilSineTruncationError, SineModeSet, sine_legendre_coefficients,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DirectSinePairingAuditConfig {
    coefficient_quadrature_order: usize,
    pairing: CompactWeilPairingConfig,
}

impl DirectSinePairingAuditConfig {
    #[inline]
    pub const fn new(
        coefficient_quadrature_order: usize,
        pairing: CompactWeilPairingConfig,
    ) -> Self {
        Self {
            coefficient_quadrature_order,
            pairing,
        }
    }

    #[inline]
    pub const fn coefficient_quadrature_order(self) -> usize {
        self.coefficient_quadrature_order
    }

    #[inline]
    pub const fn pairing(self) -> CompactWeilPairingConfig {
        self.pairing
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DirectSinePairingSample {
    parent_dimension: usize,
    max_pairing_amplitude: f64,
    max_direct_projected_pairing_residual: f64,
    max_normalized_pairing_residual: f64,
    max_parent_matrix_projection_residual: f64,
    max_pole_term_residual: f64,
    max_archimedean_term_residual: f64,
    max_prime_total_residual: f64,
    max_direct_pairing_asymmetry: f64,
    max_projected_pairing_asymmetry: f64,
    max_direct_boundary_residual: f64,
    max_projected_boundary_residual: f64,
}

impl DirectSinePairingSample {
    #[inline]
    pub const fn parent_dimension(self) -> usize {
        self.parent_dimension
    }
    #[inline]
    pub const fn max_pairing_amplitude(self) -> f64 {
        self.max_pairing_amplitude
    }
    #[inline]
    pub const fn max_direct_projected_pairing_residual(self) -> f64 {
        self.max_direct_projected_pairing_residual
    }
    #[inline]
    pub const fn max_normalized_pairing_residual(self) -> f64 {
        self.max_normalized_pairing_residual
    }
    #[inline]
    pub const fn max_parent_matrix_projection_residual(self) -> f64 {
        self.max_parent_matrix_projection_residual
    }
    #[inline]
    pub const fn max_pole_term_residual(self) -> f64 {
        self.max_pole_term_residual
    }
    #[inline]
    pub const fn max_archimedean_term_residual(self) -> f64 {
        self.max_archimedean_term_residual
    }
    #[inline]
    pub const fn max_prime_total_residual(self) -> f64 {
        self.max_prime_total_residual
    }
    #[inline]
    pub const fn max_direct_pairing_asymmetry(self) -> f64 {
        self.max_direct_pairing_asymmetry
    }
    #[inline]
    pub const fn max_projected_pairing_asymmetry(self) -> f64 {
        self.max_projected_pairing_asymmetry
    }
    #[inline]
    pub const fn max_direct_boundary_residual(self) -> f64 {
        self.max_direct_boundary_residual
    }
    #[inline]
    pub const fn max_projected_boundary_residual(self) -> f64 {
        self.max_projected_boundary_residual
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilDirectSinePairingAudit {
    modes: SineModeSet,
    parent_dimensions: Vec<usize>,
    config: DirectSinePairingAuditConfig,
    samples: Vec<DirectSinePairingSample>,
}

impl FiniteWeilDirectSinePairingAudit {
    #[inline]
    pub fn modes(&self) -> &SineModeSet {
        &self.modes
    }
    #[inline]
    pub fn parent_dimensions(&self) -> &[usize] {
        &self.parent_dimensions
    }
    #[inline]
    pub const fn config(&self) -> DirectSinePairingAuditConfig {
        self.config
    }
    #[inline]
    pub fn samples(&self) -> &[DirectSinePairingSample] {
        &self.samples
    }

    pub fn last_pairing_residual_delta(&self) -> Option<f64> {
        let [.., previous, last] = self.samples.as_slice() else {
            return None;
        };
        Some(
            (last.max_direct_projected_pairing_residual
                - previous.max_direct_projected_pairing_residual)
                .abs(),
        )
    }
}

#[derive(Debug)]
pub enum FiniteWeilDirectSinePairingError {
    EmptyParentDimensionSet,
    ZeroParentDimension,
    ParentDimensionsNotStrictlyIncreasing { previous: usize, next: usize },
    DirectSine(FiniteWeilDirectSineError),
    SineProjection(FiniteWeilSineTruncationError),
    Pairing(CompactWeilPairingError),
    ParentMatrix(FiniteWeilMatrixError),
    Boundary(WeilBoundaryError),
    NonFiniteEvaluation { stage: &'static str, value: f64 },
}

impl fmt::Display for FiniteWeilDirectSinePairingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyParentDimensionSet => {
                write!(f, "direct sine pairing audit requires a parent dimension")
            }
            Self::ZeroParentDimension => {
                write!(f, "direct sine pairing parent dimension must be positive")
            }
            Self::ParentDimensionsNotStrictlyIncreasing { previous, next } => write!(
                f,
                "direct sine pairing parent dimensions must be strictly increasing: previous={previous}, next={next}"
            ),
            Self::DirectSine(error) => write!(f, "direct sine construction failed: {error}"),
            Self::SineProjection(error) => write!(f, "sine projection failed: {error}"),
            Self::Pairing(error) => write!(f, "compact Weil pairing failed: {error}"),
            Self::ParentMatrix(error) => write!(f, "parent Legendre pairing matrix failed: {error}"),
            Self::Boundary(error) => write!(f, "projected boundary evaluation failed: {error}"),
            Self::NonFiniteEvaluation { stage, value } => {
                write!(f, "non-finite direct sine pairing diagnostic at {stage}: {value}")
            }
        }
    }
}

impl std::error::Error for FiniteWeilDirectSinePairingError {}

impl From<FiniteWeilDirectSineError> for FiniteWeilDirectSinePairingError {
    fn from(value: FiniteWeilDirectSineError) -> Self {
        Self::DirectSine(value)
    }
}
impl From<FiniteWeilSineTruncationError> for FiniteWeilDirectSinePairingError {
    fn from(value: FiniteWeilSineTruncationError) -> Self {
        Self::SineProjection(value)
    }
}
impl From<CompactWeilPairingError> for FiniteWeilDirectSinePairingError {
    fn from(value: CompactWeilPairingError) -> Self {
        Self::Pairing(value)
    }
}
impl From<FiniteWeilMatrixError> for FiniteWeilDirectSinePairingError {
    fn from(value: FiniteWeilMatrixError) -> Self {
        Self::ParentMatrix(value)
    }
}
impl From<WeilBoundaryError> for FiniteWeilDirectSinePairingError {
    fn from(value: WeilBoundaryError) -> Self {
        Self::Boundary(value)
    }
}

#[derive(Clone, Debug)]
struct ProjectedLegendreWeilFunction {
    bump: CompactArchimedeanBump,
    coefficients: Vec<f64>,
}

impl CompactWeilEvaluand for ProjectedLegendreWeilFunction {
    #[inline]
    fn bump(&self) -> CompactArchimedeanBump {
        self.bump
    }

    fn value(&self, rho: f64) -> Result<f64, WeilBoundaryError> {
        let mut total = 0.0_f64;
        for (degree, &coefficient) in self.coefficients.iter().enumerate() {
            total += coefficient * CompactWeilBasisFunction::new(self.bump, degree).value(rho)?;
        }
        Ok(total)
    }

    fn boundary_moments(
        &self,
        quadrature_order: usize,
    ) -> Result<WeilBoundaryMoments, WeilBoundaryError> {
        let mut plus_half = 0.0_f64;
        let mut minus_half = 0.0_f64;
        for (degree, &coefficient) in self.coefficients.iter().enumerate() {
            let moments = CompactWeilBasisFunction::new(self.bump, degree)
                .boundary_moments(quadrature_order)?;
            plus_half += coefficient * moments.plus_half;
            minus_half += coefficient * moments.minus_half;
        }
        Ok(WeilBoundaryMoments {
            plus_half,
            minus_half,
        })
    }
}

#[derive(Clone, Copy)]
struct SymmetricPairingTerms {
    value: f64,
    pole_term: f64,
    archimedean_term: f64,
    prime_total: f64,
    max_boundary_residual: f64,
    asymmetry: f64,
}

impl SymmetricPairingTerms {
    fn from_forward_reverse(
        forward: FiniteCompactWeilPairingAudit,
        reverse: Option<FiniteCompactWeilPairingAudit>,
    ) -> Self {
        let max_boundary_residual = forward
            .left_boundary_residual()
            .max(forward.right_boundary_residual());
        if let Some(reverse) = reverse {
            Self {
                value: 0.5 * (forward.value() + reverse.value()),
                pole_term: 0.5 * (forward.pole_term() + reverse.pole_term()),
                archimedean_term: 0.5
                    * (forward.archimedean_term() + reverse.archimedean_term()),
                prime_total: 0.5 * (forward.prime_total() + reverse.prime_total()),
                max_boundary_residual: max_boundary_residual
                    .max(reverse.left_boundary_residual())
                    .max(reverse.right_boundary_residual()),
                asymmetry: (forward.value() - reverse.value()).abs(),
            }
        } else {
            Self {
                value: forward.value(),
                pole_term: forward.pole_term(),
                archimedean_term: forward.archimedean_term(),
                prime_total: forward.prime_total(),
                max_boundary_residual,
                asymmetry: 0.0,
            }
        }
    }
}

pub fn audit_finite_weil_direct_sine_pairing(
    bump: CompactArchimedeanBump,
    modes: &SineModeSet,
    parent_dimensions: &[usize],
    config: DirectSinePairingAuditConfig,
) -> Result<FiniteWeilDirectSinePairingAudit, FiniteWeilDirectSinePairingError> {
    validate_parent_dimensions(parent_dimensions)?;

    let direct_functions = modes
        .modes()
        .iter()
        .map(|&mode| CompactSineWeilFunction::new(bump, mode))
        .collect::<Result<Vec<_>, _>>()?;
    let pairing_config = config.pairing();
    let mut samples = Vec::with_capacity(parent_dimensions.len());

    for &parent_dimension in parent_dimensions {
        let parent_matrix = audit_finite_weil_quadratic_matrix(
            bump,
            parent_dimension,
            pairing_config.correlation_order(),
            pairing_config.archimedean_order(),
            pairing_config.boundary_order(),
        )?;
        let coefficients = modes
            .modes()
            .iter()
            .map(|&mode| {
                sine_legendre_coefficients(
                    mode,
                    parent_dimension,
                    config.coefficient_quadrature_order(),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let projected_functions = coefficients
            .iter()
            .cloned()
            .map(|coefficients| ProjectedLegendreWeilFunction { bump, coefficients })
            .collect::<Vec<_>>();

        let mut max_pairing_amplitude = 0.0_f64;
        let mut max_direct_projected_pairing_residual = 0.0_f64;
        let mut max_normalized_pairing_residual = 0.0_f64;
        let mut max_parent_matrix_projection_residual = 0.0_f64;
        let mut max_pole_term_residual = 0.0_f64;
        let mut max_archimedean_term_residual = 0.0_f64;
        let mut max_prime_total_residual = 0.0_f64;
        let mut max_direct_pairing_asymmetry = 0.0_f64;
        let mut max_projected_pairing_asymmetry = 0.0_f64;
        let mut max_direct_boundary_residual = 0.0_f64;
        let mut max_projected_boundary_residual = 0.0_f64;

        for i in 0..modes.dimension() {
            for j in i..modes.dimension() {
                let direct_forward = audit_compact_weil_pairing(
                    &direct_functions[i],
                    &direct_functions[j],
                    pairing_config,
                )?;
                let projected_forward = audit_compact_weil_pairing(
                    &projected_functions[i],
                    &projected_functions[j],
                    pairing_config,
                )?;
                let direct_reverse = if i == j {
                    None
                } else {
                    Some(audit_compact_weil_pairing(
                        &direct_functions[j],
                        &direct_functions[i],
                        pairing_config,
                    )?)
                };
                let projected_reverse = if i == j {
                    None
                } else {
                    Some(audit_compact_weil_pairing(
                        &projected_functions[j],
                        &projected_functions[i],
                        pairing_config,
                    )?)
                };

                let direct =
                    SymmetricPairingTerms::from_forward_reverse(direct_forward, direct_reverse);
                let projected = SymmetricPairingTerms::from_forward_reverse(
                    projected_forward,
                    projected_reverse,
                );
                let matrix_projected =
                    projected_matrix_entry(&parent_matrix, &coefficients[i], &coefficients[j]);
                let pairing_residual = (direct.value - projected.value).abs();
                let amplitude = direct.value.abs().max(projected.value.abs());
                let normalized = if amplitude > 0.0 {
                    pairing_residual / amplitude
                } else {
                    0.0
                };

                max_pairing_amplitude = max_pairing_amplitude.max(amplitude);
                max_direct_projected_pairing_residual =
                    max_direct_projected_pairing_residual.max(pairing_residual);
                max_normalized_pairing_residual = max_normalized_pairing_residual.max(normalized);
                max_parent_matrix_projection_residual = max_parent_matrix_projection_residual
                    .max((projected.value - matrix_projected).abs());
                max_pole_term_residual =
                    max_pole_term_residual.max((direct.pole_term - projected.pole_term).abs());
                max_archimedean_term_residual = max_archimedean_term_residual
                    .max((direct.archimedean_term - projected.archimedean_term).abs());
                max_prime_total_residual = max_prime_total_residual
                    .max((direct.prime_total - projected.prime_total).abs());
                max_direct_pairing_asymmetry =
                    max_direct_pairing_asymmetry.max(direct.asymmetry);
                max_projected_pairing_asymmetry =
                    max_projected_pairing_asymmetry.max(projected.asymmetry);
                max_direct_boundary_residual =
                    max_direct_boundary_residual.max(direct.max_boundary_residual);
                max_projected_boundary_residual =
                    max_projected_boundary_residual.max(projected.max_boundary_residual);
            }
        }

        for (stage, value) in [
            ("pairing amplitude", max_pairing_amplitude),
            (
                "direct/projected pairing residual",
                max_direct_projected_pairing_residual,
            ),
            ("normalized pairing residual", max_normalized_pairing_residual),
            (
                "parent matrix projection residual",
                max_parent_matrix_projection_residual,
            ),
            ("pole term residual", max_pole_term_residual),
            ("archimedean term residual", max_archimedean_term_residual),
            ("prime total residual", max_prime_total_residual),
            ("direct pairing asymmetry", max_direct_pairing_asymmetry),
            (
                "projected pairing asymmetry",
                max_projected_pairing_asymmetry,
            ),
            ("direct boundary residual", max_direct_boundary_residual),
            ("projected boundary residual", max_projected_boundary_residual),
        ] {
            checked_finite(stage, value)?;
        }

        samples.push(DirectSinePairingSample {
            parent_dimension,
            max_pairing_amplitude,
            max_direct_projected_pairing_residual,
            max_normalized_pairing_residual,
            max_parent_matrix_projection_residual,
            max_pole_term_residual,
            max_archimedean_term_residual,
            max_prime_total_residual,
            max_direct_pairing_asymmetry,
            max_projected_pairing_asymmetry,
            max_direct_boundary_residual,
            max_projected_boundary_residual,
        });
    }

    Ok(FiniteWeilDirectSinePairingAudit {
        modes: modes.clone(),
        parent_dimensions: parent_dimensions.to_vec(),
        config,
        samples,
    })
}

fn projected_matrix_entry(
    matrix: &FiniteWeilQuadraticMatrixAudit,
    left: &[f64],
    right: &[f64],
) -> f64 {
    let dimension = matrix.dimension();
    debug_assert_eq!(left.len(), dimension);
    debug_assert_eq!(right.len(), dimension);
    let mut total = 0.0_f64;
    for (i, &left_coefficient) in left.iter().enumerate() {
        for (j, &right_coefficient) in right.iter().enumerate() {
            total += left_coefficient
                * matrix
                    .entry(i, j)
                    .expect("coefficient indices lie inside parent pairing matrix")
                * right_coefficient;
        }
    }
    total
}

fn validate_parent_dimensions(
    parent_dimensions: &[usize],
) -> Result<(), FiniteWeilDirectSinePairingError> {
    if parent_dimensions.is_empty() {
        return Err(FiniteWeilDirectSinePairingError::EmptyParentDimensionSet);
    }
    if parent_dimensions.contains(&0) {
        return Err(FiniteWeilDirectSinePairingError::ZeroParentDimension);
    }
    for pair in parent_dimensions.windows(2) {
        if pair[0] >= pair[1] {
            return Err(
                FiniteWeilDirectSinePairingError::ParentDimensionsNotStrictlyIncreasing {
                    previous: pair[0],
                    next: pair[1],
                },
            );
        }
    }
    Ok(())
}

fn checked_finite(
    stage: &'static str,
    value: f64,
) -> Result<(), FiniteWeilDirectSinePairingError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FiniteWeilDirectSinePairingError::NonFiniteEvaluation { stage, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_compact_archimedean::PositiveRational;

    fn bump() -> CompactArchimedeanBump {
        CompactArchimedeanBump::new(
            PositiveRational::new(1, 2).unwrap(),
            PositiveRational::new(7, 2).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn direct_sine_pairing_audit_records_finite_cross_route_diagnostics() {
        let modes = SineModeSet::new(vec![1, 2]).unwrap();
        let config = DirectSinePairingAuditConfig::new(
            48,
            CompactWeilPairingConfig::new(20, 20, 28),
        );
        let audit =
            audit_finite_weil_direct_sine_pairing(bump(), &modes, &[3, 5], config).unwrap();

        assert_eq!(audit.modes(), &modes);
        assert_eq!(audit.parent_dimensions(), &[3, 5]);
        assert_eq!(audit.samples().len(), 2);
        assert!(audit.last_pairing_residual_delta().unwrap().is_finite());
        for sample in audit.samples().iter().copied() {
            assert!(sample.max_pairing_amplitude().is_finite());
            assert!(sample.max_direct_projected_pairing_residual().is_finite());
            assert!(sample.max_normalized_pairing_residual().is_finite());
            assert!(sample.max_parent_matrix_projection_residual().is_finite());
            assert!(sample.max_pole_term_residual().is_finite());
            assert!(sample.max_archimedean_term_residual().is_finite());
            assert!(sample.max_prime_total_residual().is_finite());
            assert!(sample.max_direct_pairing_asymmetry().is_finite());
            assert!(sample.max_projected_pairing_asymmetry().is_finite());
            assert!(sample.max_direct_boundary_residual().is_finite());
            assert!(sample.max_projected_boundary_residual().is_finite());
            assert!(sample.max_parent_matrix_projection_residual() <= 2.0e-10);
        }
    }
}
