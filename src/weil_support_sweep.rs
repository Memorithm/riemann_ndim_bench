//! Deterministic support-window sweep for the finite Gram-normalized Weil audit.
//!
//! Each declared rational support window is evaluated independently at one
//! maximum basis dimension. Leading principal rows are then extracted from that
//! single `(A,G)` computation, reusing the dimension-sweep machinery.
//!
//! The support windows are numerical experiment parameters. No particular
//! window is assigned theoretical significance, and positivity on any finite
//! collection of windows does not imply Weil positivity or RH.

use std::fmt;

use crate::semilocal_compact_archimedean::{
    CompactArchimedeanBump, CompactArchimedeanError, PositiveRational,
};
use crate::weil_generalized_spectrum::{
    FiniteWeilGeneralizedSpectrumError, PrincipalWeilGeneralizedSpectrum,
    audit_finite_weil_generalized_spectrum,
};

/// Exact rational support window used by the sweep.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WeilSupportWindow {
    lower: PositiveRational,
    upper: PositiveRational,
}

impl WeilSupportWindow {
    #[inline]
    pub const fn new(lower: PositiveRational, upper: PositiveRational) -> Self {
        Self { lower, upper }
    }

    #[inline]
    pub const fn lower(self) -> PositiveRational {
        self.lower
    }

    #[inline]
    pub const fn upper(self) -> PositiveRational {
        self.upper
    }
}

/// One support window together with all leading-principal dimension rows.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWeilSupportWindowAudit {
    window: WeilSupportWindow,
    rows: Vec<PrincipalWeilGeneralizedSpectrum>,
    max_boundary_residual: f64,
    max_pairing_asymmetry: f64,
    max_whitened_asymmetry: f64,
}

impl FiniteWeilSupportWindowAudit {
    #[inline]
    pub const fn window(&self) -> WeilSupportWindow {
        self.window
    }

    #[inline]
    pub fn rows(&self) -> &[PrincipalWeilGeneralizedSpectrum] {
        &self.rows
    }

    #[inline]
    pub fn max_boundary_residual(&self) -> f64 {
        self.max_boundary_residual
    }

    #[inline]
    pub fn max_pairing_asymmetry(&self) -> f64 {
        self.max_pairing_asymmetry
    }

    #[inline]
    pub fn max_whitened_asymmetry(&self) -> f64 {
        self.max_whitened_asymmetry
    }
}

#[derive(Debug)]
pub enum FiniteWeilSupportSweepError {
    EmptyWindowSet,
    CompactSupport(CompactArchimedeanError),
    Generalized(FiniteWeilGeneralizedSpectrumError),
}

impl fmt::Display for FiniteWeilSupportSweepError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyWindowSet => write!(f, "finite Weil support sweep requires at least one window"),
            Self::CompactSupport(error) => write!(f, "invalid compact support window: {error}"),
            Self::Generalized(error) => write!(f, "generalized finite Weil audit failed: {error}"),
        }
    }
}

impl std::error::Error for FiniteWeilSupportSweepError {}

impl From<CompactArchimedeanError> for FiniteWeilSupportSweepError {
    fn from(value: CompactArchimedeanError) -> Self {
        Self::CompactSupport(value)
    }
}

impl From<FiniteWeilGeneralizedSpectrumError> for FiniteWeilSupportSweepError {
    fn from(value: FiniteWeilGeneralizedSpectrumError) -> Self {
        Self::Generalized(value)
    }
}

/// Evaluate every declared support window at `max_dimension`, then extract all
/// leading-principal dimension rows without recomputing pairings inside a
/// window.
pub fn audit_finite_weil_support_sweep(
    windows: &[WeilSupportWindow],
    max_dimension: usize,
    correlation_order: usize,
    archimedean_order: usize,
    boundary_order: usize,
    gram_order: usize,
) -> Result<Vec<FiniteWeilSupportWindowAudit>, FiniteWeilSupportSweepError> {
    if windows.is_empty() {
        return Err(FiniteWeilSupportSweepError::EmptyWindowSet);
    }

    let mut output = Vec::with_capacity(windows.len());
    for &window in windows {
        let bump = CompactArchimedeanBump::new(window.lower(), window.upper())?;
        let audit = audit_finite_weil_generalized_spectrum(
            bump,
            max_dimension,
            correlation_order,
            archimedean_order,
            boundary_order,
            gram_order,
        )?;
        let rows = audit.principal_sweep()?;
        output.push(FiniteWeilSupportWindowAudit {
            window,
            rows,
            max_boundary_residual: audit.pairing().max_boundary_residual(),
            max_pairing_asymmetry: audit.pairing().max_raw_pairing_asymmetry(),
            max_whitened_asymmetry: audit.max_whitened_asymmetry(),
        });
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rational(numerator: u64, denominator: u64) -> PositiveRational {
        PositiveRational::new(numerator, denominator).unwrap()
    }

    #[test]
    fn empty_window_set_is_rejected() {
        let error = audit_finite_weil_support_sweep(&[], 1, 32, 32, 48, 48).unwrap_err();
        assert!(matches!(error, FiniteWeilSupportSweepError::EmptyWindowSet));
    }

    #[test]
    fn invalid_window_order_is_rejected() {
        let windows = [WeilSupportWindow::new(rational(2, 1), rational(1, 1))];
        let error = audit_finite_weil_support_sweep(&windows, 1, 32, 32, 48, 48).unwrap_err();
        assert!(matches!(
            error,
            FiniteWeilSupportSweepError::CompactSupport(_)
        ));
    }

    #[test]
    fn support_sweep_returns_all_principal_rows_without_sign_assumptions() {
        let windows = [
            WeilSupportWindow::new(rational(1, 2), rational(7, 2)),
            WeilSupportWindow::new(rational(3, 4), rational(13, 4)),
        ];
        let audits = audit_finite_weil_support_sweep(&windows, 2, 48, 48, 64, 64).unwrap();

        assert_eq!(audits.len(), 2);
        for (audit, expected_window) in audits.iter().zip(windows) {
            assert_eq!(audit.window(), expected_window);
            assert_eq!(audit.rows().len(), 2);
            assert_eq!(audit.rows()[0].dimension(), 1);
            assert_eq!(audit.rows()[1].dimension(), 2);
            assert!(audit.max_boundary_residual().is_finite());
            assert!(audit.max_pairing_asymmetry().is_finite());
            assert!(audit.max_whitened_asymmetry().is_finite());
            for row in audit.rows() {
                assert!(row.raw_minimum_eigenvalue().is_finite());
                assert!(row.generalized_minimum_eigenvalue().is_finite());
                assert!(row.gram_minimum_eigenvalue().is_finite());
                assert!(row.gram_maximum_eigenvalue().is_finite());
                assert!(row.gram_condition_number().is_finite());
            }
        }
    }
}
