//! Exact local Fourier oracle for elementary p-adic balls.
//!
//! Under the source-locked self-dual additive character and Haar normalization
//! `vol(Z_p)=1`, the standard local identity is
//!
//! `F[1_{p^k Z_p}] = p^{-k} 1_{p^{-k} Z_p}`.
//!
//! This module records that identity symbolically.  It does not implement a
//! numerical p-adic Fourier transform, a general Bruhat--Schwartz function, or
//! the semilocal Fourier transform on `A_S`/`X_S`.

use std::fmt;

use crate::semilocal_qs::QsRational;
use crate::semilocal_trace_contract::FinitePlaceSet;

/// Ball `p^k Z_p` at one declared finite place.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PadicBall {
    prime: u64,
    exponent: i64,
}

impl PadicBall {
    /// Construct `p^k Z_p` for a prime already present in the source place set.
    pub fn new(
        prime: u64,
        exponent: i32,
        places: &FinitePlaceSet,
    ) -> Result<Self, PadicFourierError> {
        if !places.contains_prime(prime) {
            return Err(PadicFourierError::PrimeOutsidePlaceSet { prime });
        }
        Ok(Self {
            prime,
            exponent: i64::from(exponent),
        })
    }

    #[inline]
    fn from_transformed_parts(prime: u64, exponent: i64) -> Self {
        Self { prime, exponent }
    }

    #[inline]
    pub fn prime(self) -> u64 {
        self.prime
    }

    /// Integer `k` in `p^k Z_p`.
    #[inline]
    pub fn exponent(self) -> i64 {
        self.exponent
    }

    /// Exact membership of a diagonal `Q_S` sample in `p^k Z_p`.
    ///
    /// For non-zero `q`, membership is equivalent to `v_p(q) >= k`; zero
    /// belongs to every p-adic ball.
    pub fn contains_diagonal(self, q: QsRational) -> bool {
        if q.is_zero() {
            return true;
        }
        diagonal_valuation(q, self.prime) >= self.exponent
    }

    /// Symbolic self-dual local Fourier image
    /// `p^{-k} 1_{p^{-k} Z_p}`.
    pub fn fourier_transform(self) -> PadicBallFourierImage {
        PadicBallFourierImage {
            scale: PadicPowerScale {
                prime: self.prime,
                exponent: -self.exponent,
            },
            ball: Self::from_transformed_parts(self.prime, -self.exponent),
        }
    }
}

/// Exact symbolic multiplicative scale `p^e`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PadicPowerScale {
    prime: u64,
    exponent: i64,
}

impl PadicPowerScale {
    #[inline]
    pub fn prime(self) -> u64 {
        self.prime
    }

    #[inline]
    pub fn exponent(self) -> i64 {
        self.exponent
    }

    /// Whether the symbolic scale is exactly one.
    #[inline]
    pub fn is_one(self) -> bool {
        self.exponent == 0
    }

    fn multiply(self, other: Self) -> Result<Self, PadicFourierError> {
        if self.prime != other.prime {
            return Err(PadicFourierError::MismatchedScalePrimes {
                left: self.prime,
                right: other.prime,
            });
        }
        let exponent = self.exponent.checked_add(other.exponent).ok_or(
            PadicFourierError::ScaleExponentOverflow {
                left: self.exponent,
                right: other.exponent,
            },
        )?;
        Ok(Self {
            prime: self.prime,
            exponent,
        })
    }
}

/// Symbolic Fourier image `p^e 1_{p^k Z_p}` of an elementary local ball.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PadicBallFourierImage {
    scale: PadicPowerScale,
    ball: PadicBall,
}

impl PadicBallFourierImage {
    #[inline]
    pub fn scale(self) -> PadicPowerScale {
        self.scale
    }

    #[inline]
    pub fn ball(self) -> PadicBall {
        self.ball
    }

    /// Apply the local Fourier transform once more, preserving the accumulated
    /// symbolic scale exactly.
    pub fn fourier_transform(self) -> Result<Self, PadicFourierError> {
        let transformed = self.ball.fourier_transform();
        Ok(Self {
            scale: self.scale.multiply(transformed.scale)?,
            ball: transformed.ball,
        })
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PadicFourierError {
    PrimeOutsidePlaceSet { prime: u64 },
    MismatchedScalePrimes { left: u64, right: u64 },
    ScaleExponentOverflow { left: i64, right: i64 },
}

impl fmt::Display for PadicFourierError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PrimeOutsidePlaceSet { prime } => {
                write!(f, "p-adic Fourier place is not declared in S: p={prime}")
            }
            Self::MismatchedScalePrimes { left, right } => write!(
                f,
                "cannot multiply symbolic p-adic scales from different places: {left} and {right}"
            ),
            Self::ScaleExponentOverflow { left, right } => write!(
                f,
                "symbolic p-adic scale exponent overflow: {left} + {right}"
            ),
        }
    }
}

impl std::error::Error for PadicFourierError {}

fn diagonal_valuation(q: QsRational, prime: u64) -> i64 {
    valuation_u64(q.numerator_magnitude(), prime) as i64
        - valuation_u64(q.denominator(), prime) as i64
}

fn valuation_u64(mut value: u64, prime: u64) -> u32 {
    let mut valuation = 0_u32;
    while value != 0 && value.is_multiple_of(prime) {
        value /= prime;
        valuation += 1;
    }
    valuation
}
