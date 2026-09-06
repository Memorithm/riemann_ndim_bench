//! Exact `Q_S^*` unit action behind the quotient descent to `M_S`.
//!
//! Every non-zero `q in Q_S` has the source decomposition `q = u m`, with
//! `u in Q_S^*` and `m in M_S`. The finite-place part of `u` cannot be dropped:
//! for `u = sign * prod_p p^{n_p}`, multiplication of a local argument by `u`
//! transports
//!
//! `1_{p^k Z_p}(u_p x_p) = 1_{p^{k-n_p} Z_p}(x_p)`.
//!
//! At the same time the absolute archimedean scale is multiplied by
//! `prod_p p^{n_p}`. Therefore the effective exponent carried by the product
//! of archimedean and finite-place lattice scales is preserved exactly:
//!
//! `n_p + (k_p - n_p) = k_p`.
//!
//! This finite algebra is the compensation erased by the scalar shortcut
//! rejected in PR #43. It is not yet a construction of `X_S`, a proof of the
//! quotient Poisson identity, Conjecture 4.1, Weil positivity, or RH.

use std::fmt;

use crate::semilocal_factorizable_poisson::LocalBallSpec;
use crate::semilocal_qs::QsUnitMonoidDecomposition;
use crate::semilocal_trace_contract::FinitePlaceSet;

/// Exact action of the unit part `u in Q_S^*` on a factorizable local fixture.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SemilocalUnitOrbitTransport {
    unit_sign: i8,
    unit_exponents: Vec<(u64, i32)>,
}

impl SemilocalUnitOrbitTransport {
    /// Build the action from an already source-validated `q = u m`
    /// decomposition and its declared finite place set.
    pub fn from_decomposition(
        decomposition: &QsUnitMonoidDecomposition,
        places: &FinitePlaceSet,
    ) -> Result<Self, UnitOrbitError> {
        if decomposition.unit_sign() != -1 && decomposition.unit_sign() != 1 {
            return Err(UnitOrbitError::InvalidUnitSign {
                sign: decomposition.unit_sign(),
            });
        }

        if decomposition.unit_exponents().len() != places.finite_primes().len() {
            return Err(UnitOrbitError::PlaceSetMismatch);
        }
        for (&prime, &(candidate, _)) in places
            .finite_primes()
            .iter()
            .zip(decomposition.unit_exponents().iter())
        {
            if prime != candidate {
                return Err(UnitOrbitError::PlaceSetMismatch);
            }
        }

        Ok(Self {
            unit_sign: decomposition.unit_sign(),
            unit_exponents: decomposition.unit_exponents().to_vec(),
        })
    }

    /// Sign of the rational unit. The manufactured archimedean fixture is
    /// even, but the sign is retained rather than silently discarded.
    #[inline]
    pub fn unit_sign(&self) -> i8 {
        self.unit_sign
    }

    /// Exact finite-place valuation signature `n_p = v_p(u)`.
    #[inline]
    pub fn unit_exponents(&self) -> &[(u64, i32)] {
        &self.unit_exponents
    }

    /// Exact exponent `n_p` of one declared finite place.
    pub fn exponent_for(&self, prime: u64) -> Option<i32> {
        self.unit_exponents
            .binary_search_by_key(&prime, |&(candidate, _)| candidate)
            .ok()
            .map(|index| self.unit_exponents[index].1)
    }

    /// Transport one local ball under `x -> u x`: `k_p -> k_p - n_p`.
    pub fn transport_ball(&self, ball: LocalBallSpec) -> Result<LocalBallSpec, UnitOrbitError> {
        let unit_exponent =
            self.exponent_for(ball.prime())
                .ok_or(UnitOrbitError::PrimeOutsideUnitPlaces {
                    prime: ball.prime(),
                })?;
        let exponent =
            ball.exponent()
                .checked_sub(unit_exponent)
                .ok_or(UnitOrbitError::ExponentOverflow {
                    prime: ball.prime(),
                    ball_exponent: ball.exponent(),
                    unit_exponent,
                })?;
        Ok(LocalBallSpec::new(ball.prime(), exponent))
    }

    /// Transport a complete local product, requiring exactly one ball for each
    /// finite place represented by the source unit decomposition.
    pub fn transport_complete_product(
        &self,
        balls: &[LocalBallSpec],
    ) -> Result<Vec<LocalBallSpec>, UnitOrbitError> {
        if balls.len() != self.unit_exponents.len() {
            return Err(UnitOrbitError::IncompleteLocalProduct);
        }

        let mut sorted = balls.to_vec();
        sorted.sort_unstable_by_key(|ball| ball.prime());
        for (index, &(prime, _)) in self.unit_exponents.iter().enumerate() {
            if sorted.get(index).map(|ball| ball.prime()) != Some(prime) {
                return Err(UnitOrbitError::IncompleteLocalProduct);
            }
        }

        sorted
            .into_iter()
            .map(|ball| self.transport_ball(ball))
            .collect()
    }

    /// Verify the exact exponent bookkeeping at one finite place:
    /// `n_p + (k_p - n_p) = k_p`.
    pub fn compensated_exponent(
        &self,
        original_ball: LocalBallSpec,
    ) -> Result<i32, UnitOrbitError> {
        let unit_exponent = self.exponent_for(original_ball.prime()).ok_or(
            UnitOrbitError::PrimeOutsideUnitPlaces {
                prime: original_ball.prime(),
            },
        )?;
        let transported = self.transport_ball(original_ball)?;
        unit_exponent.checked_add(transported.exponent()).ok_or(
            UnitOrbitError::CompensationOverflow {
                prime: original_ball.prime(),
                unit_exponent,
                transported_exponent: transported.exponent(),
            },
        )
    }

    /// Absolute archimedean scale `|u_infinity| = prod_p p^{n_p}` evaluated in
    /// `f64`. Exact reasoning should use [`Self::unit_exponents`]; this value is
    /// only a numerical audit hook for manufactured fixtures.
    ///
    /// The implementation multiplies/divides the exact integer prime factors
    /// directly instead of taking `exp(sum n_p log p)`, avoiding unnecessary
    /// transcendental roundoff in exactly representable cases such as `9/8`.
    pub fn archimedean_absolute_scale(&self) -> Result<f64, UnitOrbitError> {
        let mut scale = 1.0_f64;
        for &(prime, exponent) in &self.unit_exponents {
            let factor = prime as f64;
            if exponent >= 0 {
                for _ in 0..exponent as u32 {
                    scale *= factor;
                    if !scale.is_finite() {
                        return Err(UnitOrbitError::ArchimedeanScaleOutOfRange);
                    }
                }
            } else {
                for _ in 0..exponent.unsigned_abs() {
                    scale /= factor;
                    if scale == 0.0 {
                        return Err(UnitOrbitError::ArchimedeanScaleOutOfRange);
                    }
                }
            }
        }
        if scale.is_finite() && scale > 0.0 {
            Ok(scale)
        } else {
            Err(UnitOrbitError::ArchimedeanScaleOutOfRange)
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum UnitOrbitError {
    InvalidUnitSign {
        sign: i8,
    },
    PlaceSetMismatch,
    PrimeOutsideUnitPlaces {
        prime: u64,
    },
    IncompleteLocalProduct,
    ExponentOverflow {
        prime: u64,
        ball_exponent: i32,
        unit_exponent: i32,
    },
    CompensationOverflow {
        prime: u64,
        unit_exponent: i32,
        transported_exponent: i32,
    },
    ArchimedeanScaleOutOfRange,
}

impl fmt::Display for UnitOrbitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidUnitSign { sign } => write!(f, "Q_S unit sign must be +/-1: {sign}"),
            Self::PlaceSetMismatch => write!(
                f,
                "Q_S unit decomposition does not match the declared finite place set"
            ),
            Self::PrimeOutsideUnitPlaces { prime } => {
                write!(f, "local ball prime is absent from the Q_S unit: p={prime}")
            }
            Self::IncompleteLocalProduct => write!(
                f,
                "local product must contain exactly one ball for every finite place"
            ),
            Self::ExponentOverflow {
                prime,
                ball_exponent,
                unit_exponent,
            } => write!(
                f,
                "local exponent overflow at p={prime}: {ball_exponent} - {unit_exponent}"
            ),
            Self::CompensationOverflow {
                prime,
                unit_exponent,
                transported_exponent,
            } => write!(
                f,
                "compensated exponent overflow at p={prime}: {unit_exponent} + {transported_exponent}"
            ),
            Self::ArchimedeanScaleOutOfRange => write!(
                f,
                "Q_S unit archimedean scale is outside finite positive f64 range"
            ),
        }
    }
}

impl std::error::Error for UnitOrbitError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semilocal_qs::QsRational;

    #[test]
    fn source_decomposition_transports_local_balls_exactly() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let q = QsRational::new(45, 8, &places).unwrap();
        let decomposition = q.unit_monoid_decomposition(&places).unwrap();
        let action =
            SemilocalUnitOrbitTransport::from_decomposition(&decomposition, &places).unwrap();

        assert_eq!(decomposition.monoid_element(), 5);
        assert_eq!(action.unit_exponents(), &[(2, -3), (3, 2)]);

        let original = [LocalBallSpec::new(2, 1), LocalBallSpec::new(3, -2)];
        let transported = action.transport_complete_product(&original).unwrap();
        assert_eq!(
            transported,
            vec![LocalBallSpec::new(2, 4), LocalBallSpec::new(3, -4)]
        );
        assert_eq!(action.compensated_exponent(original[0]).unwrap(), 1);
        assert_eq!(action.compensated_exponent(original[1]).unwrap(), -2);
    }

    #[test]
    fn archimedean_unit_scale_matches_the_exact_rational_unit() {
        let places = FinitePlaceSet::new(vec![2, 3]).unwrap();
        let q = QsRational::new(45, 8, &places).unwrap();
        let decomposition = q.unit_monoid_decomposition(&places).unwrap();
        let action =
            SemilocalUnitOrbitTransport::from_decomposition(&decomposition, &places).unwrap();

        // u = 2^-3 * 3^2 = 9/8 and m = 5.
        assert_eq!(action.archimedean_absolute_scale().unwrap(), 9.0 / 8.0);
    }
}
