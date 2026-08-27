use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Q {
    num: i128,
    den: i128,
}

impl Q {
    fn new(num: i128, den: i128) -> Self {
        assert_ne!(den, 0);

        let (num, den) = if den < 0 { (-num, -den) } else { (num, den) };
        let gcd = gcd_i128(num, den);

        Self {
            num: num / gcd,
            den: den / gcd,
        }
    }

    fn integer(value: i128) -> Self {
        Self::new(value, 1)
    }
}

fn gcd_i128(mut left: i128, mut right: i128) -> i128 {
    left = left.abs();
    right = right.abs();

    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }

    left.max(1)
}

impl Add for Q {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.num * rhs.den + rhs.num * self.den, self.den * rhs.den)
    }
}

impl Sub for Q {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl Mul for Q {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(self.num * rhs.num, self.den * rhs.den)
    }
}

impl Div for Q {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        assert_ne!(rhs.num, 0);
        Self::new(self.num * rhs.den, self.den * rhs.num)
    }
}

impl Neg for Q {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.num, self.den)
    }
}

fn q(num: i128, den: i128) -> Q {
    Q::new(num, den)
}

fn poly_sub(left: [Q; 3], right: [Q; 3]) -> [Q; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn poly_scale(poly: [Q; 3], factor: Q) -> [Q; 3] {
    [poly[0] * factor, poly[1] * factor, poly[2] * factor]
}

fn poly_add(left: [Q; 3], right: [Q; 3]) -> [Q; 3] {
    [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
}

fn universal_a_in_i(offset: Q) -> [Q; 3] {
    // A(j)=4j^2+3j+5/16 with j=i+offset.
    [
        Q::integer(4) * offset * offset + Q::integer(3) * offset + q(5, 16),
        Q::integer(8) * offset + Q::integer(3),
        Q::integer(4),
    ]
}

fn universal_b_in_i(offset: Q) -> [Q; 3] {
    // B(j)=4j^2-3j+5/16 with j=i+offset.
    [
        Q::integer(4) * offset * offset - Q::integer(3) * offset + q(5, 16),
        Q::integer(8) * offset - Q::integer(3),
        Q::integer(4),
    ]
}

fn raw_a_in_i(epsilon: i128) -> [Q; 3] {
    // d=2i+epsilon, A=(d+1/2)(d+3/2).
    let e = Q::integer(epsilon);
    [
        (e + q(1, 2)) * (e + q(3, 2)),
        Q::integer(4 * epsilon + 4),
        Q::integer(4),
    ]
}

fn raw_b_in_i(epsilon: i128) -> [Q; 3] {
    // d=2i+epsilon, B=d(d-1).
    let e = Q::integer(epsilon);
    [
        e * (e - Q::integer(1)),
        Q::integer(4 * epsilon - 2),
        Q::integer(4),
    ]
}

fn linear_product(shift_left: Q, shift_right: Q) -> [Q; 3] {
    // (j+shift_left)(j+shift_right).
    [
        shift_left * shift_right,
        shift_left + shift_right,
        Q::integer(1),
    ]
}

fn eval_poly(poly: [Q; 3], value: Q) -> Q {
    poly[0] + poly[1] * value + poly[2] * value * value
}

#[test]
fn exact_raw_to_staggered_rewrite_holds_in_both_parities() {
    for epsilon in [0_i128, 1_i128] {
        // j=i+(4*epsilon+1)/8.
        let offset = q(4 * epsilon + 1, 8);

        assert_eq!(raw_a_in_i(epsilon), universal_a_in_i(offset));
        assert_eq!(raw_b_in_i(epsilon), universal_b_in_i(offset));

        // D_raw=4d+1=8i+4epsilon+1 and D_universal=8j.
        let raw_d = [Q::integer(4 * epsilon + 1), Q::integer(8)];
        let universal_d = [Q::integer(8) * offset, Q::integer(8)];
        assert_eq!(raw_d, universal_d);
    }
}

#[test]
fn exact_first_order_solution_forces_the_second_order_rhs() {
    // Canonical recurrence:
    // A(j)(g[n+1]-g[n]) - B(j)(g[n]-g[n-1]) + mu D(j)g[n] = 0,
    // with A-B=6j, D=8j and u[n]=-(4/3)n.
    let a = [q(5, 16), Q::integer(3), Q::integer(4)];
    let b = [q(5, 16), Q::integer(-3), Q::integer(4)];
    let d = [Q::integer(0), Q::integer(8), Q::integer(0)];

    assert_eq!(
        poly_sub(a, b),
        [Q::integer(0), Q::integer(6), Q::integer(0)]
    );

    let first_order_step = q(-4, 3);
    let mu1_residual = poly_add(poly_scale(poly_sub(a, b), first_order_step), d);
    assert_eq!(mu1_residual, [Q::integer(0); 3]);

    // At order mu^2, moving D(j)u[n] to the right gives
    // -8j * (-(4/3)n) = (32/3) j n.
    let forcing_jn = -Q::integer(8) * q(-4, 3);
    assert_eq!(forcing_jn, q(32, 3));
}

#[test]
fn exact_homogeneous_increment_ratio_factorization_is_preserved() {
    let a = [q(5, 16), Q::integer(3), Q::integer(4)];
    let b = [q(5, 16), Q::integer(-3), Q::integer(4)];

    let a_factored = poly_scale(linear_product(q(1, 8), q(5, 8)), Q::integer(4));
    let b_factored = poly_scale(linear_product(q(-5, 8), q(-1, 8)), Q::integer(4));

    assert_eq!(a, a_factored);
    assert_eq!(b, b_factored);

    // Thus at mu=0 the increment equation A*delta[n+1]-B*delta[n]=0
    // has delta[n+1]/delta[n]=B/A with the exact quarter-offset factors.
}

#[test]
fn exact_left_boundary_selects_the_universal_negative_first_order_sign() {
    let a = [q(5, 16), Q::integer(3), Q::integer(4)];
    let b = [q(5, 16), Q::integer(-3), Q::integer(4)];

    for j0 in [q(1, 8), q(5, 8)] {
        assert_eq!(eval_poly(b, j0), Q::integer(0));

        let a0 = eval_poly(a, j0);
        let d0 = Q::integer(8) * j0;

        // With g0=1 and B(j0)=0, the canonical recurrence gives
        // g1=1-mu*D(j0)/A(j0), and D/A=4/3 in both parities.
        assert_eq!(d0 / a0, q(4, 3));
    }
}
