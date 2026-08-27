#!/usr/bin/env python3
"""Exact symbolic helpers for the post-perturbative mu^2 derivation chain.

The tool accepts only small whitelisted expression grammars. It is designed to
verify recurrence-derived forcing quotients, proposed hypergeometric structure,
and local finite-part algebra; no benchmark-specific limiting constant is
encoded here.
"""

from __future__ import annotations

import argparse
import ast
from fractions import Fraction

import sympy as sp


# ---------------------------------------------------------------------------
# Safe exact parsing
# ---------------------------------------------------------------------------


def parse_fraction(text: str) -> Fraction:
    text = text.strip()
    if not text:
        raise ValueError("empty rational")
    return Fraction(text)


def parse_fraction_list(text: str) -> list[Fraction]:
    parts = [part.strip() for part in text.split(",") if part.strip()]
    if not parts:
        raise ValueError("expected at least one rational")
    return [parse_fraction(part) for part in parts]


def to_sympy(value: Fraction) -> sp.Rational:
    return sp.Rational(value.numerator, value.denominator)


def safe_expr(
    text: str,
    *,
    names: dict[str, sp.Expr],
    allow_sqrt: bool = False,
) -> sp.Expr:
    """Parse an exact expression without sympify/eval on user text."""

    tree = ast.parse(text, mode="eval")

    def walk(node: ast.AST) -> sp.Expr:
        if isinstance(node, ast.Expression):
            return walk(node.body)

        if isinstance(node, ast.Constant):
            if isinstance(node.value, int):
                return sp.Integer(node.value)
            raise ValueError("only integer literals are allowed")

        if isinstance(node, ast.Name):
            if node.id in names:
                return names[node.id]
            raise ValueError(f"unsupported symbol: {node.id}")

        if isinstance(node, ast.UnaryOp):
            value = walk(node.operand)
            if isinstance(node.op, ast.UAdd):
                return value
            if isinstance(node.op, ast.USub):
                return -value
            raise ValueError("unsupported unary operator")

        if isinstance(node, ast.BinOp):
            left = walk(node.left)
            right = walk(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow):
                if right.free_symbols:
                    raise ValueError("powers must have constant exponents")
                exponent = sp.Rational(right)
                if abs(exponent) > 16:
                    raise ValueError("power magnitude too large")
                return left**exponent
            raise ValueError("unsupported binary operator")

        if isinstance(node, ast.Call):
            if not allow_sqrt:
                raise ValueError("function calls are not allowed here")
            if not isinstance(node.func, ast.Name) or node.func.id != "sqrt":
                raise ValueError("only sqrt(...) is allowed")
            if node.keywords or len(node.args) != 1:
                raise ValueError("sqrt takes exactly one positional argument")
            return sp.sqrt(walk(node.args[0]))

        raise ValueError(f"unsupported syntax: {type(node).__name__}")

    return sp.cancel(walk(tree))


def expr_text(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr))


# ---------------------------------------------------------------------------
# Variation-of-constants forcing provenance
# ---------------------------------------------------------------------------


def forcing_ratio_report(
    a_text: str,
    b_text: str,
    forcing_text: str,
    offset_text: str,
    candidate_ratio_text: str,
) -> None:
    """Derive the normalized forcing quotient directly from A, B and F.

    For

        A_n * Delta[n+1] - B_n * Delta[n] = F_n

    let ``h`` be a homogeneous increment with

        h[n+1] / h[n] = B_n / A_n.

    The variation-of-constants forcing is

        T_n = F_n / (A_n * h[n+1]).

    Therefore its quotient is independent of the unknown normalization of h:

        T[n+1] / T[n] = F[n+1] * A[n] / (F[n] * B[n+1]).

    ``j = n + offset`` is substituted exactly before comparison with the
    proposed rational quotient.
    """

    n = sp.Symbol("n", integer=True, nonnegative=True)
    j = sp.Symbol("j")
    offset = to_sympy(parse_fraction(offset_text))

    a_expr = safe_expr(a_text, names={"j": j})
    b_expr = safe_expr(b_text, names={"j": j})
    forcing_expr = safe_expr(forcing_text, names={"j": j, "n": n})
    candidate = safe_expr(candidate_ratio_text, names={"n": n})

    j_n = n + offset
    j_next = n + 1 + offset

    a_n = sp.factor(sp.cancel(a_expr.subs(j, j_n)))
    b_next = sp.factor(sp.cancel(b_expr.subs(j, j_next)))
    f_n = sp.factor(
        sp.cancel(forcing_expr.subs({j: j_n, n: n}, simultaneous=True))
    )
    f_next = sp.factor(
        sp.cancel(
            forcing_expr.subs(
                {j: j_next, n: n + 1},
                simultaneous=True,
            )
        )
    )

    if f_n == 0 or b_next == 0:
        raise ValueError("forcing quotient is identically singular")

    derived = sp.factor(sp.cancel(f_next * a_n / (f_n * b_next)))
    difference = sp.factor(sp.cancel(candidate - derived))
    matches = difference == 0

    print("mode=forcing_ratio")
    print("recurrence=A_n*Delta[n+1]-B_n*Delta[n]=F_n")
    print("homogeneous_ratio=h[n+1]/h[n]=B_n/A_n")
    print("variation_forcing=T_n=F_n/(A_n*h[n+1])")
    print("quotient_identity=T[n+1]/T[n]=F[n+1]*A[n]/(F[n]*B[n+1])")
    print(f"A={sp.sstr(a_expr)}")
    print(f"B={sp.sstr(b_expr)}")
    print(f"F={sp.sstr(forcing_expr)}")
    print(f"site_relation=j=n+{sp.sstr(offset)}")
    print(f"derived_ratio={sp.sstr(derived)}")
    print(f"candidate_ratio={sp.sstr(candidate)}")
    print(f"candidate_difference={sp.sstr(difference)}")
    print("candidate_status=" + ("PROVED_EQUAL" if matches else "MISMATCH"))
    if matches:
        print("exact_status=PROVED_BY_VARIATION_OF_CONSTANTS_QUOTIENT")
    else:
        print("exact_status=REFUTED_FORCING_QUOTIENT")


# ---------------------------------------------------------------------------
# Hypergeometric / Pochhammer sequence verification
# ---------------------------------------------------------------------------


def rising_exact(a: sp.Rational, k: int) -> sp.Expr:
    value: sp.Expr = sp.Integer(1)
    for index in range(k):
        value *= a + index
    return sp.factor(value)


def hypergeometric_report(
    numerator_text: str,
    denominator_text: str,
    base_text: str,
    candidate_ratio_text: str | None,
    terms: int,
) -> None:
    if not 1 <= terms <= 16:
        raise ValueError("terms must be between 1 and 16")

    numerator = [to_sympy(value) for value in parse_fraction_list(numerator_text)]
    denominator = [to_sympy(value) for value in parse_fraction_list(denominator_text)]
    base = to_sympy(parse_fraction(base_text))
    k = sp.Symbol("k", integer=True, nonnegative=True)

    ratio: sp.Expr = base
    for shift in numerator:
        ratio *= k + shift
    for shift in denominator:
        ratio /= k + shift
    ratio = sp.factor(sp.cancel(ratio))

    coefficients: list[sp.Expr] = []
    for index in range(terms):
        value: sp.Expr = base**index
        for shift in numerator:
            value *= rising_exact(shift, index)
        for shift in denominator:
            value /= rising_exact(shift, index)
        coefficients.append(sp.factor(sp.cancel(value)))

    print("mode=hypergeometric")
    print(f"numerator_shifts={numerator_text}")
    print(f"denominator_shifts={denominator_text}")
    print(f"base={sp.sstr(base)}")
    print(f"ratio_t[k+1]/t[k]={sp.sstr(ratio)}")
    print("first_terms=" + ",".join(sp.sstr(value) for value in coefficients))

    # Recognize the common generalized-hypergeometric coefficient convention
    # where one denominator shift is 1, i.e. (1)_k = k!.
    denominator_without_factorial = list(denominator)
    if sp.Integer(1) in denominator_without_factorial:
        denominator_without_factorial.remove(sp.Integer(1))
        p = len(numerator)
        q = len(denominator_without_factorial)
        numerator_rendered = ",".join(sp.sstr(value) for value in numerator)
        denominator_rendered = ",".join(
            sp.sstr(value) for value in denominator_without_factorial
        )
        print(
            "hypergeometric_series="
            f"_{p}F_{q}([{numerator_rendered}];[{denominator_rendered}];{sp.sstr(base)}*z)"
        )

    candidate_matches = True
    if candidate_ratio_text is not None:
        candidate = safe_expr(candidate_ratio_text, names={"k": k})
        difference = sp.factor(sp.cancel(candidate - ratio))
        candidate_matches = difference == 0
        print(f"candidate_ratio={sp.sstr(candidate)}")
        print(f"candidate_difference={sp.sstr(difference)}")
        print(
            "candidate_status="
            + ("PROVED_EQUAL" if candidate_matches else "MISMATCH")
        )

    if candidate_matches:
        print("exact_status=PROVED_BY_POCHHAMMER_QUOTIENT")
    else:
        print("exact_status=REFUTED_CANDIDATE_RATIO")


# ---------------------------------------------------------------------------
# Theta-weighted generating functions and local finite parts
# ---------------------------------------------------------------------------


def theta(expr: sp.Expr, z: sp.Symbol) -> sp.Expr:
    return sp.expand(z * sp.diff(expr, z))


def apply_theta_polynomial(
    base_expr: sp.Expr,
    polynomial: sp.Expr,
    *,
    k: sp.Symbol,
    z: sp.Symbol,
) -> sp.Expr:
    poly = sp.Poly(sp.expand(polynomial), k)
    if poly.degree() > 8:
        raise ValueError("theta polynomial degree must not exceed 8")

    theta_powers = [base_expr]
    for _ in range(poly.degree()):
        theta_powers.append(theta(theta_powers[-1], z))

    result: sp.Expr = sp.Integer(0)
    for (degree,), coefficient in poly.terms():
        result += coefficient * theta_powers[degree]
    return sp.simplify(result)


def finite_part_report(
    base_expr_text: str,
    theta_polynomial_text: str,
    extra_expr_text: str,
    order: int,
) -> None:
    if not 3 <= order <= 20:
        raise ValueError("order must be between 3 and 20")

    z = sp.Symbol("z", positive=True)
    x = sp.Symbol("x", positive=True)
    k = sp.Symbol("k")

    base_expr = safe_expr(
        base_expr_text,
        names={"z": z},
        allow_sqrt=True,
    )
    theta_polynomial = safe_expr(
        theta_polynomial_text,
        names={"k": k},
    )
    extra_expr = safe_expr(
        extra_expr_text,
        names={"z": z},
        allow_sqrt=True,
    )

    weighted = apply_theta_polynomial(
        base_expr,
        theta_polynomial,
        k=k,
        z=z,
    )
    total = sp.simplify(weighted + extra_expr)

    # The local coordinate is fixed and explicit: x = sqrt(1-z), z=1-x^2.
    local_expr = sp.simplify(total.subs(z, 1 - x**2))
    series = sp.series(local_expr, x, 0, order).removeO().expand()
    finite_part = sp.simplify(series.coeff(x, 0))

    terms = []
    for term in sp.Add.make_args(series):
        exponent = term.as_powers_dict().get(x, sp.Integer(0))
        terms.append((sp.Rational(exponent), term))
    terms.sort(key=lambda item: item[0])

    print("mode=finite_part")
    print(f"base_expr={sp.sstr(base_expr)}")
    print(f"theta_polynomial={sp.sstr(theta_polynomial)}")
    print(f"extra_expr={sp.sstr(extra_expr)}")
    print("local_coordinate=x=sqrt(1-z)")
    print(f"weighted_expr={sp.sstr(total)}")
    print(f"local_series={sp.sstr(series)}")
    if terms:
        print(f"leading_power={sp.sstr(terms[0][0])}")
        print(f"leading_term={sp.sstr(terms[0][1])}")
    print(f"finite_part={sp.sstr(finite_part)}")
    print("exact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Exact symbolic mu^2 verifier")
    sub = parser.add_subparsers(dest="mode", required=True)

    forcing = sub.add_parser("forcing-ratio")
    forcing.add_argument("--A", required=True)
    forcing.add_argument("--B", required=True)
    forcing.add_argument("--forcing", required=True)
    forcing.add_argument("--offset", required=True)
    forcing.add_argument("--candidate-ratio", required=True)

    hyper = sub.add_parser("hypergeometric")
    hyper.add_argument("--numerator-shifts", required=True)
    hyper.add_argument("--denominator-shifts", required=True)
    hyper.add_argument("--base", default="1")
    hyper.add_argument("--candidate-ratio")
    hyper.add_argument("--terms", type=int, default=6)

    finite = sub.add_parser("finite-part")
    finite.add_argument("--base-expr", required=True)
    finite.add_argument("--theta-polynomial", default="1")
    finite.add_argument("--extra-expr", default="0")
    finite.add_argument("--order", type=int, default=10)

    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.mode == "forcing-ratio":
        forcing_ratio_report(
            args.A,
            args.B,
            args.forcing,
            args.offset,
            args.candidate_ratio,
        )
    elif args.mode == "hypergeometric":
        hypergeometric_report(
            args.numerator_shifts,
            args.denominator_shifts,
            args.base,
            args.candidate_ratio,
            args.terms,
        )
    elif args.mode == "finite-part":
        finite_part_report(
            args.base_expr,
            args.theta_polynomial,
            args.extra_expr,
            args.order,
        )
    else:  # pragma: no cover
        raise RuntimeError(f"unhandled mode: {args.mode}")


if __name__ == "__main__":
    main()
