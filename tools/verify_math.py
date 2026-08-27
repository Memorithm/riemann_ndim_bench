#!/usr/bin/env python3
"""Deterministic mathematical verifier for the Riemann research harness.

The module deliberately accepts only small, whitelisted expression grammars.
It is intended to audit fragile algebra performed by an LLM; it is not a CAS
and it never evaluates model-supplied Python code.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
from math import comb
from typing import Callable

import mpmath as mp


# ---------------------------------------------------------------------------
# Shared exact helpers
# ---------------------------------------------------------------------------


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def parse_fraction(text: str) -> Fraction:
    text = text.strip()
    if not text:
        raise ValueError("empty rational")
    return Fraction(text)


def parse_fraction_list(text: str) -> list[Fraction]:
    values = [part.strip() for part in text.split(",") if part.strip()]
    if not values:
        raise ValueError("expected at least one rational argument")
    return [parse_fraction(value) for value in values]


def safe_rational_expr(text: str) -> Fraction:
    """Evaluate integer/rational arithmetic using a tiny AST grammar."""

    tree = ast.parse(text, mode="eval")

    def walk(node: ast.AST) -> Fraction:
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return Fraction(node.value)
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
                if right.denominator != 1:
                    raise ValueError("rational mode requires integer powers")
                exponent = right.numerator
                if abs(exponent) > 64:
                    raise ValueError("power too large")
                return left**exponent
            raise ValueError("unsupported binary operator")
        raise ValueError(f"unsupported rational syntax: {type(node).__name__}")

    return walk(tree)


def verify_rational(expr: str) -> None:
    value = safe_rational_expr(expr)
    print("mode=rational")
    print(f"expr={expr}")
    print(f"exact_value={fraction_text(value)}")
    print("exact_status=PROVED_BY_RATIONAL_ARITHMETIC")


# ---------------------------------------------------------------------------
# Gamma quotient reduction
# ---------------------------------------------------------------------------


def gamma_reduce_argument(argument: Fraction) -> tuple[Fraction, Fraction | None]:
    """Reduce Gamma(q) by recurrence until q is in (0, 1].

    Returns (rational factor, residual base).  base=None denotes Gamma(1)=1.
    Only positive rational arguments are accepted.
    """

    if argument <= 0:
        raise ValueError("Gamma arguments must be positive rationals")

    factor = Fraction(1)
    current = argument
    while current > 1:
        current -= 1
        factor *= current
    if current == 1:
        return factor, None
    return factor, current


def gamma_product_canonical(arguments: list[Fraction]) -> tuple[Fraction, Counter]:
    factor = Fraction(1)
    bases: Counter[Fraction] = Counter()
    for argument in arguments:
        local_factor, base = gamma_reduce_argument(argument)
        factor *= local_factor
        if base is not None:
            bases[base] += 1
    return factor, bases


def symbolic_power_text(base: str, exponent: Fraction) -> str:
    if exponent == 0:
        return ""
    if exponent == 1:
        return base
    return f"{base}^({fraction_text(exponent)})"


def simplify_quarter_half_bases(
    rational_factor: Fraction,
    bases: Counter,
) -> tuple[Fraction, Fraction, Fraction, Counter]:
    """Use exact identities for the quarter/half Gamma bases.

    Gamma(1/2) = pi^(1/2)
    Gamma(1/4) Gamma(3/4) = pi * 2^(1/2)

    Returns rational factor, exponent of 2, exponent of pi, unresolved bases.
    """

    bases = Counter({key: value for key, value in bases.items() if value})
    e14 = bases.pop(Fraction(1, 4), 0)
    e12 = bases.pop(Fraction(1, 2), 0)
    e34 = bases.pop(Fraction(3, 4), 0)

    if e14 and e34 and (e14 > 0) == (e34 > 0):
        pair = min(abs(e14), abs(e34))
        pair *= 1 if e14 > 0 else -1
        e14 -= pair
        e34 -= pair
        two_power = Fraction(pair, 2)
        pi_power = Fraction(pair) + Fraction(e12, 2)
    else:
        two_power = Fraction(0)
        pi_power = Fraction(e12, 2)

    if e14:
        bases[Fraction(1, 4)] = e14
    if e34:
        bases[Fraction(3, 4)] = e34

    return rational_factor, two_power, pi_power, bases


def gamma_numeric(arguments: list[Fraction]) -> mp.mpf:
    value = mp.mpf(1)
    for argument in arguments:
        value *= mp.gamma(mp.mpf(argument.numerator) / argument.denominator)
    return value


def verify_gamma_quotient(numerator_text: str, denominator_text: str) -> None:
    numerator = parse_fraction_list(numerator_text)
    denominator = parse_fraction_list(denominator_text)

    n_factor, n_bases = gamma_product_canonical(numerator)
    d_factor, d_bases = gamma_product_canonical(denominator)
    factor = n_factor / d_factor
    bases = n_bases.copy()
    bases.subtract(d_bases)
    bases = Counter({key: value for key, value in bases.items() if value})

    factor, two_power, pi_power, unresolved = simplify_quarter_half_bases(
        factor, bases
    )

    mp.mp.dps = 80
    numeric = gamma_numeric(numerator) / gamma_numeric(denominator)

    print("mode=gamma_quotient")
    print(f"numerator={numerator_text}")
    print(f"denominator={denominator_text}")
    print(f"rational_factor={fraction_text(factor)}")
    print(f"two_power={fraction_text(two_power)}")
    print(f"pi_power={fraction_text(pi_power)}")

    symbolic_parts = []
    if factor != 1 or (two_power == 0 and pi_power == 0):
        symbolic_parts.append(fraction_text(factor))
    if two_power:
        symbolic_parts.append(symbolic_power_text("2", two_power))
    if pi_power:
        symbolic_parts.append(symbolic_power_text("pi", pi_power))
    symbolic = " * ".join(part for part in symbolic_parts if part) or "1"
    print(f"symbolic_factor={symbolic}")

    if unresolved:
        rendered = " * ".join(
            f"Gamma({fraction_text(base)})^{exponent}"
            for base, exponent in sorted(unresolved.items())
        )
        print("exact_status=UNRESOLVED_GAMMA_BASES")
        print(f"unresolved={rendered}")
    else:
        print("exact_status=PROVED_BY_GAMMA_RECURRENCE_AND_SPECIAL_IDENTITIES")
        print("unresolved=none")

    print(f"numeric_value={mp.nstr(numeric, 60)}")


# ---------------------------------------------------------------------------
# Safe high-precision numerical identities
# ---------------------------------------------------------------------------


MP_FUNCTIONS: dict[str, Callable] = {
    "gamma": mp.gamma,
    "sqrt": mp.sqrt,
    "sin": mp.sin,
    "cos": mp.cos,
    "exp": mp.exp,
    "log": mp.log,
}


def safe_mp_expr(text: str) -> mp.mpf:
    tree = ast.parse(text, mode="eval")

    def walk(node: ast.AST):
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.Constant):
            if isinstance(node.value, int):
                return mp.mpf(node.value)
            if isinstance(node.value, float):
                return mp.mpf(str(node.value))
            raise ValueError("unsupported numeric literal")
        if isinstance(node, ast.Name):
            if node.id == "pi":
                return mp.pi
            if node.id == "e":
                return mp.e
            raise ValueError(f"unsupported name: {node.id}")
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
                return left**right
            raise ValueError("unsupported binary operator")
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name):
                raise ValueError("only direct whitelisted function calls are allowed")
            function = MP_FUNCTIONS.get(node.func.id)
            if function is None:
                raise ValueError(f"unsupported function: {node.func.id}")
            if node.keywords or len(node.args) != 1:
                raise ValueError("numeric functions take exactly one positional argument")
            return function(walk(node.args[0]))
        raise ValueError(f"unsupported numeric syntax: {type(node).__name__}")

    return walk(tree)


def verify_numeric_identity(lhs_text: str, rhs_text: str, tolerance_text: str) -> None:
    mp.mp.dps = 80
    lhs = safe_mp_expr(lhs_text)
    rhs = safe_mp_expr(rhs_text)
    tolerance = mp.mpf(tolerance_text)
    error = abs(lhs - rhs)
    scale = max(abs(lhs), abs(rhs), mp.mpf(1))
    rel_error = error / scale

    print("mode=numeric_identity")
    print("precision_digits=80")
    print(f"lhs={mp.nstr(lhs, 60)}")
    print(f"rhs={mp.nstr(rhs, 60)}")
    print(f"abs_error={mp.nstr(error, 20)}")
    print(f"rel_error={mp.nstr(rel_error, 20)}")
    print(
        "numeric_status="
        + ("MATCH_WITHIN_TOLERANCE" if error <= tolerance else "MISMATCH")
    )
    print("warning=NUMERIC_MATCH_IS_NOT_SYMBOLIC_PROOF")


# ---------------------------------------------------------------------------
# Asymptotic power comparison
# ---------------------------------------------------------------------------


def parse_mp_list(text: str) -> list[mp.mpf]:
    values = [part.strip() for part in text.split(",") if part.strip()]
    if not values:
        raise ValueError("empty numeric list")
    return [mp.mpf(value) for value in values]


def fit_limit_power(n_values: list[mp.mpf], y_values: list[mp.mpf], p: Fraction):
    exponent = mp.mpf(p.numerator) / p.denominator
    x_values = [n ** (-exponent) for n in n_values]
    count = mp.mpf(len(x_values))
    sx = mp.fsum(x_values)
    sy = mp.fsum(y_values)
    sxx = mp.fsum(x * x for x in x_values)
    sxy = mp.fsum(x * y for x, y in zip(x_values, y_values))
    denom = count * sxx - sx * sx
    if denom == 0:
        raise ValueError("singular asymptotic fit")
    slope = (count * sxy - sx * sy) / denom
    limit = (sy - slope * sx) / count
    sse = mp.fsum(
        (y - (limit + slope * x)) ** 2 for x, y in zip(x_values, y_values)
    )
    return limit, slope, sse


def verify_asymptotic_power(n_text: str, y_text: str, powers_text: str) -> None:
    mp.mp.dps = 80
    n_values = parse_mp_list(n_text)
    y_values = parse_mp_list(y_text)
    powers = parse_fraction_list(powers_text)
    if len(n_values) != len(y_values):
        raise ValueError("n and y lists must have equal length")
    if len(n_values) < 3:
        raise ValueError("at least three observations are required")

    results = []
    print("mode=asymptotic_power")
    for power in powers:
        limit, slope, sse = fit_limit_power(n_values, y_values, power)
        results.append((sse, power, limit, slope))
        print(
            f"candidate_p={fraction_text(power)} "
            f"limit={mp.nstr(limit, 30)} "
            f"slope={mp.nstr(slope, 20)} "
            f"sse={mp.nstr(sse, 20)}"
        )

    results.sort(key=lambda row: row[0])
    best_sse, best_power, best_limit, _ = results[0]
    print(f"best_power={fraction_text(best_power)}")
    print(f"best_limit={mp.nstr(best_limit, 30)}")
    print(f"best_sse={mp.nstr(best_sse, 20)}")

    if len(y_values) >= 3:
        differences = [y_values[i] - y_values[i + 1] for i in range(len(y_values) - 1)]
        ratios = []
        for left, right in zip(differences, differences[1:]):
            if left != 0:
                ratios.append(right / left)
        if ratios:
            print("successive_difference_ratios=" + ",".join(mp.nstr(v, 16) for v in ratios))

    print("warning=ASYMPTOTIC_FIT_IS_NUMERICAL_EVIDENCE_NOT_PROOF")


# ---------------------------------------------------------------------------
# Exact bivariate polynomial engine used by perturbative recurrence
# ---------------------------------------------------------------------------


Poly = dict[tuple[int, int], Fraction]


def poly_clean(p: Poly) -> Poly:
    return {monomial: coeff for monomial, coeff in p.items() if coeff}


def poly_const(value) -> Poly:
    value = Fraction(value)
    return {} if value == 0 else {(0, 0): value}


def poly_var(name: str) -> Poly:
    if name == "j":
        return {(1, 0): Fraction(1)}
    if name == "n":
        return {(0, 1): Fraction(1)}
    raise ValueError(f"unsupported polynomial variable: {name}")


def poly_add(a: Poly, b: Poly) -> Poly:
    out = dict(a)
    for monomial, coeff in b.items():
        out[monomial] = out.get(monomial, Fraction(0)) + coeff
    return poly_clean(out)


def poly_neg(a: Poly) -> Poly:
    return {monomial: -coeff for monomial, coeff in a.items()}


def poly_sub(a: Poly, b: Poly) -> Poly:
    return poly_add(a, poly_neg(b))


def poly_mul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    for (ja, na), ca in a.items():
        for (jb, nb), cb in b.items():
            monomial = (ja + jb, na + nb)
            out[monomial] = out.get(monomial, Fraction(0)) + ca * cb
    return poly_clean(out)


def poly_pow(a: Poly, exponent: int) -> Poly:
    if exponent < 0 or exponent > 12:
        raise ValueError("polynomial exponent must be between 0 and 12")
    out = poly_const(1)
    base = a
    power = exponent
    while power:
        if power & 1:
            out = poly_mul(out, base)
        base = poly_mul(base, base)
        power >>= 1
    return out


def poly_constant_value(p: Poly) -> Fraction:
    p = poly_clean(p)
    if not p:
        return Fraction(0)
    if set(p) != {(0, 0)}:
        raise ValueError("division is allowed only by an exact constant")
    return p[(0, 0)]


def exact_poly_expr(text: str) -> Poly:
    tree = ast.parse(text, mode="eval")

    def walk(node: ast.AST) -> Poly:
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return poly_const(node.value)
        if isinstance(node, ast.Name):
            return poly_var(node.id)
        if isinstance(node, ast.UnaryOp):
            value = walk(node.operand)
            if isinstance(node.op, ast.UAdd):
                return value
            if isinstance(node.op, ast.USub):
                return poly_neg(value)
            raise ValueError("unsupported polynomial unary operator")
        if isinstance(node, ast.BinOp):
            left = walk(node.left)
            right = walk(node.right)
            if isinstance(node.op, ast.Add):
                return poly_add(left, right)
            if isinstance(node.op, ast.Sub):
                return poly_sub(left, right)
            if isinstance(node.op, ast.Mult):
                return poly_mul(left, right)
            if isinstance(node.op, ast.Div):
                denominator = poly_constant_value(right)
                if denominator == 0:
                    raise ZeroDivisionError("division by zero")
                return {monomial: coeff / denominator for monomial, coeff in left.items()}
            if isinstance(node.op, ast.Pow):
                exponent = poly_constant_value(right)
                if exponent.denominator != 1:
                    raise ValueError("polynomial powers require integer exponents")
                return poly_pow(left, exponent.numerator)
            raise ValueError("unsupported polynomial binary operator")
        raise ValueError(f"unsupported polynomial syntax: {type(node).__name__}")

    return poly_clean(walk(tree))


def poly_shift_n(p: Poly, shift: int) -> Poly:
    out: Poly = {}
    for (j_power, n_power), coeff in p.items():
        for k in range(n_power + 1):
            shifted = coeff * Fraction(comb(n_power, k)) * Fraction(shift) ** (n_power - k)
            monomial = (j_power, k)
            out[monomial] = out.get(monomial, Fraction(0)) + shifted
    return poly_clean(out)


def poly_text(p: Poly) -> str:
    p = poly_clean(p)
    if not p:
        return "0"
    ordered = sorted(
        p.items(),
        key=lambda item: (-(item[0][0] + item[0][1]), -item[0][0], -item[0][1]),
    )
    parts: list[str] = []
    for (j_power, n_power), coeff in ordered:
        sign = "-" if coeff < 0 else "+"
        magnitude = abs(coeff)
        factors = []
        if magnitude != 1 or (j_power == 0 and n_power == 0):
            factors.append(fraction_text(magnitude))
        if j_power:
            factors.append("j" if j_power == 1 else f"j^{j_power}")
        if n_power:
            factors.append("n" if n_power == 1 else f"n^{n_power}")
        term = "*".join(factors) or "1"
        if not parts:
            parts.append(("-" if sign == "-" else "") + term)
        else:
            parts.append(f" {sign} {term}")
    return "".join(parts)


def perturbative_recurrence_report(a_text: str, b_text: str, d_text: str, u_text: str) -> None:
    """Expand exactly

    A(j)(g[n+1]-g[n]) - B(j)(g[n]-g[n-1]) + mu D(j) g[n] = 0

    with g[n] = 1 + mu*u[n] + mu^2*v[n] + O(mu^3).
    """

    A = exact_poly_expr(a_text)
    B = exact_poly_expr(b_text)
    D = exact_poly_expr(d_text)
    u = exact_poly_expr(u_text)

    u_plus = poly_shift_n(u, +1)
    u_minus = poly_shift_n(u, -1)
    delta_u_plus = poly_sub(u_plus, u)
    delta_u_minus = poly_sub(u, u_minus)
    order1 = poly_add(poly_sub(poly_mul(A, delta_u_plus), poly_mul(B, delta_u_minus)), D)
    mu2_known = poly_mul(D, u)
    mu2_rhs = poly_neg(mu2_known)

    print("mode=perturbative_recurrence")
    print("recurrence=A(j)*(g[n+1]-g[n])-B(j)*(g[n]-g[n-1])+mu*D(j)*g[n]=0")
    print("ansatz=g[n]=1+mu*u[n]+mu^2*v[n]+O(mu^3)")
    print("A=" + poly_text(A))
    print("B=" + poly_text(B))
    print("D=" + poly_text(D))
    print("u=" + poly_text(u))
    print("mu0_equation=0")
    print("mu1_residual=" + poly_text(order1))
    print("mu1_status=" + ("PROVED_EXACT_SOLUTION" if not order1 else "CANDIDATE_U_FAILS"))
    print(
        "mu2_equation=A(j)*(v[n+1]-v[n])-B(j)*(v[n]-v[n-1])"
        + ("+" if not poly_text(mu2_known).startswith("-") else "")
        + poly_text(mu2_known)
        + "=0"
    )
    print(
        "mu2_increment_form=A(j)*delta_v[n+1]-B(j)*delta_v[n]="
        + poly_text(mu2_rhs)
    )
    print("mu2_forcing_rhs=" + poly_text(mu2_rhs))
    print("exact_status=PROVED_BY_FORMAL_COEFFICIENT_EXTRACTION")


# ---------------------------------------------------------------------------
# Exact affine recurrence-index transform
# ---------------------------------------------------------------------------


UPoly = dict[int, Fraction]


def upoly_clean(p: UPoly) -> UPoly:
    return {power: coeff for power, coeff in p.items() if coeff}


def upoly_add(a: UPoly, b: UPoly) -> UPoly:
    out = dict(a)
    for power, coeff in b.items():
        out[power] = out.get(power, Fraction(0)) + coeff
    return upoly_clean(out)


def upoly_mul(a: UPoly, b: UPoly) -> UPoly:
    out: UPoly = {}
    for pa, ca in a.items():
        for pb, cb in b.items():
            out[pa + pb] = out.get(pa + pb, Fraction(0)) + ca * cb
    return upoly_clean(out)


def upoly_pow(a: UPoly, exponent: int) -> UPoly:
    if exponent < 0 or exponent > 12:
        raise ValueError("univariate polynomial exponent must be between 0 and 12")
    out: UPoly = {0: Fraction(1)}
    for _ in range(exponent):
        out = upoly_mul(out, a)
    return out


def parse_upoly(text: str, variable: str) -> UPoly:
    tree = ast.parse(text, mode="eval")

    def constant_value(p: UPoly) -> Fraction:
        p = upoly_clean(p)
        if not p:
            return Fraction(0)
        if set(p) != {0}:
            raise ValueError("division/exponent requires an exact constant")
        return p[0]

    def walk(node: ast.AST) -> UPoly:
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return {} if node.value == 0 else {0: Fraction(node.value)}
        if isinstance(node, ast.Name) and node.id == variable:
            return {1: Fraction(1)}
        if isinstance(node, ast.UnaryOp):
            value = walk(node.operand)
            if isinstance(node.op, ast.UAdd):
                return value
            if isinstance(node.op, ast.USub):
                return {power: -coeff for power, coeff in value.items()}
            raise ValueError("unsupported affine unary operator")
        if isinstance(node, ast.BinOp):
            left = walk(node.left)
            right = walk(node.right)
            if isinstance(node.op, ast.Add):
                return upoly_add(left, right)
            if isinstance(node.op, ast.Sub):
                return upoly_add(left, {power: -coeff for power, coeff in right.items()})
            if isinstance(node.op, ast.Mult):
                return upoly_mul(left, right)
            if isinstance(node.op, ast.Div):
                denominator = constant_value(right)
                if denominator == 0:
                    raise ZeroDivisionError("division by zero")
                return {power: coeff / denominator for power, coeff in left.items()}
            if isinstance(node.op, ast.Pow):
                exponent = constant_value(right)
                if exponent.denominator != 1:
                    raise ValueError("power must be an integer")
                return upoly_pow(left, exponent.numerator)
            raise ValueError("unsupported affine binary operator")
        raise ValueError(f"unsupported affine syntax: {type(node).__name__}")

    return upoly_clean(walk(tree))


def upoly_compose(source: UPoly, substitution: UPoly) -> UPoly:
    out: UPoly = {}
    for power, coeff in source.items():
        term = {k: coeff * v for k, v in upoly_pow(substitution, power).items()}
        out = upoly_add(out, term)
    return upoly_clean(out)


def upoly_text(p: UPoly, variable: str) -> str:
    p = upoly_clean(p)
    if not p:
        return "0"
    parts: list[str] = []
    for power in sorted(p, reverse=True):
        coeff = p[power]
        sign = "-" if coeff < 0 else "+"
        magnitude = abs(coeff)
        factors = []
        if magnitude != 1 or power == 0:
            factors.append(fraction_text(magnitude))
        if power:
            factors.append(variable if power == 1 else f"{variable}^{power}")
        term = "*".join(factors) or "1"
        if not parts:
            parts.append(("-" if sign == "-" else "") + term)
        else:
            parts.append(f" {sign} {term}")
    return "".join(parts)


def recurrence_transform_report(
    a_text: str,
    b_text: str,
    d_text: str,
    raw_var: str,
    site_var: str,
    raw_in_site_text: str,
    source_orientation: str,
) -> None:
    if raw_var == site_var:
        raise ValueError("raw_var and site_var must be distinct")
    if not raw_var.isidentifier() or not site_var.isidentifier():
        raise ValueError("variable names must be identifiers")

    raw_A = parse_upoly(a_text, raw_var)
    raw_B = parse_upoly(b_text, raw_var)
    raw_D = parse_upoly(d_text, raw_var)
    substitution = parse_upoly(raw_in_site_text, site_var)
    if max(substitution, default=0) > 1:
        raise ValueError("raw_in_site must be affine")

    A = upoly_compose(raw_A, substitution)
    B = upoly_compose(raw_B, substitution)
    D = upoly_compose(raw_D, substitution)

    allowed = {"current_minus_neighbors_equals_mu", "verifier_canonical"}
    if source_orientation not in allowed:
        raise ValueError("unsupported source orientation")

    print("mode=recurrence_transform")
    print(f"raw_variable={raw_var}")
    print(f"site_variable={site_var}")
    print(f"raw_in_site={raw_in_site_text}")
    print(f"source_orientation={source_orientation}")
    print(f"A_site={upoly_text(A, site_var)}")
    print(f"B_site={upoly_text(B, site_var)}")
    print(f"D_site={upoly_text(D, site_var)}")
    if source_orientation == "current_minus_neighbors_equals_mu":
        print(
            "canonical_rewrite=A(site)*(g[n+1]-g[n])-B(site)*(g[n]-g[n-1])"
            "+mu*D(site)*g[n]=0"
        )
        print("sign_normalization=MULTIPLY_SOURCE_EQUATION_BY_MINUS_ONE_AFTER_MOVING_RHS_LEFT")
    else:
        print("canonical_rewrite=SOURCE_ALREADY_IN_VERIFIER_CONVENTION")
        print("sign_normalization=NONE")
    print("exact_status=PROVED_BY_AFFINE_SUBSTITUTION_AND_POLYNOMIAL_ARITHMETIC")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Deterministic research math verifier")
    sub = parser.add_subparsers(dest="mode", required=True)

    rational = sub.add_parser("rational")
    rational.add_argument("--expr", required=True)

    gamma = sub.add_parser("gamma-quotient")
    gamma.add_argument("--numerator", required=True)
    gamma.add_argument("--denominator", required=True)

    identity = sub.add_parser("numeric-identity")
    identity.add_argument("--lhs", required=True)
    identity.add_argument("--rhs", required=True)
    identity.add_argument("--tolerance", default="1e-60")

    asymptotic = sub.add_parser("asymptotic-power")
    asymptotic.add_argument("--n", required=True)
    asymptotic.add_argument("--y", required=True)
    asymptotic.add_argument("--powers", default="1/2,1,3/2,2")

    perturb = sub.add_parser("perturbative-recurrence")
    perturb.add_argument("--A", required=True)
    perturb.add_argument("--B", required=True)
    perturb.add_argument("--D", required=True)
    perturb.add_argument("--u", required=True)

    transform = sub.add_parser("recurrence-transform")
    transform.add_argument("--A", required=True)
    transform.add_argument("--B", required=True)
    transform.add_argument("--D", required=True)
    transform.add_argument("--raw-var", required=True)
    transform.add_argument("--site-var", required=True)
    transform.add_argument("--raw-in-site", required=True)
    transform.add_argument(
        "--source-orientation",
        choices=["current_minus_neighbors_equals_mu", "verifier_canonical"],
        required=True,
    )

    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.mode == "rational":
        verify_rational(args.expr)
    elif args.mode == "gamma-quotient":
        verify_gamma_quotient(args.numerator, args.denominator)
    elif args.mode == "numeric-identity":
        verify_numeric_identity(args.lhs, args.rhs, args.tolerance)
    elif args.mode == "asymptotic-power":
        verify_asymptotic_power(args.n, args.y, args.powers)
    elif args.mode == "perturbative-recurrence":
        perturbative_recurrence_report(args.A, args.B, args.D, args.u)
    elif args.mode == "recurrence-transform":
        recurrence_transform_report(
            args.A,
            args.B,
            args.D,
            args.raw_var,
            args.site_var,
            args.raw_in_site,
            args.source_orientation,
        )
    else:  # pragma: no cover - argparse makes this unreachable.
        raise RuntimeError(f"unhandled mode: {args.mode}")


if __name__ == "__main__":
    main()
