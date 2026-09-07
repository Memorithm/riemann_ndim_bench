# Finite-support H1 continuity envelope for the Weil pairing

This note derives the deterministic bound implemented by
`weil_h1_continuity`.  It is a local analytic statement about the finite compact
Riemann--Weil decomposition already implemented in this repository.  It is not
a proof of the complete Weil criterion or the Riemann hypothesis.

## 1. Logarithmic coordinate

Let `u = log(rho)` and let every function be extended by zero outside one common
compact interval `I=[alpha,beta]`.  Put

`L = beta-alpha`, `R = exp(L) = upper/lower > 1`.

For compact `H1(R)` functions define the mixed correlation

`C(f,g)(t) = integral_R f(u) g(u+t) du`.

Consider target/approximant pairs `(f,f_N)` and `(g,g_N)`.  Write

`e_f = ||f-f_N||_2`, `e_g = ||g-g_N||_2`,

`N_f = max(||f||_2,||f_N||_2)`,
`N_g = max(||g||_2,||g_N||_2)`.

Bilinearity and Cauchy--Schwarz give, for every real `t`,

`|Delta C(t)| <= A := e_f N_g + N_f e_g`.

This is already sufficient for the pole and prime-power terms.

## 2. One derivative closes the real-place singularity

Let `D_f=max(||f'||_2,||f_N'||_2)` and likewise for `D_g`, and let
`e_f'=||f'-f_N'||_2`, `e_g'=||g'-g_N'||_2`, where the derivative is with
respect to `u=log(rho)`.

Differentiating the translated right factor gives

`|Delta C'(t)| <= e_f D_g + N_f e_g'`.

Integration by parts instead differentiates the left factor and gives

`|Delta C'(t)| <= e_f' N_g + D_f e_g`.

Therefore the smaller of these two valid estimates is still a valid bound:

`D = min(e_f D_g + N_f e_g', e_f' N_g + D_f e_g)`.

For the symmetrized correlation `C_sym(t)=C(t)+C(-t)`, this yields

`|Delta C_sym(0)| <= 2A`,

`|Delta C_sym(t)-Delta C_sym(0)| <= 2D t`, for `t>=0`.

This is the additional regularity missing from a pure `L2` argument.

## 3. Critical-pole term

The critical Mellin moments are

`M_+(f)=integral_I f(u) exp(u/2) du`,

`M_-(f)=integral_I f(u) exp(-u/2) du`.

Their weight norms satisfy

`||exp(u/2)||_2 ||exp(-u/2)||_2 = sqrt(R)-1/sqrt(R)`.

The mixed pole term contains two products of critical moments.  Hence

`|Delta pole| <= 2 (sqrt(R)-1/sqrt(R)) A`.

The coefficient used in code is therefore

`K_pole(R)=2(sqrt(R)-1/sqrt(R))`.

## 4. Prime-power term

The source term has coefficients `log(p)/sqrt(n)` at prime powers `n=p^k<R`
and uses the two shifts `+log(n)` and `-log(n)`.  Since each correlation
difference is bounded by `A`, and since `log(p)<=log(R)`, we can safely enlarge
the prime-power set to all integers and use

`sum_{n=2}^m 1/sqrt(n) <= 2 sqrt(m) <= 2 sqrt(R)`.

Thus

`|Delta prime| <= 4 sqrt(R) log(R) A`.

This deliberately coarse bound avoids using any finite prime enumeration.  It
is a deterministic envelope, not a sharp constant.

## 5. Real-place term

The source-locked implementation uses

`W_R = 1/2 C_sym(0) c_R
       + integral_0^L [exp(t/2) C_sym(t)-C_sym(0)]/(2 sinh(t)) dt`,

where

`c_R = gamma + log(4*pi*(R-1)/(R+1))`.

The code evaluates the numerator in an algebraically equivalent cancellation-
reduced form.  For the difference between target and approximant pairings,

`|Delta C_sym(0)| <= 2A`

and

`|Delta C_sym(t)-Delta C_sym(0)| <= 2Dt`.

Using, for `0<=t<=L`,

`sinh(t) >= t`,

`exp(t/2)-1 <= (t/2) exp(L/2)`,

`exp(t/2) <= exp(L/2)=sqrt(R)`,

we obtain

`|Delta W_R| <= [|c_R| + (L/2)sqrt(R)] A + L sqrt(R) D`.

Thus

`K_arch,value(R)=|c_R| + (L/2)sqrt(R)`,

`K_arch,derivative(R)=L sqrt(R)`.

The derivative term is exactly why the pure `L2` probe in the preceding work is
not by itself an analytic continuity theorem for the full real-place term.

## 6. Complete finite pairing envelope

The finite source decomposition is

`psi = pole - W_R - prime`.

The triangle inequality gives

`|Delta psi| <= K_value(R) A + K_derivative(R) D`,

where

`K_value = K_pole + 4 sqrt(R) log(R) + K_arch,value`,

`K_derivative = L sqrt(R)`.

`FiniteWeilH1ContinuityBound` retains the three component bounds separately and
also reports their sum.

For the repository's manufactured support `(1/2,7/2)`, `R=7`.  The binary64
regression values currently used by the test suite are approximately

- `K_pole = 4.535573676110727`,
- prime value coefficient `= 20.593577312307954`,
- archimedean value coefficient `= 5.394755003457536`,
- archimedean derivative coefficient `= 5.148394328076988`.

These values instantiate the displayed formulas; they are not experimentally
fitted constants.

## 7. Exact scope of the result

The inequality is mathematically valid when the supplied quantities are genuine
upper bounds for the stated `L2` and logarithmic-derivative norms of compact
`H1` functions.  The Rust implementation evaluates the resulting closed
formula in binary64.  It does not automatically certify norm bounds obtained by
numerical quadrature.

Consequently this result does **not** establish:

- certified bounds for the sine/Legendre approximation errors currently
  measured numerically;
- density or completeness of those finite families;
- continuity on the complete test-function space required by the Weil
  criterion;
- positivity of the full Weil quadratic form;
- the semilocal `L2(X_S)` identification;
- Conjecture 4.1;
- the Riemann hypothesis.

The immediate next engineering/mathematical bridge is to produce validated or
certified upper bounds for the `H1` approximation quantities of concrete test
families and compare the resulting analytic envelope with the directly measured
pairing residuals.
