# Analytically certified envelopes for the actual compact Weil matrix

Evidence class: **rigorous elementary analytic inequalities + exact checked
integer arithmetic** for the ideal smooth compact Legendre-weighted
family. Implementation: [src/weil_certified_enclosure.rs](../src/weil_certified_enclosure.rs).

These are the **first source-locked full-entry enclosures constructed by the
bench without assuming that a numerical quadrature was correctly rounded**.
They are deliberately broad, sign-indeterminate enclosures. They do **not**
certify the existing binary64 numerical values, the true near-self-dual Sonin
complement, uniform positivity, or the Riemann hypothesis.

## The precise source object

The actual finite pairing code in
[src/weil_quadratic_matrix.rs](../src/weil_quadratic_matrix.rs) computes the
compact Connes–Consani Weil decomposition

    A_ij = pole_ij - W_R(theta_ij) - sum_(p^k) W_(p^k)(theta_ij),

    theta_ij(s) = integral_R h_i(e^u)*h_j(e^(u+s)) du,

with h_j=Qg_j, Q=-d_x^2+1/4 for x=log(rho), and

    g_j(rho) = exp(-1/(t*(1-t))) * P_j(2*t-1),
    t=(rho-a)/(b-a),  a<rho<b,

extended by exact zero outside (a,b). The support endpoints a,b are
the stored positive rational numbers. The finite mixed pole terms are
**exactly zero**, by the Mellin integration-by-parts argument in
[EXACT_WEIL_STRUCTURAL_TERMS.md](EXACT_WEIL_STRUCTURAL_TERMS.md).

The source archimedean coefficient and integral are exactly the
mathematical expression implemented in mixed_archimedean_term:

    W_R(theta) = theta(0)*c + integral_0^log(R) I(s) ds,
    c = gamma + log(4*pi*(R-1)/(R+1)), R=b/a,

    I(s) = ((exp(s/2)-1)*2*theta(0)
          +exp(s/2)*(theta(s)+theta(-s)-2*theta(0)))
          /(2*sinh(s)).

Connes & Consani, *Weil positivity and Trace formula, the archimedean
place*, arXiv:2006.13771, Lemma 3.3 for Qg boundary conditions and
the compact explicit formula: https://arxiv.org/abs/2006.13771 .

## Proven global value and derivative bounds

Set d=t(1-t), y=1/d, u=1-2t and B(t)=exp(-y).
Since exp(y)>=y^k/k! for y>0, the exact inequalities
exp(-y)*y^k<=k! hold for every nonnegative integer k.
Direct differentiation yields

    B_t    = B*u*y^2
    B_tt   = B*(-2*y^2-2*u^2*y^3+u^2*y^4)
    B_ttt  = B*(12*u*y^3+6*u^3*y^4-6*u*y^4
                -6*u^3*y^5+u^3*y^6).

As |u|<=1, these establish, globally and at the smooth zero
extension, the exact bounds

    |B|<=1,  |B_t|<=2,  |B_tt|<=40,  |B_ttt|<=1800.

The classical Legendre bound |P_j(z)|<=1 on [-1,1] and the
iterated Markov polynomial inequality imply that for n=j^2

    |P_t|<=2*n,  |P_tt|<=4*n^2,  |P_ttt|<=8*n^3.

Leibniz then bounds the derivatives of g_j=B*P_j:

    |g_j|<=1,
    |g_j,t|<=G1=2+2*n,
    |g_j,tt|<=G2=40+8*n+4*n^2,
    |g_j,ttt|<=G3=1800+240*n+24*n^2+8*n^3.

Let A=ceil(b/(b-a)); this is computed using exact rational
cross-products. With d_x=(rho/(b-a))*d_t, the actual basis h_j=Qg_j obeys

    |h_j(e^x)| <= H0_j=1+A*G1+A^2*G2,
    |d_x h_j(e^x)| <= H1_j=2*A*G1+3*A^2*G2+A^3*G3.

Every arithmetic operation uses checked i128, including powers,
products and sums. If a support ratio or polynomial degree exceeds
the representable range, construction fails rather than creating a
falsely narrow interval.

## Full source pairing bound

Let

    T = ceil(R-1) >= log(R),   C=T+1>=R,
    K = 17+ceil((R+1)/(R-1)).

The bounds use only exact rational support ratios. Because
0<gamma<1, 3<pi<4 and |log z|<=max(z,1/z) for z>0,
|c|<=K, without approximating gamma, pi or log.

Writing X=H0_i*H0_j and
Y=H1_i*H0_j+H0_i*H1_j, the convolution and its shift derivative
satisfy

    |theta_ij(s)| <= T*X,
    |d_s(theta_ij(s)+theta_ij(-s))| <= T*Y.

For the derivative, differentiate the smooth compact convolution;
one factor can alternatively be differentiated through integration
by parts, giving two bounds. Adding them bounds the sum of two shift
derivatives.

At every s>0, sinh(s)>=s, exp(s/2)<=R<=C and
exp(s/2)-1 <= C*s/2. The apparent singularity of I(s)
therefore cancels in the **absolute estimate**:

    |I(s)| <= C*T*(X+Y).

After integrating over [0,log R] and adding the coefficient,

    |W_R(theta_ij)| <= K*T*X + C*T^2*(X+Y).

For every eligible prime power n=p^k<R, its source weight is
log(p)/sqrt(n). This is <=log R<=T; the sum of the
two shifted correlations has modulus <=2*T*X. There are at most

    N=max(0, ceil(R)-2)

integers in the range 2<=n<R. (Counting ALL integers is
conservative; only prime powers actually occur.) Thus

    |sum_(p^k<R) W_(p^k)(theta_ij)| <= 2*N*T^2*X.

Since pole_ij=0 exactly, the FULL ideal mixed Weil entry lies in

    -B_ij <= A_ij <= B_ij,
    B_ij = K*T*X + C*T^2*(X+Y) + 2*N*T^2*X.

All B_ij are explicit checked nonnegative integers. Matrix
intervals are symmetric and have common denominator S=1, directly
matching the exact interval sign gate.

The zero-shift correlation is also rigorously enclosed:

    theta_ij(0) in [-T*X,T*X] (i != j),
    theta_ii(0) in [0,T*H0_i^2].

The nonnegative diagonal follows from a real square under the
logarithmic Haar measure. It does not imply positive Weil diagonal.

## What this proves and does not prove

**Proved for the ideal smooth family:** these analytic inequalities
enclose *every exact entry* of the finite compact Weil matrix under
the cited decomposition, irrespective of the binary64 implementation.
They constitute valid, albeit coarse, input to
[weil_rational_gate.rs](../src/weil_rational_gate.rs).

**Not proved:** a tight enclosure around computed Gauss–Legendre values,
an effective quadrature remainder for the nested correlations, the
rounding error of exp/log/sinh/Legendre evaluation, or any positive or
negative Weil sign. The wide symmetric intervals cannot decide a
nonzero direction: the sign gate returns Undetermined on the initial
fixtures. A finite sign, even if someday certified, would not
imply full Weil positivity or RH.

## Reproducibility

    cargo fmt --all -- --check
    cargo clippy --locked --all-targets --all-features -- -D warnings
    cargo test --locked weil_certified_enclosure
    cargo test --locked --all-targets
    cargo run --locked --release

Regression tests cover ratio R=2 (prime contribution identically
zero), strict rational prime-window endpoints, interval symmetry,
nonnegative zero-shift diagonals, exact sign-gate compatibility,
conservative agreement with current numerical matrix observations,
monotonic degree envelopes and fail-closed i128 overflow. The
numerical observations are **tests, not proof of the bounds**.

Next: replace coarse global integral envelopes with mathematically
derived local derivative/variation bounds plus *certified* elementary
function evaluations, then narrow intervals enough to assess a
finite sign. Do not infer quadrature error from differences between
two integration orders.
