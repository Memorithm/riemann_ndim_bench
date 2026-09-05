# Weil compact support -> finite places

Status: executable source fact, separated explicitly from Conjecture 4.1.

Primary source:

- Alain Connes, Caterina Consani, *The Scaling Hamiltonian*, Fact 3.2 and Conjecture 4.1, arXiv:1910.14368.

## Fact 3.2

The Riemann--Weil explicit formula applied to a compactly supported test function involves only finitely many places.

At a finite prime `p`, the local contribution samples the test function on non-zero powers of `p`. For a compact support bounded away from both zero and infinity, sufficiently large primes and all of their non-zero powers lie outside the support.

This is a `SOURCE FACT`.

## Multiplicative support bookkeeping

If

```text
supp(f) subset [a,b],
supp(g) subset [c,d],
```

then multiplicative convolution satisfies

```text
supp(f*g) subset [ac, bd].
```

For the involution used in the Weil convolution square,

```text
h1*(rho) = overline(h1(rho^-1)),
```

one has the support envelope

```text
supp(h1*) subset [1/b, 1/a].
```

Hence

```text
supp(h1 * h1*) subset [a/b, b/a].
```

## The q-window

Conjecture 4.1 uses the open support window

```text
h1 support subset (q^-1/2, q^1/2).
```

Pure support arithmetic then gives

```text
h1 * h1* support subset (q^-1, q).
```

Therefore, for every prime `p >= q`, all positive powers satisfy

```text
p^n >= p >= q,
```

and all negative powers satisfy

```text
p^-n <= p^-1 <= q^-1.
```

For an actual compact support strictly inside the open window, those values are outside support. Thus primes `p >= q` are absent from the explicit formula.

This explains the finite source set

```text
S(q) = {infinity} union {p prime : p < q}.
```

The exclusion of the omitted primes is support bookkeeping plus Fact 3.2. It is **not** Conjecture 4.1.

## What Conjecture 4.1 adds

Conjecture 4.1 makes a much stronger statement: the semilocal operator-theoretic framework associated with `S(q)` should suffice to prove the Weil inequality for every test function in the declared support window.

RiemannBench must keep the two statements separate:

```text
compact support => only primes in a finite set can contribute       SOURCE FACT

those places suffice to establish Weil positivity on the window     SOURCE CONJECTURE
```

The first is executable here. The second remains open.

## Repository API

`src/weil_support.rs` provides:

- `WeilSupportWindow` for the open source window parameter `q`;
- source `h1` and convolution-window endpoints;
- strict containment tests for actual compact supports;
- inverse-support and convolution-support envelopes;
- a sufficient exact exclusion test for a finite prime;
- deterministic construction of `{p : p < q}` for integer `q`.

The finite place set is represented by the `FinitePlaceSet` contract introduced with the semilocal trace formula.

## Scientific boundary

No positivity statement is returned by this API.
No finite crossing is identified with a zeta zero.
No version of Conjecture 4.1 is assumed or marked proved.
No RH proof claim is made.
