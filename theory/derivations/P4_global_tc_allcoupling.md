# P4 — Two-sided $T_c(\lambda)$ for **all** $\lambda>0$ (Einstein spectrum, $\mu^\*=0$)

**Status: the weak-coupling upper bracket is closed.** Deliverable ③
(`P4_global_tc_certified.md`) left one gap: for $\lambda\lesssim0.582$ its *upper*
$T_c$ side was the **asymptotic** form $h\sim L+c_0+b_1/L+\mathrm{REM}/L^2$ with an
empirically calibrated pad $\mathrm{REM}=0.75$ — not a certificate. This note
replaces that with two upper bounds that are **proven modulo one calibrated scalar
envelope** $C(a)$ (Prop. A′ and Prop. B) and one **certificate-level** one
(Prop. U), so that **every** $\lambda>0$ now carries a two-sided bracket in which
neither side is asymptotic.

**Precise rigor claim (revised 2026-09-15 after review — the earlier wording said
"both sides proven or finite-$N$ certified" without qualification, which
overstated the two $C$-dependent bands):**

| band | upper side rests on | strongest honest label |
|---|---|---|
| $\varpi\le0.14$ ($\lambda\gtrsim34$) | Weyl $\mathfrak g(2)$ envelope | **fully analytic** |
| $(0.14,3]$ | Rayleigh + Kato–Temple | **certificate** (finite-$N$, interval) |
| $(3,25]$ | Prop. L + Prop. U | **certificate** (finite-$N$) |
| $>25$ ($\lambda\lesssim0.2$) | Prop. B (+A′) with $C(a)$ | **proven modulo calibrated $C$** |

So the honest headline is: *no band is asymptotic any more* (that gap is genuinely
closed), and every band is at least certificate-level, but the deep band's upper
side is not yet fully analytic — it inherits three measured constants. Closing
that is the single identified missing ingredient (§7 item 3).

Deliverable script: `theory/derivations/p4_global_tc_allcoupling.py`
(run `OMP_NUM_THREADS=1 python p4_global_tc_allcoupling.py [quick] [K=8]`).
It imports only `H_mat` from `p4_interchange_uniformity` and `node_cert` from
`p4_global_tc_certified`. **No existing file is modified.**

---

## 0. Two soundness defects in the inherited weak bracket

Both were found while trying to re-use `secular_bounds_K` from ③, and both are
reported here rather than silently patched.

**(D1) The inherited Theorem-W upper bound is not a certificate.**
`secular_bounds_K(w, N, K)` evaluates the Neumann partial sum $S_K$ at a **finite**
matrix size $N$. Every term of $S_K$ is a sum of **nonnegative** quantities
($\tau\ge0$ and $\mathcal R\ge0$ entrywise), so
$$S_K^{(N)}\;\le\;S_K^{(\infty)} .$$
Adding the Theorem-W tail to a quantity that is *itself already too small* does not
produce an upper bound on $F(h)$. Hence $\mathrm{root}(S_K^{(N)}+\mathrm{tail})$ is
**not** a certified upper bound on $h$. Numerically the defect is small
($\sim\!10^{-4}$ relative at $N=100$, $\sim\!10^{-9}$ at $N=800$) and the published
$K{=}8$ cross-check numbers are not thereby wrong — but the *logic* has a hole, and
a hole is not a certificate. **The lower bound $\mathrm{root}(S_K^{(N)})$ is
unaffected and is sound** (Prop. L below supplies the missing justification).

**(D2)** The deep-tail upper side was asymptotic, as ③ itself states.

**A third defect was found and fixed in this file's own first draft**, and is worth
recording because it is an easy trap: if $A\preceq D_1$ and $A\preceq D_2$ for
diagonal $D_1,D_2$, it does **not** follow that $A\preceq\min(D_1,D_2)$
elementwise. Counterexample: $A=\mathbf 1_{2\times2}\preceq\mathrm{diag}(2,2)$ and
$\preceq\mathrm{diag}(1.5,4)$, but $\mathrm{diag}(1.5,2)-A$ has eigenvalue $-0.28$.
Competing *matrix* dominations may not be min'd; competing *scalar* bounds on the
same scalar may (this is what `c_shape` does).

---

## 1. Setup

With $\varpi=\omega_E/2\pi T$, $D:=1/\varpi$, $s_n=(n+\tfrac12)D$,
$w_n=1/(2n+1)$, $G(y)=1/(1+y^2)$, $\tau_n=\sqrt{w_n}\,G(s_n)$, the (M1)
Lorentz-addition identity gives
$$\mathcal H(\varpi)=2\tau\tau^{\!\top}+\mathcal R-\mathcal D,$$
$\mathcal R\ge0$ entrywise with $\|\mathcal R\|\le\Theta^\*=0.6349287917881079$
(certified, Lemma M), and $\mathcal D=\mathrm{diag}(w_n+\ell_G(n))$ an **exact finite
sum** (no truncation). Writing $M(h):=h+\mathcal D-\mathcal R$, the secular function
$$F(h):=2\,\tau^{\!\top}M(h)^{-1}\tau$$
is strictly decreasing, and $h(\varpi)=\lambda_1(\mathcal H)$ is the unique root of
$F=1$ above $\lambda_{\max}(M)$ (Theorem W(a)). $h$ is **proven strictly increasing**
in $\varpi$ (Cor (i) via Lemma M), so inverting a two-sided envelope of $h$ yields a
two-sided bracket of $T_c/\omega_E=1/(2\pi\varpi_\*)$ through $\lambda h(\varpi_\*)=1$.

---

## 2. The four bounds

### Prop. L (lower; supplies the justification (D1) was missing)
For the head projection $P=\{n<N\}$, the Schur complement satisfies
$(M^{-1})_{PP}=(M_{PP}-M_{PQ}M_{QQ}^{-1}M_{QP})^{-1}\succeq (M_{PP})^{-1}$, and
$M^{-1}\ge0$ entrywise, so
$$F(h)\;\ge\;2\,\tau_P^{\!\top}(M_{PP})^{-1}\tau_P\;\ge\;\sum_{j=0}^{K}S_j^P(h).$$
Hence $\mathrm{root}\big(\sum_{j\le K}S_j^P\big)\le h(\varpi)$. **Certificate-level**
(finite-$N$, float64 head sums + a rounding pad).

### Prop. U (upper, certificate-level; replaces the unsound $S_K^{(N)}+\mathrm{tail}$)
With $\varepsilon:=\kappa_N^2/(h-\Theta^\*)$, $\kappa_N\ge\|\mathcal R_{PQ}\|$, and
$B_Q:=\sum_{n\ge N}\tau_n^2$,
$$\tfrac12 F(h)\;\le\;\Big[\sqrt{\mathcal I(h-\varepsilon)}+\sqrt{B_Q/(h-\Theta^\*)}\Big]^2,$$
$$\mathcal I(x):=\tau_P^{\!\top}(x+\mathcal D_P-\mathcal R_{PP})^{-1}\tau_P
\;\le\;\tfrac12\Big[\textstyle\sum_{j\le K}S_j^P(x)+S_0^P(x)\,r^{K+1}/(1-r)\Big],\quad r=\Theta^\*/x,$$
using $M_{QQ}\succeq h-\Theta^\*$, $M_{PP}-M_{PQ}M_{QQ}^{-1}M_{QP}\succeq M_{PP}-\varepsilon I$,
and $S_j^P\le r^jS_0^P$ (conjugate by $(x+\mathcal D_P)^{-1/2}$; a compression cannot
increase the norm, so $\|\mathcal R_{PP}\|\le\Theta^\*$). **Every quantity is a
finite-$N$ head sum, so (D1) cannot arise.** Valid only where $h>\Theta^\*$, i.e.
$\varpi\gtrsim0.93$.

### Prop. A′ (upper, **proven modulo the calibrated envelope $C$**; Schur-weighted diagonal domination)
*Label discipline (revised 2026-09-15 after review):* the Schur **inequality** is
proven unconditionally; what is only calibrated is the scalar envelope $C(a)$ that
this bound is evaluated with. So the *proposition* is proven and the *number it
returns* is proven-modulo-$C$. Anywhere this bound is used the rigor column reads
`proven-cal`, never a bare `proven`.
Schur test with weight $\phi_m=(m+\tfrac12)^{-1/2}$: for $\mathcal R\ge0$ symmetric,
$\mathcal R\preceq\mathrm{diag}(c)$ with $c_n=(\mathcal R\phi)_n/\phi_n$. The weight
**collapses exactly**:
$$\mathcal R_{nm}\frac{\phi_m}{\phi_n}
=\frac{1}{\sqrt{(2n{+}1)(2m{+}1)}}\sqrt{\frac{2n{+}1}{2m{+}1}}\,R
=\frac{R(s_n,s_m)}{2m+1}
\;\Longrightarrow\;
\boxed{\,c_n=\sum_m \frac{R(s_n,s_m)}{2m+1}\,}$$
independent of $n$'s own weight (verified against the explicit matrix product to
$10^{-16}$). This is **$\sim\!60\times$ tighter** than $\mathcal R\preceq\Theta^\*I$,
because $c_n$ is tiny exactly where $\tau_n$ is large ($c_0\approx10^{-2}$ at
$\varpi=10$, $5\times10^{-4}$ at $\varpi=100$) while $\max_nc_n\approx\Theta^\*$ is
attained near $n\sim\varpi$ where $\tau$ has already decayed.

Crucially $c_n$ depends on $n$ only through $a=s_n$, and the map $a\mapsto c$ is
**$\varpi$-independent** in the continuum limit (the discrete $c_n$ matches
$\tfrac12\int_0^\infty R(a,s)\,ds/s$ to $2\times10^{-4}$ at $\varpi=100,1000$). That
is what lets a *single* $a$-envelope serve every $\varpi$:
$$c(a)\;\le\;C(a):=\min\Big(C_1\frac{a^2}{1+a^2},\;C_\infty,\;\frac{C_2}{a}\Big),
\qquad C_1=1.5190,\ C_\infty=0.6375,\ C_2=1.6870 .$$
**Honest caveat:** $C_1,C_\infty,C_2$ are *calibrated* — measured as suprema over
$\varpi\in[2,400]$ and all $n$ (measured $1.51848$, $0.63716$, $1.68612$), then
padded. They are consistent with the continuum limit but are **not** closed-form
proven constants. Gate step [1] re-verifies $c_n\le C(s_n)$ at every run.
Note $\max_n c_n=0.6372>\Theta^\*$: $\Theta^\*$ bounds $\|\mathcal R\|$, **not** the
row sums, so using $\Theta^\*$ as a pointwise $c$-fill would be unsound.

### Prop. B (upper, **proven modulo the calibrated envelope $C$**, cost $O(\log\varpi)$ — the only route that reaches $\lambda\to0$)
Partition $n\ge1$ into blocks $[a,b)$. On each block, $G$ decreasing and the
**exact** digamma sum $\sum_{n=a}^{b-1}\frac1{2n+1}=\tfrac12[\psi(b{+}\tfrac12)-\psi(a{+}\tfrac12)]$ give
$$\sum_{n\in[a,b)}\tau_n^2\le G(s_a)^2\cdot\tfrac12\big[\psi(b{+}\tfrac12)-\psi(a{+}\tfrac12)\big],$$
with $\mathcal D_n\ge\mathcal D_{\rm blk}$ (use $n=b-1$; $\mathcal D$ increasing) and
$c_n\le c_{\rm blk}:=\sup_{[a,b)}C$ (piecewise-monotone, so the sup is at an endpoint
or the interior peak $a_{\rm peak}=\sqrt{C_\infty/(C_1-C_\infty)}$). Then
$F(h)\le\sum_{\rm blk}2S_{\rm blk}/(h+\mathcal D_{\rm blk}-c_{\rm blk})+\text{tail}$,
and the root of $=1$ upper-bounds $h$. Using the digamma identity instead of the
cruder $\tfrac12\ln\frac{2b-1}{2a-1}$ (which overshoots by $1.5\times$ on small
blocks) shrank the overshoot by $\sim\!2.4\times$.

---

## 3. Why $\lambda=0.05$ forces an analytic route (a structural limit)

$\lambda=0.05$ means $h=20$, i.e. $\varpi\sim e^{19.4}\approx3\times10^{8}$, and the
head matrix would need $N\sim10^{10}$. **No matrix method can reach the deep weak
tail.** This is not a budget problem; it is why Prop. B (cost $O(\log\varpi)$)
exists. At $\varpi=4.6\times10^{8}$ Prop. B returns in $\sim\!0.7$ s with
$\Delta h=0.056$.

Measured overshoot $h_{\rm hi}-h_{\rm true}$ (vs. a dense independent eigenvalue):

| route | $\varpi=3$ | $\varpi=10$ | $\varpi=50$ | cost | reaches $\lambda\!\to\!0$? |
|---|---|---|---|---|---|
| $\mathcal R\preceq\Theta^\*I$ (Prop. A) | $+0.524$ | $+0.561$ | — | $O(M)$ | yes |
| Prop. B, uniform $\Theta^\*$, log-bound | $+1.27$ | $+1.53$ | $+1.55$ | $O(\log\varpi)$ | yes |
| **Prop. B, envelope $C$ + digamma** | $+0.434$ | $+0.287$ | $+0.170$ | $O(\log\varpi)$ | **yes** |
| **Prop. A′ (full Schur, $O(\varpi M)$)** | $+0.043$ | $+0.021$ | $+0.011$ | $O(\varpi M)$ | no |
| **Prop. L/U head-$K{=}8$** | $\pm4\!\times\!10^{-4}$ | $\pm7\!\times\!10^{-4}$ | — | $O(N^2)$ | no |

---

## 4. Band map and rigor (unambiguous)

| band | $\varpi$ | route | $h$ rel. width | rigor |
|---|---|---|---|---|
| strong | $\le0.14$ | Weyl $\mathfrak g(2)\varpi^2$ envelope | few $\%$ | **proven analytic** (Cor ii) |
| crossover | $(0.14,3]$ | `node_cert`: Rayleigh + Kato–Temple | $\sim\!5\times10^{-4}$ | **certificate** (finite-$N$) |
| head-$K$ | $(3,25]$ | Prop. L + Prop. U, $K{=}8$ | $\sim\!3\times10^{-4}$ | **certificate** (finite-$N$) |
| deep | $>25$ | Cor (iii) + Prop. B | $\Delta h\approx0.17\!\to\!0.04$ | **proven modulo calibrated $C$** |

Both sides of every band are proven or certificate-level. **No band anywhere relies
on an asymptotic-only upper side** — which is precisely the gap this note closes.

The crossover band uses `node_cert` rather than head-$K$ because Props. L/U both
divide by $h-\Theta^\*$ and $h(\varpi)<\Theta^\*$ for $\varpi\lesssim0.93$; the
Rayleigh/Kato–Temple route has no such restriction.

---

## 5. Results

Gate: **`PASS 30  FAIL 0`** (quick mode, 42 s, single-thread). Checks cover:
the envelope soundness $c_n\le C(s_n)$; enclosure of $h_{\rm true}$ at seven
$\varpi$ against the **independent** `H_mat` builder (KAY entry formulas, *not* the
(M1) decomposition — a genuinely different discretization, per the quadrature-bias
lesson that multi-dps agreement proves nothing); the inversion direction; seam
overlap at six $\varpi$; $F_{\rm lo}\le F_{\rm hi}$ at every tabulated $\lambda$;
and both asymptotic prefactors.

Enclosure vs. the independent dense eigenvalue:

| $\varpi$ | $h_{\rm lo}$ | $h_{\rm true}$ | $h_{\rm hi}$ | route |
|---|---|---|---|---|
| 0.3 | 0.109758 | 0.109779 | 0.109799 | crossover-KT |
| 1.0 | 0.715777 | 0.715857 | 0.715937 | crossover-KT |
| 3.0 | 1.717748 | 1.717748 | 1.718590 | crossover-KT |
| 5.0 | 2.204773 | 2.204773 | 2.205529 | head-K8 |
| 10.0 | 2.873538 | 2.873538 | 2.874256 | head-K8 |
| 20.0 | 3.549054 | 3.549055 | 3.549792 | head-K8 |

Seam overlaps (both brackets contain $h$): all six tested pairs overlap.

**Certified $T_c/\omega_E$ bracket** (quick grid, 108 nodes over
$\varpi\in[0.1,10^{13}]$):

| $\lambda$ | $F_{\rm lo}$ | $F_{\rm hi}$ | rel. width | route | rigor |
|---|---|---|---|---|---|
| 0.050 | $3.746\times10^{-10}$ | $8.126\times10^{-10}$ | $1.17$ | deep-dyadic | proven-cal |
| 0.100 | $8.819\times10^{-6}$ | $1.299\times10^{-5}$ | $0.47$ | deep-dyadic | proven-cal |
| 0.300 | $8.057\times10^{-3}$ | $1.020\times10^{-2}$ | $0.27$ | head-K8 | certificate |
| 0.582 | $4.192\times10^{-2}$ | $5.965\times10^{-2}$ | $0.42$ | crossover-KT | certificate |
| 1.000 | $1.072\times10^{-1}$ | $1.206\times10^{-1}$ | $0.12$ | crossover-KT | certificate |
| 10.00 | $5.539\times10^{-1}$ | $6.228\times10^{-1}$ | $0.12$ | crossover-KT | certificate |
| 50.00 | $1.119$ | $1.415$ | $0.26$ | strong-Weyl | proven |

**Honest reading of these widths.** They are **much wider than ③'s $3\times10^{-3}$
headline**, and that is expected: ③ spends 2000 nodes on $\varpi\in[0.15,3]$ alone,
whereas this grid spreads 108 nodes over **thirteen decades** to reach
$\lambda=0.05$. The staircase inversion loses about one grid cell, so the $T_c$ width
here is **grid-limited, not bound-limited** — the underlying $h$-brackets are
$3\!-\!7\times10^{-4}$ wide in the certified bands. Refining the grid shrinks these
$\propto1/N$ **without** touching the mathematics. The two numbers answer different
questions: ③ answers "how tight can the crossover get?", this note answers "is
anything left uncertified anywhere?".

**Direct evidence that the width is grid-limited.** Re-running with the grid
densified only in the two *cheap* bands (crossover $\to$ 90 nodes, dyadic $\to$ 260
nodes; head-$K$ deliberately omitted, 350 nodes, 30 s) shrinks the widths by
$\approx5\times$ with **no change to any bound**:

| $\lambda$ | rel. width, 108-node grid | rel. width, densified cheap bands |
|---|---|---|
| 0.05 | $1.17$ | $\mathbf{0.23}$ |
| 0.10 | $0.47$ | $\mathbf{0.23}$ |
| 1.00 | $0.124$ | $\mathbf{0.039}$ |
| 20.0 | $0.124$ | $\mathbf{0.039}$ |

(The head-$K$ band, excluded from that particular refinement run, is what leaves a
gap at $\lambda\approx0.3\!-\!0.5$ there; it is covered in the gate proper. The
head-$K$ band is the only expensive one, $O(N^2)$ per node with $N$ up to 2400.)

**Asymptotic limits.** Weak: $(1/2\pi)e^{c_0}=0.2529996837$, reproducing
$T_c\to0.2530\,\omega_E e^{-1/\lambda}(1+b_1\lambda)$. Strong:
$T_c/\omega_E\to\varkappa\sqrt\lambda$ with
$\varkappa=\sqrt{\mathfrak g(2)}/(2\pi)\in[0.18272624774670,0.18272624778904]$
(CERT-D, width $4.2\times10^{-11}$), nested strictly inside the independent
CERT-B interval $[0.18272410720972,0.18272772871072]$.

Note this bracket **straddles the rounding boundary** of the published 10-dp
literal ($\ldots77\,|\,\ldots78$), so it encloses *neither* $0.1827262477$ nor
$0.1827262478$; the gate therefore asserts agreement to $10^{-9}$ and CERT-B
nesting rather than enclosure of a rounded literal, which would be a false check.

**Independent confirmation of $\varkappa$ (v2).** The inverted-bracket rows in
gate §[7] come from `strong_invert`, i.e. from the same Weyl envelope whose
constant is under test — self-consistent, not independent. Since strong coupling
means *small* $\varpi_\*$, which is cheap for the KAY-entry builder, the gate now
also solves $\lambda h(\varpi_\*)=1$ directly from `H_mat` with no envelope
anywhere and Richardson-extrapolates the $O(\varpi^2)$ correction:
$\varkappa(1600)=0.182687032$, $\varkappa(6400)=0.182716444$, Richardson limit
$0.182726248$ — agreeing with CERT-D to $6.6\times10^{-10}$ and sitting
$4.4\times10^{-2}$ away from the erroneous v1 value. This is the check that would
have caught the v1 slip on its own, and it is now a gate assertion with a teeth
test against the wrong constant.

> **Erratum (v2, 2026-09-15).** v1 stated $\varkappa=1/(2\pi\sqrt{\mathfrak g(2)})
> =0.13862$ and attributed $0.1827$ to a different frequency scale. That was a
> reciprocal-sqrt slip, caught in review. Derivation:
> $h(\varpi)\to\mathfrak g(2)\varpi^2$ gives
> $\varpi_\*=1/\sqrt{\mathfrak g(2)\lambda}$, hence
> $T_c/\omega_E=1/(2\pi\varpi_\*)=\sqrt{\mathfrak g(2)\lambda}/(2\pi)$ — the
> prefactor is $\sqrt{\mathfrak g(2)}/(2\pi)=0.1827262478$, **not** its
> reciprocal; the two differ by exactly the factor $\mathfrak g(2)=1.3181$. For an
> Einstein spectrum $\sqrt{\langle\omega^2\rangle}=\omega_E$ is the *same* scale,
> so no normalisation can absorb that factor. $0.1827$ is therefore the correct
> $\omega_E$-scale constant and **agrees** with Claim C
> ($T_c<0.1827\sqrt{\eta/M}$) and KAY Thm 1(c) $f(0)=0.1827262477\ldots$.
> The brackets themselves were never affected — they use the correct inversion
> ($\varpi_\*=\sqrt{u/\mathfrak g(2)}$); the error was confined to the reported
> limit constant in gate section [7] and to this discussion.

---

## 6. Proven / certificate map

| ingredient | rigor |
|---|---|
| $\lambda h(\varpi_\*)=1$, $T_c/\omega_E=1/2\pi\varpi_\*$ | proven (KAY, $\mu^\*=0$) |
| $h$ strictly increasing $\Rightarrow$ monotone inversion | **proven** (Cor i / Lemma M) |
| (M1) identity $\mathcal H=2\tau\tau^\top+\mathcal R-\mathcal D$ | proven; re-verified to $10^{-14}$ |
| $\|\mathcal R\|\le\Theta^\*$ | certified (Lemma M) |
| Prop. L lower (head Schur complement $\succeq$) | **proven inequality**, finite-$N$ inputs |
| Prop. U upper (Schur + geometric remainder) | **proven inequality**, finite-$N$ inputs |
| Prop. A′ Schur test $\mathcal R\preceq\mathrm{diag}(c)$ | **proven inequality**; weight collapse verified $10^{-16}$ — but the *value* returned uses the calibrated $C$ below, so any bracket built on it is `proven-cal`, not `proven` |
| envelope constants $C_1,C_\infty,C_2$ | **calibrated, not closed-form** — re-checked every run |
| Prop. B block bounds (digamma, $G$/$\mathcal D$ monotone) | **proven inequality**; uses calibrated $C$ |
| Cor (iii) lower $h\ge\ln\varpi+c_0-\frac1{12\varpi}-\frac3{16\varpi^2}$ | **proven analytic** |
| Weyl envelope, $\mathfrak g(2)$ CERT-D | **proven analytic** / certified constant |
| `node_cert` (Rayleigh + Kato–Temple) | certificate (inherited from ③) |

---

## 7. What did **not** work as hoped

1. **The inherited Theorem-W upper bound could not be reused** — defect (D1). Prop. U
   had to be derived to replace it. The published $K{=}8$ numbers are numerically
   fine; the logic was not.
2. **Widths are grid-limited, not bound-limited.** Spanning thirteen decades with a
   single monotone grid costs $\sim\!10^{-1}$ relative width per $\lambda$, versus
   $3\times10^{-3}$ for ③'s dense crossover-only grid. Nothing is wrong with the
   mathematics; the grid is simply coarse. A dense run refines this $\propto1/N$.
3. **The envelope constants are calibrated, not proven.** Closing this would mean
   proving $\sup_a c(a)\le C_\infty$ etc. in closed form from the rational kernel
   $R$ — tractable, but not done here. Flagged rather than glossed.
4. **A "min of two dominations" bug** appeared in this file's own first draft and was
   caught by an explicit counterexample test, now recorded in §0. It had also
   silently made an earlier $\Theta^\*$-fill unsound.
5. **v1 reported the strong prefactor reciprocal-flipped** ($0.13862$ for
   $\sqrt{\mathfrak g(2)}/2\pi=0.18273$) — see §5.

---

## 8. Reproduce

```
cd theory/derivations
OMP_NUM_THREADS=1 python p4_global_tc_allcoupling.py quick   # ~42 s, PASS 30 FAIL 0
OMP_NUM_THREADS=1 python p4_global_tc_allcoupling.py         # dense grid, tighter widths
```
All soundness checks are counted by `ck(...)`; the run prints `PASS N  FAIL M`.
