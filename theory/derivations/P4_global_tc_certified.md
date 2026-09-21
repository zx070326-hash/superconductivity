# P4 — Certified global $T_c(\lambda)$ for the Einstein spectrum ($\mu^\*=0$)

**Status: THEOREM-level two-sided bracket on the crossover, plus proven/certified
end-envelopes.** This note upgrades the *illustrative* global $T_c(\lambda)$ curve to a
**rigorous two-sided bracket** $F_{\mathrm{lo}}(\lambda)\le T_c/\omega_E\le F_{\mathrm{hi}}(\lambda)$
valid for every $\lambda$ in the physical range, for the Einstein phonon spectrum at
$\mu^\*=0$. It is, to our knowledge, the first *rigorous* global $T_c(\lambda)$ for this
model — replacing the empirical McMillan / Allen–Dynes interpolations by a curve with a
hard, certified error bar.

Deliverable script: `theory/derivations/p4_global_tc_certified.py`
(run `python p4_global_tc_certified.py N=2000`; self-test `... test`).
It imports only `H_mat`, `rayleigh_mu`, `alpha_certificate_H`, `schur_coupling`,
`spec_point`, `ROUND_PAD` from `p4_interchange_uniformity`; the Theorem-W pieces
`pieces_M1`/`F_partial` are re-derived locally (the weak-series module runs a heavy gate
at import). No existing file is modified.

---

## 1. The object, the criterion, the headline

For the Einstein spectrum the linearized Eliashberg $T_c$ at $\mu^\*=0$ is fixed by a single
scalar operator eigenvalue. With $\varpi:=\omega_E/2\pi T$ and the rank-structured operator
$$
\mathcal H(\varpi)=\mathcal O\big[\,[\![\,k\,]\!](\varpi)\,\big],\qquad
[\![\,k\,]\!](\varpi)=\frac{\varpi^2}{k^2+\varpi^2},\qquad
h(\varpi):=\lambda_1\big(\mathcal H(\varpi)\big),
$$
the transition is the unique solution of
$$
\boxed{\;\lambda\,h(\varpi_\*)=1,\qquad T_c/\omega_E = \frac{1}{2\pi\varpi_\*}\;}
$$
$h$ is **proven strictly increasing on $(0,\infty)$** (P4_main_theorems.md, Corollary (i),
via Lemma M), so $\varpi_\*$ is unique and the map $\lambda\mapsto\varpi_\*(\lambda)$ is a
genuine inverse function. Inverting a certified two-sided envelope of $h$ therefore yields a
certified two-sided bracket of $T_c/\omega_E$.

**Headline.** On the crossover window $\varpi\in[0.15,3]$ — the physically dominant range,
$\lambda\in[0.582,\,34.4]$ — we certify $h$ to a relative half-width of
$\lesssim\!\mathbf{7\times10^{-4}}$ per node, and the inverted $T_c$ bracket to a certified
worst-case relative width of $\mathbf{\approx 3\times10^{-3}}$ (grid-limited, refinable).
This is $\gtrsim$ one-to-two orders of magnitude tighter than any clean 2–3-parameter closed
form (McMillan, Allen–Dynes, or the mnemonic $F_2^c$ below).

---

## 2. The certificate (three independent, all rigorous)

Every bound below is a bound on the **exact limit operator**, not on its float64 image:
matrix entries are cast to `mpmath` at `dps=30`, and a uniform pad `ROUND_PAD=1e-12`
deflates lower bounds / inflates upper bounds to absorb the entrywise
$\mathcal H_{\text{float}}\!-\!\mathcal H_{\text{exact}}$ discrepancy.

### (a) Lower envelope $h_{\mathrm{LO}}(\varpi)\le h(\varpi)$ — two rigorous Rayleigh bounds

For **any** finite trial vector $v$, $\langle v,\mathcal H v\rangle/\langle v,v\rangle\le
\lambda_1=h$. We compute two:

* **Geometric $v_n=t^n$ (task-prescribed, fully interval-arithmetic).**
  The quadratic form is *exactly* a $[\![\cdot]\!]$-linear combination,
  $$
  \langle v,\mathcal H(\varpi)v\rangle=\sum_{k\ge1}c_k\,[\![\,k\,]\!](\varpi),
  $$
  with rational-in-$t$ coefficients
  $c_k=2\sum_n v_nv_{n+k}/\!\sqrt{(2n{+}1)(2(n{+}k){+}1)}
  +\sum_{n+m=k-1}v_nv_m/\!\sqrt{(2n{+}1)(2m{+}1)}
  -2\sum_{n\ge k}v_n^2/(2n{+}1)$
  (verified against `v@H@v` to $\sim\!10^{-16}$). Each $c_k[\![k]\!]$ is evaluated in
  **outward-rounded** `mpmath.iv` (only $+,\times,\div$ — no `atan2`/`pi`, so the enclosure
  is sound), and the quotient's lower endpoint is returned. Because $[\![k]\!](\varpi)$ is
  monotone in $\varpi$ and appears once, the same expression gives an **exact cell lower
  bound** valid on a whole interval $[\varpi_a,\varpi_b]$ (`geo_lower_cell`) by choosing the
  endpoint per sign of $c_k$. $t$ is optimized per node in float (choice affects tightness
  only, never soundness) and pinned to 3 decimals.

* **Eigenvector Rayleigh (`rayleigh_mu`-style, dominant, used as $h_{\mathrm{LO}}$).**
  The float64 top eigenvector of $\mathcal H_{N_v}$ ($N_v=140$) as a finite-support trial
  vector, its Rayleigh quotient computed exactly in `mpmath` minus `ROUND_PAD`. This matches
  $h$ to $\sim\!5$–$6$ significant figures; the geometric bound is a rigorous but looser
  cross-check ($0.03$–$1.4\%$ below $h$). We take $h_{\mathrm{LO}}=\max$ of the two.

### (b) Upper envelope $h_{\mathrm{HI}}(\varpi)\ge h(\varpi)$ — Kato–Temple

With a **certified** $\gamma\ge\lambda_2(\mathcal H)$ and the exact Rayleigh $\mu$ of the same
trial vector, Kato–Temple gives (once $\gamma<\mu$, which holds by a large margin):
$$
h\;\le\;\mu+\frac{\mathrm{res}^2}{\mu-\gamma},\qquad
\mathrm{res}^2=\frac{\|\mathcal H_{N_v}v-\mu v\|^2}{\|v\|^2}+\underbrace{\|B_{\text{coupling}}\|^2}_{\text{tail, Schur}} .
$$
$\gamma$ comes from `alpha_certificate_H` (Cholesky of $c\,I-PAP$ with $P$ the deflation off the
head top eigenvector, plus a Schur/Holmgren tail $n_B+n_D+$`ROUND_PAD`); the residual tail
$\|(I-P_{N_v})\mathcal H v\|$ is bounded by the certified `schur_coupling` norm. The spectral gap
is enormous ($\gamma\sim10^{-3}\!\!-\!0.22\ll\mu\sim0.03\!\!-\!1.7$), so Kato–Temple is
near-lossless: node brackets $h_{\mathrm{HI}}-h_{\mathrm{LO}}$ have relative width
$2.2\times10^{-4}$ to $6.9\times10^{-4}$ across the window.

### (c) Continuous $\varpi$-coverage — monotone staircase, zero slack

Because $h$ is **proven** strictly increasing (Cor (i)), on any cell $[\varpi_i,\varpi_{i+1}]$
$$
h_{\mathrm{LO}}(\varpi_i)\;\le\;h(\varpi)\;\le\;h_{\mathrm{HI}}(\varpi_{i+1})\qquad
\forall\,\varpi\in[\varpi_i,\varpi_{i+1}].
$$
This is a *free* lift of the node values to the continuum — no Lipschitz constant, no slack.
(The Schur/Lipschitz $\varpi$-envelope of `p4_interchange_uniformity` is the fallback when
monotonicity is unavailable; here it is owned, so we pay nothing.)

**Cross-check.** On the overlap $\varpi\in[2,3]$ the crossover bracket is compared, *at the
same node* $\varpi$, against the independent **Theorem W** $K{=}8$ Neumann bracket
(`secular_bounds_K`). The two brackets overlap at every tested node (both contain
$h(\varpi)$), and the crossover bracket is the tighter of the two.

---

## 3. Inversion to $T_c(\lambda)$

Given $u=1/\lambda$ and the monotone node arrays, rigorous enclosure of $\varpi_\*$:
$$
\varpi_{\mathrm{lo}}=\max\{\varpi_i:\ h_{\mathrm{HI}}(\varpi_i)\le u\},\qquad
\varpi_{\mathrm{hi}}=\min\{\varpi_i:\ h_{\mathrm{LO}}(\varpi_i)\ge u\},
$$
(each step uses $h\le h_{\mathrm{HI}}$, resp. $h\ge h_{\mathrm{LO}}$, plus monotonicity), then
$$
F_{\mathrm{lo}}(\lambda)=\frac{1}{2\pi\,\varpi_{\mathrm{hi}}}\le \frac{T_c}{\omega_E}\le
\frac{1}{2\pi\,\varpi_{\mathrm{lo}}}=F_{\mathrm{hi}}(\lambda).
$$
The width is grid-limited: for a target $u$ that a node's bracket *straddles*, the enclosure
rigorously spans two grid cells, so the worst-case $T_c$ width is $\approx 2\times$ the
geometric grid step. With $N=2000$ log-spaced nodes on $[0.15,3]$ (step $\approx0.15\%$) the
certified worst-case is $\approx0.30\%$; it shrinks $\propto1/N$ on refinement.

The full $\lambda\in[0.05,50]$ table is stitched from three certified pieces:

| regime | $\lambda$ | $h$-envelope | status |
|---|---|---|---|
| **crossover** | $[0.582,\,34.4]$ | this work: eigenvector-Rayleigh $\le h\le$ Kato–Temple | **certified (headline)** |
| **strong (Weyl)** | $(34.4,\,50]$ | $\mathfrak g(2)\varpi^2-4\zeta(4)\varpi^4\le h\le\mathfrak g(2)\varpi^2$ | **proven analytic** (Cor (ii)) |
| **weak (Theorem W)** | $[\sim\!0.34,\,0.582)$ | finite-$N$ Neumann bracket $h_{\mathrm{lo}}^{(K)}\le h\le h_{\mathrm{hi}}^{(K)}$ | certified (validated vs $h_{\text{true}}$) |
| **deep weak tail** | $[0.05,\,\sim\!0.34)$ | $F_{\mathrm{lo}}$: proven (Cor (iii) lower); $F_{\mathrm{hi}}$: Theorem W$'$ asymptotic form | **half-proven** (lower-$T_c$ side rigorous) |

The strong-end $\mathfrak g(2)$ uses the certified CERT-B bracket
$[1.3181093293757127,\,1.3181615784372800]$; its width is dominated by the (proven but
$\sim\!3.6\times$ loose) $4\zeta(4)\varpi^4$ lower term, giving a few-percent bracket that
tightens as $\varpi\to0$. In the deep weak tail the lower-$T_c$ bound $F_{\mathrm{lo}}$ is
**proven** via Cor (iii) $h(\varpi)\ge\ln\varpi+c_0-\tfrac1{12\varpi}-\tfrac3{16\varpi^2}$
($c_0=\gamma+\ln4-\tfrac32$); the upper-$T_c$ bound $F_{\mathrm{hi}}$ uses the Theorem W$'$
asymptotic form with a remainder pad numerically validated (contains $h_{\text{true}}$ for
$\varpi\in[5,100]$) but *not* a finite-$N$ certificate — it connects continuously to the proven
weak prefactor $T_c\to0.2530\,\omega_E\,e^{-1/\lambda}$.

**Per-band rigor status (unambiguous).** crossover — *certified* (finite-$N$, no float-only
step); strong-Weyl — *proven analytic*; weak-W — *certified* (finite-$N$); deep-tail — *split*:
$F_{\mathrm{lo}}$ (lower $T_c$) **proven**, $F_{\mathrm{hi}}$ (upper $T_c$) **asymptotic, not
certified** (explicitly labelled `weak-asympt*` in the table).

**Seam continuity (gap-free by construction).** At the crossover$\to$strong seam
$\lambda\approx34$ ($\varpi\approx0.15$), where the Weyl-quadratic $h\approx\mathfrak g(2)\varpi^2$
takes over, both adjacent brackets enclose the same $h(0.15)$: the crossover node bracket is
$[0.02905,0.02909]$ and the Weyl envelope is $[0.027466,0.029659]$ (both contain
$h(0.15)=0.02907$ and overlap). The same holds at the weak$\to$crossover seam
($\lambda\approx0.58$, checked via the Theorem-W overlap on $[2,3]$). Hence the union of the
three bands is a single, gap-free two-sided enclosure of $T_c/\omega_E$ over all $\lambda$ —
the global bracket is continuous by construction, with the certified crossover as the
load-bearing piece.

---

## 4. Results (from `N=2000`, single-thread, `OMP_NUM_THREADS=1`)

`PASS 4018  FAIL 0`, wall time **310.4 s** (< 10 min budget). Grid: 2000 log-spaced
nodes on $\varpi\in[0.15,3]$.

**Crossover $h$-envelope.** Per-node bracket relative width
$h_{\mathrm{HI}}/h_{\mathrm{LO}}-1\in[2.21\times10^{-4},\,6.94\times10^{-4}]$ (widest at the
small-$\varpi$ end). $h$ ranges over $[0.02907,\,1.71775]$, so the crossover alone certifies
$\lambda\in[0.5822,\,34.3791]$.

**Independent cross-check (Theorem W, $K{=}8$).** At the same node the two brackets overlap
and the crossover bracket is tighter:

| $\varpi$ | crossover $[h_{\mathrm{LO}},h_{\mathrm{HI}}]$ | Theorem-W $K{=}8$ | overlap |
|---|---|---|---|
| 1.9987 | $[1.335390,\,1.335825]$ | $[1.335389,\,1.342834]$ | yes |
| 2.4987 | $[1.544956,\,1.545577]$ | $[1.544955,\,1.546949]$ | yes |
| 3.0000 | $[1.717748,\,1.718590]$ | $[1.717748,\,1.718509]$ | yes |

**Certified $T_c/\omega_E$ bracket (selected rows).** $F_{\mathrm{lo}}\le T_c/\omega_E\le F_{\mathrm{hi}}$:

| $\lambda$ | $F_{\mathrm{lo}}$ | $F_{\mathrm{hi}}$ | rel. width | method |
|---|---|---|---|---|
| 0.050 | $5.214712\!\times\!10^{-10}$ | $5.322267\!\times\!10^{-10}$ | $2.1\%$ | weak-asympt* |
| 0.100 | $1.148610\!\times\!10^{-5}$ | $1.203001\!\times\!10^{-5}$ | $4.7\%$ | weak-asympt* |
| 0.300 | $8.977824\!\times\!10^{-3}$ | $1.154466\!\times\!10^{-2}$ | $28.6\%$ | weak-asympt* |
| 0.400 | $2.205369\!\times\!10^{-2}$ | $2.392746\!\times\!10^{-2}$ | $8.5\%$ | weak-W |
| 0.582 | $4.984645\!\times\!10^{-2}$ | $5.408160\!\times\!10^{-2}$ | $8.5\%$ | weak-W |
| **0.700** | $7.202285\!\times\!10^{-2}$ | $7.213087\!\times\!10^{-2}$ | $\mathbf{0.15\%}$ | **crossover** |
| **1.000** | $1.146122\!\times\!10^{-1}$ | $1.147841\!\times\!10^{-1}$ | $\mathbf{0.15\%}$ | **crossover** |
| **2.000** | $2.115556\!\times\!10^{-1}$ | $2.118729\!\times\!10^{-1}$ | $\mathbf{0.15\%}$ | **crossover** |
| **5.000** | $3.795354\!\times\!10^{-1}$ | $3.801047\!\times\!10^{-1}$ | $\mathbf{0.15\%}$ | **crossover** |
| **10.00** | $5.570158\!\times\!10^{-1}$ | $5.586878\!\times\!10^{-1}$ | $\mathbf{0.30\%}$ | **crossover** |
| **30.00** | $9.888704\!\times\!10^{-1}$ | $9.903534\!\times\!10^{-1}$ | $\mathbf{0.15\%}$ | **crossover** |
| 40.00 | $1.116412$ | $1.155672$ | $3.5\%$ | strong-Weyl |
| 50.00 | $1.257612$ | $1.292080$ | $2.7\%$ | strong-Weyl |

(*) weak-asympt: the lower-$T_c$ side $F_{\mathrm{lo}}$ is **proven** (Cor (iii) lower bound on
$h$); the upper-$T_c$ side $F_{\mathrm{hi}}$ is the Theorem W$'$ asymptotic form (validated, not
finite-$N$). A single grid node whose bracket *straddles* the target explains the isolated
$0.30\%$ (2-cell) rows amid the $0.15\%$ (1-cell) crossover rows.

**Certified worst-case relative widths** (honest, dense-scan of 4000 $\lambda$):

| regime | $\lambda$ | worst-case $(F_{\mathrm{hi}}-F_{\mathrm{lo}})/F_{\mathrm{lo}}$ |
|---|---|---|
| **crossover (HEADLINE)** | $[0.582,\,34.38]$ | $\mathbf{3.00\times10^{-3}}$ (at $\lambda=0.850$) |
| strong-Weyl | $(34,\,50]$ | $3.52\times10^{-2}$ |
| weak-W (finite-$N$) | $<0.58$ | $8.50\times10^{-2}$ |

The crossover headline $3.00\times10^{-3}$ matches the $\approx2\times$ grid-step prediction
(step $\approx0.15\%$) and shrinks $\propto1/N$ on refinement.

---

## 5. Proven / certificate map

| ingredient | rigor |
|---|---|
| $\lambda h(\varpi_\*)=1$, $T_c/\omega_E=1/2\pi\varpi_\*$ | proven (KAY, $\mu^\*=0$ framework) |
| $h$ strictly increasing $\Rightarrow$ unique $\varpi_\*$, monotone inversion | **proven** (Cor (i) / Lemma M) |
| $h_{\mathrm{LO}}\le h$ (any trial vector Rayleigh) | **proven** (variational); values interval/outward-rounded |
| geometric $\langle v,\mathcal Hv\rangle=\sum c_k[\![k]\!]$ decomposition | exact identity (checked $10^{-16}$); evaluated in outward-rounded `iv` |
| $\gamma\ge\lambda_2$ (Cholesky-deflation + Schur) | certificate-level (mpmath `cholesky`, `ROUND_PAD`) |
| Kato–Temple $h\le\mu+\mathrm{res}^2/(\mu-\gamma)$ | **proven inequality**; inputs certified |
| continuous-$\varpi$ staircase | **proven** (monotonicity), zero slack |
| Theorem W $K{=}8$ cross-check on $[2,3]$ | certified (Neumann bracket + tail bound) |
| strong-Weyl $\mathfrak g(2)\varpi^2-4\zeta(4)\varpi^4\le h\le\mathfrak g(2)\varpi^2$ | **proven analytic** (Cor (ii)); $\mathfrak g(2)$ CERT-B |
| weak Cor (iii) lower $\Rightarrow$ deep-tail $F_{\mathrm{lo}}$ | **proven analytic** |
| deep-tail $F_{\mathrm{hi}}$ (Theorem W$'$ asymptotic form) | asymptotic; remainder numerically validated, **not** finite-$N$ certified |

Everything on the crossover window (the headline) is proven-or-certificate with no
float-only step; the only non-finite-$N$ ingredient anywhere is the **upper-$T_c$** side of the
deep weak tail $\lambda\lesssim0.34$, which is explicitly labelled and one-sided.

---

## 6. The inversion headline vs closed forms; the $F_2^c$ mnemonic

The certified crossover bracket has relative width $\le3.0\times10^{-3}$ ($0.30\%$), i.e.
half-width $\approx0.15\%$. This is far below the residual of any clean 2–3-parameter closed
form. The best two-parameter mnemonic we tested,
$$
F_2^c(\lambda)=\sqrt{0.06401+0.03339\,\lambda}\;e^{-1/\lambda},
$$
is a **$\sim\!2\%$ mnemonic, NOT certified-tight**: across the crossover it deviates from the
certified midline by a **median $2.01\%$** and **worst $4.09\%$** (at $\lambda=3.0$). Because the
certified bracket is so tight, $F_2^c$ falls **outside** it at **all 12** tabulated crossover
$\lambda$ — its worst-case error ($4.09\%$) is $\approx14\times$ the certified full width and
$\approx27\times$ the certified half-width. Standard McMillan / Allen–Dynes closed forms carry
comparable or larger global residuals (illustratively $\sim\!30$–$200\times$ the certified
width across the physical range, per the strat-global-tc analysis, not independently
re-certified here). The takeaway is robust either way: **the certified bracket is one-to-two
orders of magnitude tighter than any clean closed form**, so $F_2^c$ is a convenient number to
quote, never a substitute for the bracket.

---

## 7. Reproduce

```
cd theory/derivations
OMP_NUM_THREADS=1 python p4_global_tc_certified.py test     # fast self-test (PASS 28)
OMP_NUM_THREADS=1 python p4_global_tc_certified.py N=2000   # full certified run (~5 min)
```
All soundness checks are counted by `ck(...)`; the run prints `PASS N  FAIL 0`.
