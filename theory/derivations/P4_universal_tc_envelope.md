# P4 — A universal, certified $T_c$ upper envelope over **all** phonon spectra ($\mu^\*=0$)

**Status: THEOREM (one-sided), load-bearing on the certified Einstein bracket (③).**
For *every* coupling $\lambda>0$ and *every* phonon spectrum $P$ (probability measure on
$\mathbb R_+$ with finite second moment $\langle\omega^2\rangle$, $P\neq\delta_0$), we
certify a single scalar upper bound
$$
\boxed{\;T_c(P;\lambda)\;\le\;F_{\mathrm{hi}}^{E}(\lambda)\,\sqrt{\langle\omega^2\rangle}\;}
$$
where $F_{\mathrm{hi}}^{E}(\lambda)$ is the **certified upper curve of $T_c/\omega_E$ vs
$\lambda$ for the single-frequency (Einstein) spectrum** — deliverable ③
(`p4_global_tc_certified.{md,py}`). The bound is **strictly sharper at every finite
$\lambda$** than the asymptotic ceiling $T_c\le C_\infty\sqrt{\lambda\langle\omega^2\rangle}$,
$C_\infty=(1/2\pi)\sqrt{\mathfrak g(2)}=0.1827262477\ldots$, which is tight only as
$\lambda\to\infty$.

Deliverable script: `theory/derivations/p4_universal_tc_envelope.py`
(run `OMP_NUM_THREADS=1 python p4_universal_tc_envelope.py`; prints `PASS 66  FAIL 0`,
$\approx 73$ s). It imports `H_mat` from `p4_interchange_uniformity` (committed base) and the
certified inversion pieces (`node_cert`, `invert_envelope`, `tc_from_w`, `strong_invert`,
`_weak_invert`, `_weak_asympt_invert`, `G2_LO`, `G2_HI`, `TWO_PI`) from
`p4_global_tc_certified` (deliverable ③, commit `1184684`, **inherited certificate-level**).
No existing file is modified.

---

## 1. The exact statement — with the correct moment

The extremal input is **Theorem 2 of `P4_main_theorems.md`** (Einstein extremality, all
$\lambda$), quoted verbatim (lines 65–68):

> **定理 2（Einstein 极值，全 λ）.** 对所有 $P\in\mathcal P$、$T>0$:
> $\;\mathfrak k(P,T)\le h(\varpi_{\mathrm{rms}})$, $\varpi_{\mathrm{rms}}:=\sqrt{\langle\omega^2\rangle}/2\pi T$, **等号 $\Longleftrightarrow$ $P$ 是点质量**。
> 等价：固定 $\lambda$ 与 $\langle\omega^2\rangle$，Einstein 谱 $\delta_{\sqrt{\langle\omega^2\rangle}}$ 是线性化 Eliashberg $T_c$ 的极大者；$T_c(\lambda,P)\le T_c^E(\lambda,\sqrt{\langle\omega^2\rangle})$ 对**所有** $\lambda$。
> *证*（**核级 Jensen**）：记 $x:=\varpi^2$。$F_T(k)=\int x/(k^2+x)\,P(d\omega)$，而
> $x\mapsto x/(k^2+x)$ 于 $[0,\infty)$ 严格凹 $\Longrightarrow$ 逐 $k$ 有
> $F_T(k)\le\langle x\rangle/(k^2+\langle x\rangle)=[\![k]\!](\varpi_{\mathrm{rms}})=:G(k)$ …
> 引理 M $\Longrightarrow \mathfrak k(P,T)=\lambda_1(\mathcal O[F_T])\le\lambda_1(\mathcal O[G])=h(\varpi_{\mathrm{rms}})$；严格性 $\Longleftrightarrow G(1)>F(1)\Longleftrightarrow x$ 非 $P$-a.s. 常值 $\Longleftrightarrow P$ 非点质量。∎

**The fixed moment is *exactly* the second moment $\langle\omega^2\rangle$.** The Jensen
step averages $x=\varpi^2=(\omega/2\pi T)^2$, so $\langle x\rangle=\langle\omega^2\rangle/(2\pi T)^2$
and the maximizing Einstein spectrum sits at the **root-mean-square frequency**
$$
\omega_E \;=\; \sqrt{\langle x\rangle}\,(2\pi T)\;=\;\sqrt{\langle\omega^2\rangle}\;=\;\omega_{\mathrm{rms}} .
$$
No $\langle\omega\rangle$, no logarithmic-average $\mathfrak g(2)$-form: the matched moment is
$\langle\omega^2\rangle$ and nothing else. (The definitions used — $\mathcal O[f]$ elementwise
linear, $[\![k]\!](\varpi)=\varpi^2/(k^2+\varpi^2)$, $h(\varpi)=\lambda_1(\mathcal O[[\![\cdot]\!](\varpi)])$,
$\mathfrak k(P,T)=\lambda_1(\mathcal O[F_T])$ — are `P4_main_theorems.md` §0, lines 16–18.)

### 1.1 Proof of the universal envelope

Combine three facts, each already established:

1. **Extremality (analytic).** Theorem 2 ⟹ $T_c(P;\lambda)\le T_c^{E}(\lambda;\sqrt{\langle\omega^2\rangle})$
   for every $P\in\mathcal P$ and every $\lambda>0$. (The crossing $T_c$ is where
   $\lambda\,\mathfrak k=1$; since $\mathfrak k(P,\cdot)\le\mathfrak k_E(\cdot)$ pointwise and
   both are strictly decreasing in $T$ — Theorem 1(b) — the crossing for $P$ occurs at
   $T\le$ the crossing for the matched Einstein spectrum.)
2. **Scale form of the Einstein $T_c$.** For the Einstein spectrum, $T_c^E$ is homogeneous of
   degree $1$ in $\omega_E$: $T_c^E(\lambda;\omega_E)=F^E(\lambda)\,\omega_E$, where
   $F^E(\lambda)=T_c/\omega_E$ solves $\lambda\,h(\varpi_\*)=1$, $F^E=1/(2\pi\varpi_\*)$.
3. **Certified upper bracket (③).** $F^E(\lambda)\le F_{\mathrm{hi}}^E(\lambda)$, the certified
   upper curve of deliverable ③, obtained by inverting a certified two-sided envelope of $h$
   through the **proven** strict monotonicity of $h$ (Corollary (i), via Lemma M).

Chaining, for every $P$ and $\lambda$,
$$
T_c(P;\lambda)\;\overset{(1)}{\le}\;T_c^E(\lambda;\sqrt{\langle\omega^2\rangle})
\;\overset{(2)}{=}\;F^E(\lambda)\sqrt{\langle\omega^2\rangle}
\;\overset{(3)}{\le}\;F_{\mathrm{hi}}^E(\lambda)\sqrt{\langle\omega^2\rangle}. \qquad\blacksquare
$$

Steps (1)–(2) are **analytic** (Theorem 2 + Theorem 1(b), both proven with no
computer-assisted component). Step (3) is the **inherited certificate ③**. Nothing in the
reduction is fudged: the moment matching is exact, and the whole $P$-dependence collapses to
the single factor $\sqrt{\langle\omega^2\rangle}$.

The gate validates the load-bearing inequality (1) **directly** on genuine multi-frequency
$P$ (check group [A]): for discrete $P=\{(p_i,\varpi_i)\}$, elementwise linearity gives
$\mathcal O[F_T]=\sum_i p_i\,\mathcal O[[\![\cdot]\!](\varpi_i)]=\sum_i p_i\,H_{\mathrm{mat}}(\varpi_i)$,
and we confirm $\mathfrak k(P,T)=\lambda_1(\mathcal O[F_T])\le h(\varpi_{\mathrm{rms}})$ with a
large margin for spread spectra and *exact equality* on a point mass.

---

## 2. Sharpening vs the $C_\infty$ ceiling — analytic, then certified

**Theorem 1(c)** (`P4_main_theorems.md` line 55) gives, for all $P,\lambda$,
$T_c(P;\lambda)<\tilde T_c=(1/2\pi)\sqrt{\mathfrak g(2)\langle\omega^2\rangle\lambda}
=C_\infty\sqrt{\lambda\langle\omega^2\rangle}$, with
$C_\infty=(1/2\pi)\sqrt{\mathfrak g(2)}$. This ceiling is a valid bound for *all* $\lambda$ but
is **tight only as $\lambda\to\infty$**; it is the $G(k)=\langle\varpi^2\rangle/k^2$ limit of
Theorem 2 (the Weyl envelope $h(\varpi)\le\mathfrak g(2)\varpi^2$).

The universal envelope is strictly below this ceiling **at every finite $\lambda$**, and the
reason is analytic: for the Einstein spectrum $T_c/\omega_E=1/(2\pi\varpi_\*)$ with
$\lambda h(\varpi_\*)=1$, so
$$
\frac{F^E(\lambda)}{C_\infty\sqrt{\lambda}}
=\frac{1/(2\pi\varpi_\*)}{(1/2\pi)\sqrt{\mathfrak g(2)\lambda}}
=\frac{1}{\varpi_\*\sqrt{\mathfrak g(2)\lambda}}
=\sqrt{\frac{r(\varpi_\*)}{\mathfrak g(2)}},\qquad r(\varpi):=\frac{h(\varpi)}{\varpi^2}.
$$
By **Corollary (ii)** (`P4_main_theorems.md` line 73): $r(\varpi)<\mathfrak g(2)$ strictly, and
$r$ is strictly decreasing on $(0,\infty)$ with $r(0^+)=\mathfrak g(2)$. Hence the ratio is
$<1$ for every finite $\lambda$ and $\to1$ only as $\lambda\to\infty$ ($\varpi_\*\to0$). **The
$C_\infty$ ceiling is loose at finite $\lambda$; the universal envelope closes that gap.**
Note $\pi$ **cancels** in the ratio — no transcendental enters the comparison.

**Certified test (sound).** The gate certifies $F_{\mathrm{hi}}^E(\lambda)\le C_\infty\sqrt\lambda$
without any $\sqrt{}$ or $\pi$: with $\varpi_{\mathrm{lo}}\le\varpi_\*$ (certified by ③) and
$G_{2,\mathrm{LO}}=1.3181093293757127$ (CERT-B lower on $\mathfrak g(2)$),
$$
\varpi_{\mathrm{lo}}^2\,G_{2,\mathrm{LO}}\,\lambda\ \ge\ 1
\quad\Longrightarrow\quad
F_{\mathrm{hi}}^E(\lambda)=\frac{1}{2\pi\varpi_{\mathrm{lo}}}
\le\frac{1}{2\pi}\sqrt{G_{2,\mathrm{LO}}\,\lambda}
\le C_\infty\sqrt\lambda,
$$
using only $\times$ and a comparison of exact `mpf` quantities (the last step because the true
$\mathfrak g(2)\ge G_{2,\mathrm{LO}}$). This is the most conservative form of the claim.

### 2.1 Certified table (gate output)

`ratio` $=F_{\mathrm{hi}}^E/(C_\infty\sqrt\lambda)$ (a certified **upper** bound, using
$G_{2,\mathrm{LO}}$); `sharpen` $=(C_\infty\sqrt\lambda)/F_{\mathrm{hi}}^E$ (a certified
**lower** bound on how much looser the ceiling is). All $\lambda\le30$ satisfy the certified
sharpening test $\varpi_{\mathrm{lo}}^2G_{2,\mathrm{LO}}\lambda\ge1$.

| $\lambda$ | $F_{\mathrm{lo}}^E$ | $F_{\mathrm{hi}}^E$ (envelope) | $C_\infty\sqrt\lambda$ (ceiling) | ratio $\le$ | sharpen $\ge$ | band |
|---:|---:|---:|---:|---:|---:|:--|
| 0.30 | $8.977824\times10^{-3}$ | $1.154466\times10^{-2}$ | $1.000833\times10^{-1}$ | $0.11535$ | $8.669$ | weak-asympt\* |
| 0.50 | $3.902908\times10^{-2}$ | $4.234515\times10^{-2}$ | $1.292070\times10^{-1}$ | $0.32773$ | $3.051$ | weak-W |
| 1.00 | $1.144820\times10^{-1}$ | $1.148264\times10^{-1}$ | $1.827262\times10^{-1}$ | $0.62841$ | $1.591$ | crossover |
| 2.00 | $2.114872\times10^{-1}$ | $2.121234\times10^{-1}$ | $2.584139\times10^{-1}$ | $0.82088$ | $1.218$ | crossover |
| 5.00 | $3.794728\times10^{-1}$ | $3.806143\times10^{-1}$ | $4.085883\times10^{-1}$ | $0.93155$ | $1.074$ | crossover |
| 10.0 | $5.571087\times10^{-1}$ | $5.587846\times10^{-1}$ | $5.778311\times10^{-1}$ | $0.96705$ | $1.034$ | crossover |
| 30.0 | $9.881502\times10^{-1}$ | $9.911227\times10^{-1}$ | $1.000833\times10^{0\phantom{-}}$ | $0.99031$ | $1.010$ | crossover |

**Reading.** At $\lambda=1$ the $C_\infty$ ceiling overstates the true universal bound by
$\approx59\%$; the overstatement falls to $\approx7.4\%$ ($\lambda=5$), $\approx3.4\%$
($\lambda=10$), $\approx1.0\%$ ($\lambda=30$), and $\to0$ only as $\lambda\to\infty$. At weak
coupling the gap is enormous ($\times 8.7$ at $\lambda=0.3$) because $T_c/\omega_E$ is
exponentially small, $\sim e^{-1/\lambda}$, far below the power-law ceiling.

**Note on the reported numbers.** The gate recomputes $\varpi_\*$ with a *local* certified grid
(fast); deliverable ③'s $N{=}2000$ run gives marginally tighter $F_{\mathrm{hi}}^E$ (e.g.
$\lambda{=}1$: ③ publishes $1.147841\times10^{-1}$ vs the gate's $1.148264\times10^{-1}$, a
$0.004\%$ difference). Both are certified upper bounds; the gate confirms its bracket
**overlaps** ③'s published bracket at all six shared $\lambda$ (fidelity to the inherited
certificate). Using the tighter ③ value only *strengthens* the sharpening; the table above
quotes the (looser, self-contained) gate output to keep every number reproducible from one
script.

### 2.2 Where the sharpening ends (honest boundary)

For $\lambda\gtrsim34.4$ the certified $h$-envelope is the strong-Weyl bracket
$\mathfrak g(2)\varpi^2-4\zeta(4)\varpi^4\le h\le\mathfrak g(2)\varpi^2$ (proven, Cor (ii)).
Its *upper* branch uses $G_{2,\mathrm{HI}}$, so inversion gives
$\varpi_{\mathrm{lo}}=\sqrt{(1/\lambda)/G_{2,\mathrm{HI}}}$, i.e.
$\varpi_{\mathrm{lo}}^2G_{2,\mathrm{HI}}\lambda=1$ **exactly** and
$F_{\mathrm{hi}}^E=(1/2\pi)\sqrt{G_{2,\mathrm{HI}}\lambda}\ge C_\infty\sqrt\lambda$. The gate
confirms this at $\lambda=45$: $\varpi_{\mathrm{lo}}^2G_{2,\mathrm{HI}}\lambda=1.0$ and
$\varpi_{\mathrm{lo}}^2G_{2,\mathrm{LO}}\lambda=0.99996<1$. So at large $\lambda$ the certified
envelope has **merged into the $C_\infty$ ceiling to within the $\mathfrak g(2)$ bracket
width**: certified sharpening is a strictly finite-$\lambda$ phenomenon (as it must be — the
ceiling is the $\lambda\to\infty$ limit). The representative table stops at $\lambda=30$, where
a genuine certified margin ($\times1.010$) still remains.

---

## 3. Honest scope — the envelope is **one-sided**

Theorem 2 says the Einstein spectrum **maximizes** $\mathfrak k(P,T)$ (hence $T_c$) at fixed
$\langle\omega^2\rangle$ and $\lambda$: this is an **upper** envelope, valid for all $P$. There
is **no matching nondegenerate lower envelope.** At fixed second moment,
$$
\inf_{P\in\mathcal P,\;\langle\omega^2\rangle\ \text{fixed}} T_c(P;\lambda)\;=\;0,
$$
approached only by **degenerate** spectra whose mass escapes to $\{0,\infty\}$. Concretely,
$P_\Omega=(1-p)\,\delta_{\varepsilon}+p\,\delta_{\Omega}$ with
$p=(\sigma^2-\varepsilon^2)/(\Omega^2-\varepsilon^2)$ keeps $\langle\omega^2\rangle=\sigma^2$
fixed while $\varepsilon\to0$, $\Omega\to\infty$ ($p\to0$); then $F_T(k)\to0$ for every $k$
(the near-zero atom contributes $[\![k]\!](\varepsilon/2\pi T)\to0$; the far atom carries
vanishing weight $p\to0$), so $\mathfrak k(P_\Omega,T)\to0$ and the crossing $T_c\to0$. The gate
(check group [D]) exhibits this at $\sigma^2=4$: with the Einstein maximum $h=1.336$, the
sequence gives $\mathfrak k=0.353,\,0.115,\,0.0358,\,0.00742,\,0.00234,\,6.66\times10^{-4}$ as
$\Omega=5,10,20,50,100,300$ — monotonically collapsing toward $0$, all below the Einstein
maximum. Every individual $P\in\mathcal P$ still has $T_c>0$ (Theorem 1(a)); the infimum is
simply not attained and the minimizing "spectrum" is degenerate. **We therefore claim the
upper envelope only.**

---

## 4. Proven / certificate map

| ingredient | rigor |
|---|---|
| Einstein extremality $\mathfrak k(P,T)\le h(\varpi_{\mathrm{rms}})$, $\varpi_{\mathrm{rms}}=\sqrt{\langle\omega^2\rangle}/2\pi T$ | **proven analytic** (Theorem 2, kernel-level Jensen + Lemma M) |
| moment match $\omega_E=\sqrt{\langle\omega^2\rangle}$ (second moment, exact) | **proven** (the Jensen average is over $x=\varpi^2$) |
| reduction $T_c(P;\lambda)\le T_c^E(\lambda;\sqrt{\langle\omega^2\rangle})=F^E(\lambda)\sqrt{\langle\omega^2\rangle}$ | **proven analytic** (Theorem 2 + monotone crossing, Theorem 1(b); homogeneity of $T_c^E$ in $\omega_E$) |
| certified upper bracket $F^E(\lambda)\le F_{\mathrm{hi}}^E(\lambda)$ | **inherited certificate-level** (deliverable ③: crossover finite-$N$ Kato–Temple; strong-Weyl proven; weak-W finite-$N$) |
| sharpening $F^E(\lambda)<C_\infty\sqrt\lambda$ at finite $\lambda$, $\to$ tight as $\lambda\to\infty$ | **proven analytic** (Corollary (ii): $r(\varpi)=h/\varpi^2<\mathfrak g(2)$, strictly $\downarrow$, $r(0^+)=\mathfrak g(2)$) |
| certified numeric test $\varpi_{\mathrm{lo}}^2G_{2,\mathrm{LO}}\lambda\ge1\Rightarrow F_{\mathrm{hi}}^E\le C_\infty\sqrt\lambda$ | **certified** (sound: only $\times$ and comparison on `mpf`; $\varpi_{\mathrm{lo}}$ from ③, $G_{2,\mathrm{LO}}$ CERT-B) |
| $C_\infty=(1/2\pi)\sqrt{\mathfrak g(2)}=0.1827262477\ldots$ | inherited (KAY; $\mathfrak g(2)$ CERT-B/CERT-D, analytic §4.1) |
| one-sided scope: $\inf_P T_c=0$ at fixed $\langle\omega^2\rangle$, degenerate | **proven analytic** (extreme-point / Jensen argument), numerically illustrated |
| deep-weak $F_{\mathrm{hi}}^E$ at $\lambda\lesssim0.34$ (the $\lambda=0.3$ row) | **asymptotic, NOT finite-$N$ certified** (Theorem W$'$ form, inherited from ③; flagged `weak-asympt*`) |

**Summary of rigor.** The *reduction* to the Einstein curve and the *sharpening* over the
$C_\infty$ ceiling are fully **analytic** (Theorems 1(b),1(c),2 and Corollary (ii), all proven
in `P4_main_theorems.md` with no computer-assisted step). The only certificate-level ingredient
is the numeric value of the Einstein curve $F_{\mathrm{hi}}^E(\lambda)$, inherited from ③; and
the single non-finite-$N$ (asymptotic) datum is the *upper* side of the deep-weak tail
$\lambda\lesssim0.34$ (the $\lambda=0.3$ row), explicitly flagged.

---

## 5. Source lines relied upon (`P4_main_theorems.md`)

- **§0 definitions**, lines 16–18: $\mathcal P$, $\langle\omega^2\rangle$,
  $[\![k]\!](\varpi)=\varpi^2/(k^2+\varpi^2)$; $\mathcal O[f]$ elementwise linear;
  $h(\varpi)=\lambda_1(\mathcal H(\varpi))$, $\mathfrak k(P,T)=\lambda_1(\mathcal O[F_T])$;
  $\lambda_1:=\max\mathrm{spec}$.
- **Theorem 1(c)**, line 55 (ceiling): $T_c<\tilde T_c=(1/2\pi)\sqrt{\mathfrak g(2)\langle\omega^2\rangle\lambda}$;
  stronger, $\mathfrak k(P,T)<\mathfrak g(2)\langle\omega^2\rangle/(2\pi T)^2$ for all $T$.
- **Theorem 1(b)**, line 54 (strict monotone, unique crossing) — underpins step (1) of the
  reduction.
- **Theorem 2**, lines 65–68 (Einstein extremality + the equivalent
  $T_c(\lambda,P)\le T_c^E(\lambda,\sqrt{\langle\omega^2\rangle})$ for all $\lambda$; the
  kernel-level Jensen proof, showing the matched frequency is $\omega_{\mathrm{rms}}$).
- **Corollary (i)**, line 72 ($h$ strictly increasing → unique $\varpi_\*$, monotone inversion)
  — underpins ③'s bracket and the well-definedness of $F^E$.
- **Corollary (ii)**, line 73 ($h(\varpi)<\mathfrak g(2)\varpi^2$; $r(\varpi)=h/\varpi^2$
  strictly decreasing, $r(0^+)=\mathfrak g(2)$) — the analytic engine of the sharpening.
- $C_\infty$ value and its Allen–Dynes / KAY attribution: lines 60–63, 330.

---

## 6. Reproduce

```
cd theory/derivations
# inherited certificate ③ (present in this worktree; commit 1184684):
OMP_NUM_THREADS=1 python p4_global_tc_certified.py test        # PASS 28  FAIL 0
# this deliverable:
OMP_NUM_THREADS=1 python p4_universal_tc_envelope.py           # PASS 66  FAIL 0  (~73 s)
```

The gate checks: [A] Einstein extremality on 6 genuine multi-frequency $P$ + 3 point-mass
equalities (the load-bearing Theorem 2); [B] the envelope table, the certified sharpening test
at 7 representative $\lambda$, monotonicity of the ratio toward 1, and interval-overlap
consistency with ③'s published $N{=}2000$ brackets; [C] the $C_\infty$ enclosure and the
strong-$\lambda$ boundary where sharpening ends; [D] the degenerate infimizing sequence showing
the envelope is one-sided.
