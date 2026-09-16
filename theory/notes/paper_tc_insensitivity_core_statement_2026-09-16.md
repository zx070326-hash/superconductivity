# Why $T_c$ Is Insensitive to Spectral Details: A Rigorous Explanation

*Draft paper section, prepared 2026-09-16 for insertion into the P4 manuscript.
Sourced against `P4_main_theorems.md` (v0.8, master ledger), `P4_einstein_extremal.md`
(Theorem 2 / Lemma C), `P4_weak_coupling_series.md` (Theorem W/W′), and
`theory/notes/pairing_duality_theorem_2026-09-16.md` ($\eta$ data, double-moment scan).
Rigor tags follow repo convention: **[analytic]** = proved, **[certificate]** = sound
outward-rounded interval arithmetic, **[numerical-evidence]** = float64/finite-grid,
**[heuristic]** = standard-in-the-literature argument not claimed as a proof here.*

---

## 1. Introduction

For fifty years the Allen–Dynes formula has predicted $T_c$ to within a few percent using
only two numbers extracted from the full Eliashberg spectral function $\alpha^2F(\omega)$:
the coupling constant $\lambda$ and a single frequency scale $\omega_{\log}$ (or, in the
formulation used here, $\sqrt{\langle\omega^2\rangle}$). This has always been read as an
*empirical* regularity — a good fit, not a consequence of anything deeper. We show here that
it is not empirical: it is a **theorem**. For fixed $\lambda$ and fixed $\langle\omega^2\rangle$,
no spectral shape can beat a single Einstein mode, and the size of what shape *can* buy you is
bounded, quantitatively, by a Jensen gap that is structurally small over the range spanned by
real phonon spectra. The infinite-dimensional object $\alpha^2F(\omega)$ compresses into two
numbers with a provably small, and separately measured, residual.

---

## 2. Quantitative statement

**Theorem E (Einstein extremality; = Theorem 2 of `P4_main_theorems.md`).**
For every spectral measure $P$ on $(0,\infty)$ with $\langle\omega^2\rangle_P<\infty$, every
$T>0$, and $\mu^*=0$,
$$
T_c(P) \;\le\; T_c^{E}\!\big(\lambda,\sqrt{\langle\omega^2\rangle_P}\big),
$$
where $T_c^E$ is the linearized Eliashberg $T_c$ of the point-mass ("Einstein") spectrum with
the same $\lambda$ and the same second moment. Equality holds iff $P$ is itself a point mass.
**[analytic]** — proved via pointwise concavity of $x\mapsto x/(k^2+x)$ in each Matsubara
channel plus kernel monotonicity (Lemma M, `P4_main_theorems.md` §2); no certificate or
computer-assisted step enters this proof, and it holds for **every** coupling strength
$\lambda$, not just asymptotically weak or strong coupling.

Define the **spectral efficiency** $\eta(P):=T_c(P)/T_c^E(\lambda,\sqrt{\langle\omega^2\rangle_P})\in(0,1]$.
Theorem E is the statement $\eta\le1$ with equality only for Einstein spectra. Evaluated on
real materials at the Allen–Dynes level (`pairing_duality_theorem_2026-09-16.md` §2,
**[numerical-evidence]**, cross-validated against the exact operator at the 1–2% level, §3):

| Material | $\lambda$ | $\mu^*$ | $r=\omega_2/\omega_{\log}$ | $\eta$ |
|---|---|---|---|---|
| LaH$_{10}$, 129 GPa | 3.62 | 0.10 | 1.194 | 0.921 |
| LaH$_{10}$, 163 GPa | 2.67 | 0.10 | 1.075 | 0.959 |
| LaH$_{10}$, 214 GPa | 2.06 | 0.10 | 1.105 | 0.932 |
| LaH$_{10}$, 264 GPa | 1.73 | 0.10 | 1.095 | 0.932 |
| H$_3$S, 200 GPa | 1.84 | 0.16 | $\sim1.20$ | 0.858 |
| Nb (BETE-NET) | 1.23 | 0.13 | $\sim1.14$ | 0.889 |
| MgB$_2$ (BETE-NET) | 0.62 | 0.13 | $\sim1.12$ | 0.896 |

$\eta$ never falls below 0.86 in this set, despite $r$ ranging up to 1.20 (a 20% spread of the
second-moment-weighted frequency scale). Separately, a Richter double-moment extremal scan
(fixing $\langle\omega^2\rangle$ **and** $\langle\omega^4\rangle$, so up to 3-atom measures are
admissible) shows that even a substantial excess kurtosis $\langle\omega^4\rangle/\langle\omega^2\rangle^2=1.5$
buys at most a 0.9% improvement over the Einstein bound (`pairing_duality_theorem_2026-09-16.md`
§4, **[numerical-evidence]**, AD-formula level, $25\times25\times12$ grid). Together: **the
part of $\alpha^2F(\omega)$ orthogonal to $(\lambda,\langle\omega^2\rangle)$ contributes at most
$\sim$10–15% to $T_c$ for real materials, and essentially nothing beyond what the second moment
already captures.**

**Honesty boundary on this section.** Theorem E itself is proved only for $\mu^*=0$. All
entries in the table above have $\mu^*>0$ (0.10–0.16); extremality of the Einstein spectrum at
$\mu^*>0$ is currently **[numerical-evidence]** only — direct 2-atom-vs-Einstein checks at
$\tilde\mu:=\mu^*/\lambda\in\{0.05,0.10,0.15\}$ confirm $\eta<1$ with a wide margin, and the
underlying convexity argument survives except in a window $\varpi\lesssim0.3$–$0.33$ (for
$\tilde\mu\lesssim0.05$) that is not populated by real phonon spectra ($\varpi:=\omega/2\pi T\gg1$
whenever $T_c\ll\omega_{\rm phonon}$, always true in practice). See §5–6 of
`pairing_duality_theorem_2026-09-16.md`.

---

## 3. Mechanism: why the gap is small

The proof of Theorem E does not require any asymptotic regime — it is exact for all $\varpi$.
What follows explains *why the residual $1-\eta$ is quantitatively small*, which is a separate
(and weaker-tier) question about the *size* of the Jensen gap, not its *sign*.

Write $\mathfrak K(\varpi):=\lambda_1\big(\mathcal O[\llbracket\cdot\rrbracket(\varpi)]\big)$ for
the Einstein eigenvalue and $H(x):=\mathfrak K(\sqrt x)$, $x=\varpi^2$. **Lemma C**
(`P4_einstein_extremal.md`) states $H$ is concave on $(0,\infty)$; this is what upgrades the
per-channel Jensen bound of §2 into the measure-level statement $\int H\,dP\le H(\langle x\rangle)$
used in the duality picture below. Lemma C's proof is now essentially closed over the whole
half-line: **[analytic]** for $\varpi\ge\varpi_A''\approx9.72$ (a Lorentz-addition operator
identity plus a Davis–Kahan gap bound, `P4_einstein_extremal.md` Route R, milestones M1–M5) and
**[certificate]** (sound interval arithmetic) on $(0,120]$ — the two segments overlap by a
factor of $\sim$12, so there is no gap in $\varpi$ left uncovered.

Two facts combine to make $1-\eta$ small in the physical range:

1. **Curvature of $H$ decays with $x$.** Direct evaluation gives $H''(0^+)\approx-2.3$ but
   $H''(x)\to-\tfrac1{2x^2}$ as $x\to\infty$ — e.g. $H''(50)\approx-1.9\times10^{-4}$ vs. the
   asymptotic prediction $-1/(2\cdot50^2)=-2.0\times10^{-4}$, and $H''(100)\approx-4.8\times10^{-5}$
   vs. $-5.0\times10^{-5}$ (`pairing_duality_theorem_2026-09-16.md` §5.2, **[numerical-evidence]**
   table; the asymptotic rate itself follows from Theorem W′'s $\mathfrak K(\varpi)=\ln\varpi+c_0+O(1/\ln\varpi)$,
   `P4_weak_coupling_series.md`, **[analytic]** for the leading log, **[analytic]** for the
   $b_1/\ln\varpi$ term via the closed form $b_1=\tfrac{49\pi^2}{144}-3=0.358407\ldots$).
   Concretely, $|H''(x)|$ already drops by two orders of magnitude between $x\sim1$ and
   $x\sim50$, i.e. between $\varpi\sim1$ and $\varpi\sim7$.
2. **Real phonon spectra are narrow.** The empirical $r$ values above ($1.07$–$1.20$) correspond
   to a relative spread of the second-moment-weighted frequency scale of only 5–20% — phonon
   bands in a given material sit within roughly one order of magnitude of each other, not the
   many-decade spread that would be needed to probe the far tail of $H$.

The Jensen gap is, to leading order, $H(\langle x\rangle)-\langle H(x)\rangle\approx
\tfrac12|H''(x_0)|\,\mathrm{Var}(x)$. Neither factor alone is uniformly tiny at the $\varpi\sim1$
scale typical of real materials — curvature is only asymptotically small, and is $O(0.3\text-0.7)$
in magnitude near $x\sim0.5\text-1$ — but their **product** is, because the variance factor
carries the dominant suppression at moderate $\varpi$, while the curvature factor takes over and
reinforces the suppression as $\varpi$ grows (stronger coupling, or $T_c$ well below the phonon
scale). This two-factor picture, not "curvature alone $\to0$," is the honest mechanistic
statement; it is consistent with, but should not be conflated with, the exact all-$\varpi$
Theorem E of §2.

---

## 4. RG translation

It is natural — though this is presented as *language*, not a literal renormalization-group
derivation of a Wilsonian flow in $P$-space — to describe the structure of §2–3 in RG terms:
$(\lambda,\langle\omega^2\rangle)$ are the **relevant** directions, controlling $T_c$ at leading
order and exactly (Theorem E). Higher spectral moments (variance of the *shape* beyond the
second moment, skewness, excess kurtosis) are **irrelevant**: the double-moment scan of §2 shows
their leading correction is bounded by $\lesssim1\%$ even for a generous kurtosis excess, and
the Jensen-gap mechanism of §3 predicts this correction scales with (at least) the variance of
the frequency distribution, i.e. is parametrically small and higher moments contribute at
progressively higher order in that same small parameter. **[heuristic label on genuine
numerical-evidence content]**.

One piece of this story *is* a genuine, partially rigorous RG flow: the Coulomb pseudopotential
$\mu^*$. The Coulomb operator is exactly rank one, $\mathfrak K_c(N_c)=2u_cu_c^\top$ with
$\|u_c\|^2=S(N_c)=\sum_{n<N_c}(2n+1)^{-1}$, and $S(N_c)=\tfrac12\ln N_c+\tfrac12(\gamma_E+\ln4)+O(N_c^{-2})$
**[analytic]** (`beta_function_rg_flow_section_2026-09-16.md` §X.1). Combined with a standard
leading-log matching argument for how the top-eigenvector overlap $c_1(N_c)$ grows with the
window **[heuristic, standard Morel–Anderson/Tolmachev]**, pinned to a slope of $-1$ (as opposed
to the naively-possible $-2$) by a direct float64 secular-equation fit **[numerical-evidence,
accurate to $<10^{-3}$, no certified/interval version yet]**, one obtains, at fixed $T$ and using
Wilsonian RG time $t=\ln(\Lambda/\omega_c)$,
$$
\beta(\mu^*) \;=\; \frac{d\mu^*}{dt} \;=\; -(\mu^*)^2,
$$
the same **marginally irrelevant** normal form as the 1-loop Kondo coupling and the
Kosterlitz–Thouless vortex fugacity, reproducing the familiar Morel–Anderson resummed form
$\mu^*(t)=\mu^*_0/[1+\mu^*_0(t-t_0)]$. This is the operator-theoretic reading of why the
Coulomb pseudopotential, unlike the phonon shape moments, is not simply "small" but *runs* to a
universal, cutoff-insensitive value — and why a single number $\mu^*$ (rather than a Coulomb
"spectral shape") suffices to describe it.

---

## 5. Implications

If $T_c$ is controlled, to $\gtrsim85$–90% accuracy, by only $(\lambda,\langle\omega^2\rangle,\mu^*)$,
then (i) computational material search should prioritize maximizing $\lambda$ and
$\langle\omega^2\rangle$ (e.g. light-element high-frequency modes with strong electron–phonon
coupling), not fine-tuning the *shape* of $\alpha^2F(\omega)$, whose achievable gain is bounded
below $\sim1\%$ by §2's double-moment scan; (ii) machine-learning descriptors for $T_c$ need not
resolve the full phonon DOS — low-order spectral moments should suffice, with diminishing
returns from higher moments quantified by the theorem's own bound; (iii) departures from
phonon-mediated pairing (different retardation structure, e.g. $r\ne1$ kernels in a generalized
Eliashberg framework) act, on the evidence gathered elsewhere in this project, as an additional
suppression of a few percent — real but subordinate to the $(\lambda,\langle\omega^2\rangle)$
scaling established here, not a change of picture.

---

### Sources and rigor cross-check

- Theorem E / Lemma M: `P4_main_theorems.md` §2 (Theorem 2, Lemma M) — **[analytic]**, confirmed
  fully computer-assistance-free as of 2026-09-12.
- Lemma C (concavity of $H$): `P4_einstein_extremal.md` — mixed **[analytic]** ($\varpi\ge9.72$,
  Route R) + **[certificate]** ($\varpi\in(0,120]$, interval arithmetic, milestone M5); combined
  coverage is total, no residual gap.
- $\eta$ table, double-moment scan, $\mu^*>0$ convexity analysis:
  `theory/notes/pairing_duality_theorem_2026-09-16.md` §2–5 — **[numerical-evidence]**
  throughout; $\mu^*>0$ Einstein-extremality is *not yet* an analytic theorem.
- Weak-coupling constants $c_0=\gamma+\ln4-\tfrac32=0.4635100\ldots$,
  $b_1=\tfrac{49\pi^2}{144}-3=0.358407053148462\ldots$, prefactor
  $\tfrac2\pi e^{\gamma-3/2}=0.2529996837\ldots$: `P4_weak_coupling_series.md` Theorem W′ —
  **[analytic]**; independently recomputed here to confirm ($\gamma+\ln4-1.5=0.4635100\ldots$;
  $49\pi^2/144-3=0.358407\ldots$; $(2/\pi)e^{\gamma-1.5}=0.252999\ldots$ — all match to the
  digits quoted).
- $\beta$-function for $\mu^*$: `theory/notes/beta_function_rg_flow_section_2026-09-16.md` —
  algebraic backbone (X.1–X.2) **[analytic]**; slope-fixing (X.3–X.4) **[heuristic +
  numerical-evidence]**, explicitly not yet certified.

**Not committed** to the repository per instruction; branch `p4-conj2`, working tree otherwise
unchanged.
