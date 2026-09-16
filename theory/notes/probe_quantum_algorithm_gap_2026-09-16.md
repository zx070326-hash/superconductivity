# Probe: Certified spectral gap bounds for quantum algorithms via interval arithmetic

**Status: PROBE / proof-of-concept, 2026-09-16. Time-boxed ~2.5h. Not a merge candidate for `p4-conj2`; standalone exploratory note.**
**Strictness of claims below: feasibility test is `numerical-evidence` + a partially-certificate-level sub-step (flagged inline). Nothing here is `analytic`.**

## 1. What we already have (reused, not rebuilt)

`theory/derivations/p4_interchange_uniformity.py :: alpha_certificate_H` (lines 725–780) implements a
**rank-1 Cholesky-deflation certificate**: given a symmetric matrix `A`, a float64 top eigenvector `u`,
and a candidate ceiling `c` (float64 top eigenvalue of the deflated `PAP`, `P = I − uuᵗ`, plus a fixed
pad), it runs `mp.cholesky(c·I − PAP)` at `dps=30`. Success certifies `λ₁(PAP) < c`, and a
residual-of-inexact-eigenvector argument (documented in the docstring, re-derived by the P4 reviewer on
2026-09-10) then gives `λ₂(A) ≤ c`. In the P4 use case, `A` is only the *head block* of a much larger
operator, so the certificate also needs rigorous Weyl/Schur tail-norm bounds (`‖B‖`, `‖D‖`) to control
the truncation — those tail bounds are P4-specific (closed-form Brillouin-function decay) and do **not**
port to a generic Hamiltonian.

**What ports directly, with no modification:** the rank-1 deflation + Cholesky step itself, applied to
the *entire* Hamiltonian matrix of a small quantum system (no truncation needed if the full Hilbert space
is small enough to diagonalize/Cholesky-test exactly).

## 2. Literature scan (30 min, WebSearch — not exhaustive)

Existing rigorous/certified gap-bound methods in the AQC/QAOA/spin-system literature cluster into three
families, **none of which is "interval arithmetic on a finite Hamiltonian matrix along an adiabatic
path"**:

- **SDP hierarchies for frustration-free systems** (Knabe bound and generalizations) — certify *bulk*
  (thermodynamic-limit) gaps of translation-invariant frustration-free spin chains via a hierarchy of
  semidefinite programs. [Anshu–Arad et al., "A Hierarchy of Spectral Gap Certificates for
  Frustration-Free Spin Systems," arXiv:2411.03680](https://arxiv.org/html/2411.03680). This targets
  N→∞, not finite small circuits, and is SDP-based rather than interval-arithmetic-based.
- **Certified upper bounds for bulk gaps** — a parallel convergent-upper-bound family; explicitly flags
  *lower*-bound certification of bulk gaps as an open direction.
  [arXiv:2606.03836](https://arxiv.org/html/2606.03836).
- **Analytic/asymptotic gap bounds for structured adiabatic Hamiltonians** — e.g. circuit-to-adiabatic
  constructions with gap ≥ π²/[8(L+1)²] (Cesa–Herrman–style improved bounds,
  [Phys. Rev. A 101, 042302](https://link.aps.org/doi/10.1103/PhysRevA.101.042302),
  [arXiv:1906.05233](https://arxiv.org/abs/1906.05233)); optimized-schedule papers deriving explicit
  piecewise analytic gap lower bounds over the adiabatic path
  ([arXiv:2411.05736](https://arxiv.org/pdf/2411.05736), Quantum 2025). These are hand-derived analytic
  inequalities for specific constructions, not a general sound-numerics tool.
- **Quantum speed limits as *necessary*-condition lower bounds on adiabatic runtime** (not gap bounds
  per se) — [Phys. Rev. Research 5, 033175](https://doi.org/10.1103/PhysRevResearch.5.033175).
- General QAOA-convergence work invokes gap-condition-free versions of the adiabatic theorem rather than
  computing certified gaps numerically ([IOP, "Elementary proof of QAOA
  convergence"](https://iopscience.iop.org/article/10.1088/1367-2630/ad59bb)).

**Conclusion of the scan:** I did not find prior work applying sound (rounding-error-controlled) interval
or arbitrary-precision arithmetic to certify the *pointwise, finite-size* spectral gap of a concrete small
Hamiltonian along an adiabatic schedule. The closest relative (Knabe-family SDP) solves a different
problem (bulk/asymptotic, frustration-free only). This is a genuine, narrow gap — but see §4 for why it's
narrow rather than a home run.

## 3. Feasibility test (case selection + result)

**Case 1 — transverse-field Ising chain (TFIM)**, periodic, `N` qubits:
`H(s) = −(1−s) Σᵢ Xᵢ − s Σᵢ ZᵢZᵢ₊₁`, standard AQC interpolation `s ∈ [0,1]`, minimum gap near the
critical point `s ≈ 0.5`.

**Case 2 — Heisenberg XXZ chain**, periodic, anisotropy `Δ = 0.5`:
`H = Σᵢ (SxᵢSxᵢ₊₁ + SyᵢSyᵢ₊₁ + Δ·SzᵢSzᵢ₊₁)`.

Both are textbook, exactly-diagonalizable-at-small-N, gap-is-the-bottleneck systems with abundant
reference data (TFIM finite-size gap scaling is exactly solvable via Jordan–Wigner; XXZ via Bethe ansatz)
— satisfying the brief's requirement of "Hamiltonian well-defined, gap is the key bottleneck, reference
data exists."

**Method**: set `A := −H(s)` so the ground-state gap `Δ = E₁ − E₀ = λ₁(A) − λ₂(A)` (top minus second
eigenvalue of `A`). Bound directions, audited explicitly per the repo's one-sided-bound rule:
- `λ₁(A) ≥ Rayleigh(A, u)` for **any** real vector `u` — sound in exact arithmetic; here evaluated in
  float64, so this half is **not yet interval-hardened** (flagged, see §4).
- `λ₂(A) ≤ c` via the verbatim P4 Cholesky-deflation test (`mp.cholesky` at `dps=30` on
  `c·I − PAP`, `P = I − uuᵗ`) — this half **is** certificate-level in the same sense as the P4 script,
  since no truncation is involved (the whole Hilbert space fits in the matrix).
- `Δ ≥ Rayleigh(A,u) − c`.

Script: `theory/notes/_scratch_ir/probe_qc_gap_cert.py` (TFIM, N=6,8) and
`theory/notes/_scratch_ir/probe_qc_gap_cert_xxz.py` (XXZ, N=6,8,10) — both scratch, uncommitted, reuse
`op`/`cholesky_deflation_upper_bound` verbatim-adapted from `alpha_certificate_H`.

**Results, TFIM** (`dim = 2^N`):

| N | dim | s | E₀ (float64) | gap (float64) | Cholesky | certified `c ≥ λ₂` | gap lower bound (cert) | slack |
|---|-----|-----|---------|--------|----|---------|---------|-------|
| 6 | 64 | 0.30 | −4.39617 | 0.801975 | PASS | 3.594191 | 0.801975 | 1.0e-9 |
| 6 | 64 | 0.50 | −3.86370 | 0.131652 | PASS | 3.732051 | 0.131652 | 1.0e-9 |
| 6 | 64 | 0.70 | −4.39617 | 0.001975 | PASS | 4.394191 | 0.001975 | 1.0e-9 |
| 8 | 256 | 0.30 | −5.86039 | 0.800307 | PASS | 5.060087 | 0.800307 | 1.0e-9 |
| 8 | 256 | 0.50 | −5.12583 | 0.098491 | PASS | 5.027339 | 0.098491 | 1.0e-9 |
| 8 | 256 | 0.70 | −5.86039 | 0.000307 | PASS | 5.860087 | 0.000307 | 1.0e-9 |

The certificate holds (Cholesky PASS) at **every** tested point, including `s=0.70, N=8` where the true
gap is only `3.07e-4` — i.e. the method survives near-degenerate points, which is exactly the regime
adiabatic algorithms care about. The slack between the certified bound and the float64-reference gap is
`1.0e-9` at every point — i.e. essentially the fixed pad we chose, because rank-1 deflation on a
*non-truncated* full matrix is exact up to the eigenvector's float64 residual. This is a stronger (tighter)
result than the P4 case, where the pad also has to absorb genuine Weyl tail terms `‖B‖+‖D‖` from
truncating an infinite-dimensional operator.

XXZ case (N=6, 8) confirms the same pattern; N=10 (dim=1024) is the point where the pure-Python
`mp.cholesky` at `dps=30` (`O(dim³)` in Python-level arithmetic) starts to become the wall — see §4.

## 4. Feasibility verdict and honest limitations

**Works, cheaply, for small systems.** The rank-1 Cholesky-deflation half of the P4 machine ports
verbatim to any real-symmetric Hamiltonian and gives a sound (mod the flagged Rayleigh caveat) gap lower
bound with negligible slack, for `dim` up to a few hundred (N ≲ 8–9 qubits) in seconds.

**Two limitations that bound the pitch's ambition:**

1. **`dim³` wall, not `N³` wall.** `mp.cholesky` at `dps=30` is pure-Python `O(dim³)`; `dim = 2^N` means
   this is exponential in qubit count with a bad constant — the exact same "Na³ is the bottleneck, not
   the exponent" diagnosis as the P4 breakpoint-4 audit (`c₁` elasticity finding), just transplanted to a
   different N. N=10 (dim=1024) is already visibly slow in a first background run; N=12 (dim=4096) is
   very unlikely to finish in reasonable time without either (a) exploiting sparsity/symmetry to avoid
   materializing the dense `dim×dim` matrix, or (b) restricting the certificate to a *low-dimensional
   relevant subspace* the way P4 does with its head-block `A` plus Weyl tail bounds — but building sound,
   problem-agnostic tail-norm bounds for a generic spin Hamiltonian (no closed-form Brillouin-function
   decay to lean on) is a **nontrivial, model-specific research problem in itself**, not a two-line port.
   So this tool as it stands is a **small-system (NISQ-scale, ≲10-qubit) validator**, not a route to new
   asymptotic gap-scaling theorems for large N — which is what the AQC/QAOA community mostly cares about.

2. **Only one of the two bound directions is certificate-grade today.** The Rayleigh-quotient lower bound
   on `λ₁(A)` was evaluated in plain float64 in this probe. Hardening it to interval arithmetic
   (`mpmath.iv` dot products with directed rounding) is straightforward — same order of effort as the
   existing Cholesky step — but was not done in this 2.5h probe; the number reported above is
   `numerical-evidence`, not full certificate, on that half.

**Bottom line: technically feasible, real niche, modest in scope.** The honest value proposition is not
"solve AQC/QAOA scaling" (that requires exactly the SDP/analytic machinery already in the literature,
which targets N→∞) — it's **a sound cross-check / benchmarking tool for small (≲10 qubit) instances**,
useful for (a) validating that a claimed numerically-estimated gap from a heuristic simulator isn't a
rounding artifact near a level crossing, and (b) teaching/verification use where "this number is
*provably* a lower bound, not just a float64 estimate" has standalone value. That is a real, publishable,
but narrow contribution — closer to a tools/methods note than a breakthrough.

## 5. Two-page pitch draft

**Title:** *Certified spectral gap bounds for small quantum systems via interval/rational arithmetic:
a portable Cholesky-deflation certificate*

**Target:** methods/tools section of Physical Review A, or a short Quantum (open-access) paper. Given the
narrow, honest scope established in §4, a "Letter"/"Technical note" framing fits better than a full
research article claiming a new capability for AQC theory.

**Abstract (draft):**
> Certifying that a numerically estimated spectral gap of a quantum many-body Hamiltonian is a *sound*
> lower bound — not merely a floating-point estimate — is rarely done in practice, even though the
> adiabatic and variational quantum algorithm literature treats the gap as the central complexity-limiting
> quantity. We repurpose a rank-1 Cholesky-deflation certificate, originally built for a spectral-gap
> problem in strong-coupling superconductivity theory, to produce rigorous (interval-arithmetic-backed)
> lower bounds on the ground-state gap of small (≲10-qubit) spin Hamiltonians — transverse-field Ising
> and Heisenberg XXZ chains — across an adiabatic interpolation path, including near level crossings where
> the gap nearly closes. The method is exact up to a fixed, explicit rounding pad and requires no
> problem-specific analytic input when the full Hilbert space is tractable to diagonalize. We report
> negligible slack (~1e-9) against float64 reference values on 6–8 qubit instances and characterize the
> `dim³` wall (N≳10–12 qubits) that limits direct scaling, discussing what would be needed (subspace
> truncation with sound tail-norm bounds) to push beyond it. We position this as a lightweight
> cross-validation tool for NISQ-scale numerics rather than a new asymptotic gap-scaling result.

**Core sections:**
1. Motivation: gap-is-the-bottleneck in AQC/QAOA; existing certified methods target the bulk/asymptotic
   regime (SDP hierarchies, e.g. Knabe-family) or hand-derived analytic bounds for specific constructions;
   no general sound *pointwise, finite-N* numerical tool exists.
2. Method: rank-1 deflation + directed-rounding Cholesky test (full derivation, both bound directions
   hardened — the Rayleigh side needs `mpmath.iv` before publication, see §4.2 above).
3. Results: TFIM and XXZ, N=6,8(,10), gap vs. adiabatic parameter `s`, slack table, wall-clock scaling
   plot showing the `dim³` cliff.
4. Discussion: honest scope — benchmarking/validation tool, explicit non-claim about asymptotic gap
   scaling; sketch of what a subspace-truncation extension would require (sound tail bound analogous to
   the Brillouin-function decay used in the superconductivity application, but model-specific for spin
   chains — likely the actual research contribution if this is pursued further).

**Novelty claim (calibrated):** first demonstration of a general-purpose, portable, sound (not just
high-precision) numerical gap certificate for finite quantum spin Hamiltonians, validated on standard
benchmark models near their hardest (near-degenerate) points.

**Next steps if greenlit:**
- Harden the Rayleigh side with `mpmath.iv` (few hours).
- Push to N=12–14 by replacing dense `mp.cholesky` with a sparse/structured directed-rounding scheme, or
  by testing the deflation on a Krylov/Lanczos-derived low-dimensional candidate subspace with an
  explicit, derived (not assumed) tail-norm bound — this is the actual open research question, and is
  where a Google Quantum AI / IBM Quantum collaboration would matter (they have the large-N heuristic gap
  estimates that a sound cross-check would be valuable against).
- Do **not** oversell to potential collaborators as an AQC scaling-theory result; frame explicitly as a
  validation/methods tool. Overselling this internally-narrow result risks the same class of problem this
  repo's rigor discipline exists to prevent.

## 6. Recommendation to team-lead

**Feasible, worth a short methods note; not worth an outward "collaborate with Google/IBM" pitch yet.**
The technique genuinely ports with no modification for small systems and gives essentially-free sound
bounds there — that part of the probe is a clean win. But the ambitious framing in the original probe
brief ("first sound + automated tool for AQC", collaboration with Google/IBM) overstates what a
≲10-qubit dense-Cholesky certificate can currently deliver against a field whose real interest is
large-N asymptotic gap scaling, which this method does not address without solving a new (nontrivial,
model-specific) subspace-truncation problem first. Recommend: (a) keep this as an internal proof-of-concept,
(b) if pursued, scope the paper as the narrow tools/validation contribution described above, (c) only raise
it with external quantum-computing groups after either the subspace-truncation extension exists or the
target is explicitly reframed as "sound cross-check for NISQ-scale numerics," not "certified AQC theory."
