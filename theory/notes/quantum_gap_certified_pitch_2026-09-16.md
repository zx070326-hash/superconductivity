# Certified spectral gaps for small quantum spin chains via interval arithmetic

**Status: methods/tools note draft, 2026-09-16. Extends the 2026-09-16 probe
(`theory/notes/probe_quantum_algorithm_gap_2026-09-16.md`) by hardening the previously
float64 Rayleigh-quotient side to `mpmath.iv` interval arithmetic — both bound directions
are now certificate-level. Not a merge candidate for `p4-conj2`; standalone note.**
**Strictness: every number in the tables below is `certificate`-level (interval Rayleigh
lower bound + extended-precision-Cholesky upper bound with an explicit, tested pad).
Nothing here is `analytic`.**

## 1. Motivation

Adiabatic quantum computing (AQC) and QAOA both treat the instantaneous spectral gap as
the central complexity-limiting quantity, yet gap values reported in the literature and in
simulators are almost always plain floating-point numbers — not certified lower bounds.
Near a level crossing (exactly where AQC runtime bounds are tightest and QAOA landscapes
are hardest), a float64 estimate cannot itself distinguish "gap is small" from "gap is a
rounding artifact of a near-degenerate diagonalization." Certified methods that exist in
the literature target the opposite regime: SDP hierarchies (Knabe-family, e.g. Anshu–Arad
et al., arXiv:2411.03680) certify *bulk* (N→∞) gaps of frustration-free translation-invariant
chains; hand-derived analytic bounds exist for specific adiabatic constructions (e.g.
Cesa–Herrman, PRA 101, 042302). Neither is a general, sound, *pointwise finite-N* numerical
tool. We report such a tool, repurposed with no modification from a spectral-gap certificate
built for a strong-coupling superconductivity problem in this repository.

## 2. Method

Convention: `A := -H`, so `λ₁(A) = -E₀`, `λ₂(A) = -E₁`, gap `Δ = λ₁(A) - λ₂(A)`. Both bound
directions are audited per this repo's one-sided-bound-direction rule (a lower estimate of
an upper bound only ever supports "certified", never "falsified"):

- **`λ₂(A) ≤ c`** — rank-1 Cholesky-deflation certificate, ported verbatim from
  `theory/derivations/p4_interchange_uniformity.py :: alpha_certificate_H`. Deflate the
  float64 top eigenvector `u` (`P = I - uuᵀ`), certify `c·I - PAP ≻ 0` via `mp.cholesky` at
  `dps=30`, `c` = float64 top eigenvalue of `PAP` plus an explicit pad. Because `A` here is
  the *whole* Hamiltonian (no truncation), this is exact modulo (a) the eigenvector's
  float64 residual and (b) `A` itself being built from float64-rounded model parameters
  rather than the exact decimal value — both effects are bounded and covered by the pad
  (see §3, pad-floor scan).
- **`λ₁(A) ≥ lb`** (**new in this note**) — interval Rayleigh quotient
  `⟨v|A|v⟩/⟨v|v⟩` computed with `mpmath.iv` outward rounding. `A`'s matrix elements are built
  directly from the *decimal-string* parameter (`iv.mpf("0.7")`, enclosing the intended
  exact value, not a float64 rounding of it) via a direct bit-operation action on the state
  vector — no dense interval matrix is ever materialized, giving `O(dim·N)` instead of
  `O(dim²)` for this step. `v` is the float64 top eigenvector cast into `iv.mpf`
  component-wise, which is *lossless* (a float64 value is an exact dyadic rational, so
  `mpmath`'s base-2 mantissa represents it with zero extra rounding). By the Rayleigh–Ritz
  variational principle `λ₁(A) ≥ ⟨v|A|v⟩/⟨v|v⟩` for **any** real `v`, and since every interval
  operation outward-rounds, the left endpoint `lb := R.a` is a sound certificate-level lower
  bound on `λ₁(A)`. This closes the gap flagged in the original probe, where only the
  `λ₂` side was certificate-grade.

`gap_cert_lb := lb - c` is then a fully certificate-level lower bound on the true gap.

**Pad choice.** The pad must cover (a) the eigenvector-residual term (established at
`~1e-9` in the original P4 review) and (b) a Weyl-inequality bound on the parameter-rounding
mismatch between the float64-built `A` (used for the Cholesky side) and the exact-decimal
`A` (used for the Rayleigh side): `|λ₂(A_float) - λ₂(A_exact)| ≤ ‖A_float - A_exact‖₂ ≲
√(nnz)·10⁻¹⁶ ≲ 10⁻¹⁴` for the sizes tested here. A pad-floor scan
(`theory/notes/_scratch_ir/quantum_gap_certified_padscan.txt`) confirmed `mp.cholesky`
still **passes** (i.e. the certificate is still sound — Cholesky success is itself the
proof, not an estimate) down to pad `1e-12` on the near-critical instances tested; we adopt
`pad = 1e-11` throughout as a safety margin above that observed floor.

## 3. Results

Script: `theory/notes/_scratch_ir/quantum_gap_certified.py` (full run log:
`theory/notes/_scratch_ir/quantum_gap_certified_output.txt`, both TFIM and XXZ, N=6,8).
All 32 (model, N, parameter) cases below: **Cholesky certificate PASSes at every point**
(0/32 fails), including gaps as small as `~1e-8`–`1e-11`.

**TFIM** (periodic, `H(s) = -(1-s)ΣXᵢ - sΣZᵢZᵢ₊₁`):

| N | s | E₀ (float64) | gap (float64) | lb(λ₁) [certificate] | c ≥ λ₂ [certificate] | gap certified lb | slack |
|---|-----|----------|-----------|-----------|-----------|-----------|-------|
| 6 | 0.50 | −3.86370 | 0.131652 | 3.863703 | 3.732051 | 0.131652 | 1.0e-11 |
| 6 | 0.70 | −4.39617 | 0.001975 | 4.396166 | 4.394191 | 0.001975 | 1.0e-11 |
| 6 | 0.90 | −5.41668 | 8.29e-07 | 5.416680 | 5.416679 | 8.29e-07 | 1.0e-11 |
| 6 | 0.95 | −5.70395 | 9.93e-09 | 5.703948 | 5.703948 | 9.92e-09 | 1.0e-11 |
| 6 | 0.99 | −5.94015 | 5.15e-13 | 5.940152 | 5.940152 | **inconclusive** | (below pad floor) |
| 6 | 1.00 | −6.00000 | 0.0 (exact) | 6.000000 | 6.000000 | **inconclusive** | exact/near-exact degeneracy |
| 6 | 1.01 | −6.06015 | 4.68e-13 | 6.060149 | 6.060149 | **inconclusive** | (below pad floor) |
| 6 | 1.05 | −6.30357 | 6.02e-09 | 6.303572 | 6.303572 | 5.02e-09 | 1.0e-11 |
| 6 | 1.10 | −6.61364 | 3.05e-07 | 6.613644 | 6.613643 | 3.05e-07 | 1.0e-11 |
| 6 | 1.50 | −9.25227 | 0.000967 | 9.252270 | 9.251303 | 0.000967 | 1.0e-11 |
| 8 | 0.50 | −5.12583 | 0.098491 | 5.125831 | 5.027339 | 0.098491 | 1.0e-11 |
| 8 | 0.70 | −5.86039 | 0.000307 | 5.860395 | 5.860087 | 0.000307 | 1.0e-11 |
| 8 | 0.90 | −7.22224 | 8.71e-09 | 7.222239 | 7.222239 | 8.70e-09 | 1.0e-11 |
| 8 | 0.95 | −7.60526 | 2.34e-11 | 7.605264 | 7.605264 | 1.34e-11 | 1.0e-11 |
| 8 | 0.99 | −7.92020 | 3.55e-15 | 7.920202 | 7.920202 | **inconclusive** | (below pad floor) |
| 8 | 1.00 | −8.00000 | 0.0 (exact) | 8.000000 | 8.000000 | **inconclusive** | exact/near-exact degeneracy |
| 8 | 1.01 | −8.08020 | 5.33e-15 | 8.080198 | 8.080198 | **inconclusive** | (below pad floor) |
| 8 | 1.05 | −8.40476 | 1.16e-11 | 8.404763 | 8.404763 | 1.6e-12 | 1.0e-11 |
| 8 | 1.10 | −8.81819 | 2.14e-09 | 8.818191 | 8.818191 | 2.13e-09 | 1.0e-11 |
| 8 | 1.50 | −12.33576 | 9.10e-05 | 12.335761 | 12.335670 | 9.10e-05 | 1.0e-11 |

**XXZ** (periodic, `H = Σᵢ [SxᵢSxᵢ₊₁ + SyᵢSyᵢ₊₁ + Δ·SzᵢSzᵢ₊₁]`):

| N | Δ | E₀ (float64) | gap (float64) | lb(λ₁) [certificate] | c ≥ λ₂ [certificate] | gap certified lb | slack |
|---|-----|----------|-----------|-----------|-----------|-----------|-------|
| 6 | 0.0 | −2.00000 | 0.267949 | 2.000000 | 1.732051 | 0.267949 | 1.0e-11 |
| 6 | 0.5 | −2.36803 | 0.454018 | 2.368034 | 1.914016 | 0.454018 | 1.0e-11 |
| 6 | 0.9 | −2.71063 | 0.634683 | 2.710628 | 2.075945 | 0.634683 | 1.0e-11 |
| 6 | 1.0 | −2.80278 | 0.684742 | 2.802776 | 2.118034 | 0.684742 | 1.0e-11 |
| 6 | 1.5 | −3.30085 | 0.529836 | 3.300849 | 2.771013 | 0.529836 | 1.0e-11 |
| 6 | 2.0 | −3.85577 | 0.406283 | 3.855773 | 3.449490 | 0.406283 | 1.0e-11 |
| 8 | 0.0 | −2.61313 | 0.198912 | 2.613126 | 2.414214 | 0.198912 | 1.0e-11 |
| 8 | 0.5 | −3.08699 | 0.340015 | 3.086994 | 2.746980 | 0.340015 | 1.0e-11 |
| 8 | 0.9 | −3.53117 | 0.482099 | 3.531165 | 3.049066 | 0.482099 | 1.0e-11 |
| 8 | 1.0 | −3.65109 | 0.522674 | 3.651093 | 3.128419 | 0.522674 | 1.0e-11 |
| 8 | 1.5 | −4.30355 | 0.375722 | 4.303545 | 3.927823 | 0.375722 | 1.0e-11 |
| 8 | 2.0 | −5.03943 | 0.258764 | 5.039429 | 4.780665 | 0.258764 | 1.0e-11 |

XXZ has no near-degenerate case in this parameter range, so all 12 rows certify cleanly
with slack pinned at the pad (`1.0e-11`).

**Wall-clock.** N=6 (dim 64): <0.3 s/case. N=8 (dim 256): 5.5–14 s/case, dominated by the
`O(dim³)` pure-Python `mp.cholesky` at `dps=30` (the interval-Rayleigh step is negligible by
comparison, `O(dim·N)`). This is the same `dim³`-not-`N³` bottleneck already diagnosed in
this repo's P4 breakpoint-4 audit, transplanted to a different Hamiltonian family.

## 4. Discussion

**The near-`s=1` TFIM cases are the interesting honest result, not a bug.** At `s ∈
{0.99, 1.0, 1.01}` the true float64 gap for this finite periodic chain is `≲1e-13` or
exactly machine-zero — no pad shrinkage certifies positivity there (a pad-floor scan down
to `1e-12` still leaves these three inconclusive; see
`quantum_gap_certified_padscan.txt`). This is very likely an exact (or near-exact)
finite-size symmetry-protected degeneracy at the self-dual point, not a precision artifact
of the tool — the certificate correctly reports "cannot conclude" rather than either
falsely asserting a gap or silently returning a meaningless float64 difference of two
numbers at the rounding floor. A symmetry-sector-resolved version of the certificate (block
per momentum/parity sector) would likely resolve this, but that is Hamiltonian-specific
follow-up work, not a two-line fix.

**Scope limit: `dim³` wall.** `mp.cholesky` at `dps=30` is exponential in qubit count with a
non-trivial constant; N=8 (dim=256) already takes up to ~14 s/case, and N≈10–12 is expected
to become impractical without sparsity-aware or subspace-truncated Cholesky. This tool is
therefore a **small-system (≲10-qubit) validator**, not a route to new asymptotic
gap-scaling results — that remains the province of the SDP/analytic literature (§1), which
targets N→∞ and is not competing with this tool's niche.

**Both bound directions are now certificate-grade** (this note's contribution over the
2026-09-16 probe): the Rayleigh side is interval-hardened via `mpmath.iv`, matching the
Cholesky side's rigor, and the cross-term between "float64-rounded `A`" (Cholesky side) and
"exact-decimal `A`" (Rayleigh side) is explicitly bounded via Weyl's inequality rather than
assumed negligible.

## 5. Code availability

- `theory/notes/_scratch_ir/quantum_gap_certified.py` — full script (both bound directions
  hardened; TFIM + XXZ builders; pad rationale in module docstring).
- `theory/notes/_scratch_ir/quantum_gap_certified_output.txt` — full run log (32 cases).
- `theory/notes/_scratch_ir/quantum_gap_certified_padscan.txt` — pad-floor sensitivity scan.
- Predecessor probe (float64-Rayleigh only): `theory/notes/_scratch_ir/probe_qc_gap_cert.py`,
  `probe_qc_gap_cert_xxz.py`, write-up `theory/notes/probe_quantum_algorithm_gap_2026-09-16.md`.
- Reused core primitive: `theory/derivations/p4_interchange_uniformity.py :: alpha_certificate_H`.

All scripts are scratch/uncommitted per repo convention; none of this is merged into
`p4-conj2` or any other theory branch.

## 6. Target and framing

**Target:** *Quantum* (open access) technical note, or PRA brief report — same calibrated
framing recommended in the predecessor probe: a **sound cross-check / benchmarking tool for
NISQ-scale (≲10-qubit) numerics**, not a new AQC/QAOA asymptotic-scaling result. The
near-critical TFIM findings (§4) — a certificate that correctly refuses to conclude exactly
where a finite-size symmetry degeneracy sits — are a genuine, citable illustration of what a
*sound* (vs. merely high-precision) numerical gap bound buys you, and are a good hook for
the note's introduction.
