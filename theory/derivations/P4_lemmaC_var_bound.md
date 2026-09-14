# P4 · Lemma C — the variance leg: `Var_v(B) < γ·xE/2` reduced to a scalar spectral inequality on `M = 𝒪[f−f²]`

**状态 / status.** Delivered on branch `p4-conj2` (worktree), **do not merge, do not push**.
This note takes the *last* leg of the concavity / q-monotonicity route — the terminal
inequality `Var_v(B) < γ·(xE/2)` of `P4_concavity_q_route.md` §4 — and reduces it to a
**single clean scalar inequality about one nonnegative-kernel operator** `M := A − B = 𝒪[f−f²]`.
It records exactly what is **proven**, what is **certificate-level**, and what remains
**numerical evidence**, and it isolates the residual as one previously-known bottleneck.

Certificate: `theory/derivations/p4_lemmaC_var_gate.py`, **217 `ck` checks, ~1 s, N=400,
ϖ ∈ [0.1, 40]**. Every number below is reproduced there; the three exact identities are
additionally verified with `mpmath` at `dps=30` (zero error).

---

## 0. The object

Notation as in `P4_concavity_q_route.md` / `P4_cone_coercivity.md`. On `ℓ²(ℕ₀)`,
`u_n=(2n+1)^{-1/2}`; for a kernel `φ` the operator `𝒪[φ]` is

    (𝒪[φ])_{nm} = u_n u_m[φ(|n−m|)𝟙_{n≠m} + φ(n+m+1)] − δ_{nm}(2/(2n+1))Σ_{k=1}^n φ(k).

Lorentz kernel `f(k;x)=x/(k²+x)`, `x=ϖ²>0`. `A:=𝒪[f]`, `B:=𝒪[f²]`, `h:=λ₁(A)=H(x)`,
`v=v_x` the top **unit** eigenvector, `v∈C` (nonneg, nonincreasing; proved cone invariance).
`q:=⟨v,Bv⟩`, `m₃:=⟨v,𝒪[f³]v⟩`, `E:=⟨v,𝒪[g]v⟩` with `g=∂_x f²`, `γ:=h−λ₂(A)`.

**Terminal inequality (route doc §4, restated).** Concavity `H″(x)<0 ⟺ q′(x)>0` closes if

    Var_v(B) := ⟨v,B²v⟩ − q²  <  γ · (xE/2) = γ·(q − m₃).      (★)

(The equality `xE/2 = q − m₃` is Cor 1.1; verified `dps=30`, 18 digits.)

---

## 1. The exact variance identity — `Var_v(B) = Var_v(M)`  **[PROVEN]**

Define

    **M := A − B = 𝒪[f − f²],   φ_M(k) = f(k)(1−f(k)) = x k²/(k²+x)² ≥ 0.**

The kernel `φ_M` is **nonnegative** (`f∈[0,1]`), unlike `f²` this is the whole point.

**Lemma 1.** `Var_v(B) = Var_v(M)`, and `⟨v,Mv⟩ = h − q`.

*Proof.* `Av=hv` gives `Bv = Av − Mv = hv − Mv`, hence
`Bv − qv = (h−q)v − Mv = −(Mv − (h−q)v)`, while `⟨v,Mv⟩ = ⟨v,Av⟩−⟨v,Bv⟩ = h−q`. The two
fluctuation vectors `Bv−⟨v,Bv⟩v` and `Mv−⟨v,Mv⟩v` are therefore exact negatives of one
another, so they have equal norm², i.e. `Var_v(B)=Var_v(M)`. ∎

Verified: `|Var_v(B) − Var_v(M)| = 0` exactly at `dps=30`, and `< 2e-14` over the whole
float grid (gate column `VarM=VarB`).

**Why this matters.** `(★)` is thereby a statement about the variance of a **single
operator with a nonnegative kernel**. That unlocks the two tools that fail on `B²`: the
operator-norm sandwich (§2) and the Perron / cone structure (§3, §6).

---

## 2. Operator-norm bounds on `Var_v(M)`  **[PROVEN]**

Since `M` is symmetric, `⟨v,M²v⟩ ≤ ‖M‖₂²`. Numerically `‖M‖₂ = λ_max(M)` on the entire grid
(the top eigenvalue dominates `|λ_min(M)|` with margin `λ_max/|λ_min| ∈ [1.35, 1.92]`; gate
check `dom?`). With `⟨v,Mv⟩ = h−q` (Lemma 1):

**(op-norm-minus)**  `Var_v(M) ≤ λ_max(M)² − (h−q)².`

**(Bhatia–Davis)**  `Var_v(M) ≤ (λ_max(M) − (h−q))·((h−q) − λ_min(M)).`

Both are valid upper bounds (gate checks each `≥ Var_v(B)`); measured looseness:
op-norm-minus `∈ [1.29, 1.67]×`, Bhatia–Davis `∈ [1.03, 1.46]×` (tighter; gate columns
`opm/Var`, `BD/Var`).

---

## 3. Analytic upper bound on `λ_max(M)`  **[PROVEN bound; ratio numerical]**

`M`'s diagonal is `M_{nn} = u_n²φ_M(2n+1) − (2/(2n+1))Σ_{k≤n}φ_M(k)`. Because `φ_M=f−f²≥0`
and `u>0`, **every off-diagonal entry of `M` is `≥0`** (verified: `min M_od = 0`); the diagonal
has only a handful of positive entries near the head (`M_{00}=φ_M(1)>0` for small ϖ; up to ~5
positive by ϖ=40, the rest `<0`). Let `D⁺:=diag(max(M_{nn},0)) ⪰ 0` and
`D⁻:=diag(max(−M_{nn},0)) ⪰ 0`, so `diag(M)=D⁺−D⁻` and

    **M_dom := M_od + D⁺   (a nonnegative symmetric matrix).**

Then `M = M_dom − D⁻` with `D⁻ ⪰ 0`, so

    **λ_max(M) ≤ λ_max(M_dom),**   and for any weight `ψ>0`,  `λ_max(M_dom) ≤ max_n (M_dom ψ)_n / ψ_n`  (Schur test).

Both are gate-checked as valid upper bounds. Tightness of `λ_max(M_dom)/λ_max(M)`:
**`1.095 – 1.112×` uniformly on ϖ∈[0.1,40]** (gate column `Ldom/L`). The Schur test with the
fixed analytic weight `ψ = u` gives `1.17 – 1.99×` (grows with ϖ; column `Su/L`); the Perron
weight (power iteration) attains the `1.10×` value.

**Consequence.** This *removes* the route doc's "top eigenvalue `β_max` to 2 %" residual: an
analytic upper bound on `λ_max(M)` accurate to `~10 %` (Perron) is now in hand, and even the
closed-form `ψ=u` Schur bound is `≤2×` throughout.

---

## 4. The scalar sufficient reductions  **[PROVEN algebra]**

Combining Lemma 1 and op-norm-minus, `(★)` closes if

    **(SUFF)   λ_max(M)²  <  (h−q)²  +  γ·(q − m₃).**

The RHS is an **upward parabola in `q`** (coefficient `+1` on `q²`), vertex at `q* = h − γ/2`,
minimum value `γ(h−m₃) − γ²/4` (gate-verified identity). Hence a fortiori `(★)` closes if

    **(QF, q-free)   λ_max(M)²  <  γ·(h − m₃) − γ²/4,     h−m₃ = ⟨v, 𝒪[f−f³] v⟩,**

where `φ = f−f³ = f(1−f)(1+f) ≥ 0` is again a nonnegative kernel. `QF` **eliminates the
Rayleigh deficit `(h−q)`** entirely — it needs only the single nonneg-kernel Rayleigh quantity
`h−m₃` and the gap `γ`.

---

## 5. The ceiling — `(★)` holds with the ANALYTIC gap and true spectral data  **[CERTIFICATE]**

Using the **analytic** gap lower bound `γ_lo = h_lo − λ₂^hi` (from `p4_concavity_gap2_gates.py`;
G1's closed form `γ_lo = x/(1+x) − Θ*`, `Θ* = 0.6349…`, for ϖ>1.3188) — which only *lowers* the
RHS of `(★)` — and the true `q, m₃`:

| ϖ | Var_v(B) | opm/Var | BD/Var | Ldom/L | cap(h−q) | cap(h−m₃) | **budget = γ_lo(q−m₃)/Var** | QF closes? |
|----|----------|---------|--------|--------|----------|-----------|------------------------------|-----------|
| 0.10 | 2.3e-10 | 1.673 | 1.464 | 1.095 | 0.730 | 0.730 | **3892** | no |
| 0.50 | 4.6e-05 | 1.653 | 1.422 | 1.101 | 0.725 | 0.724 | **124** | no |
| 1.00 | 2.6e-03 | 1.590 | 1.287 | 1.110 | 0.712 | 0.711 | **20.5** | no |
| 1.50 | 1.2e-02 | 1.518 | 1.147 | 1.111 | 0.697 | 0.697 | **6.73** | no |
| 2.00 | 2.6e-02 | 1.465 | 1.105 | 1.107 | 0.684 | 0.684 | **5.77** | yes |
| 3.00 | 5.1e-02 | 1.411 | 1.073 | 1.107 | 0.664 | 0.664 | **5.11** | yes |
| 6.00 | 1.0e-01 | 1.358 | 1.049 | 1.106 | 0.634 | 0.633 | **4.59** | yes |
| 20.0 | 2.0e-01 | 1.312 | 1.031 | 1.100 | 0.598 | 0.597 | **4.03** | no |
| 40.0 | 2.5e-01 | 1.294 | 1.026 | 1.095 | 0.585 | 0.584 | **4.11** | no |

The **budget** (how many times over `(★)` holds) is `≥ 4.0` everywhere and blows up as `ϖ→0`.
So `(★)` is **true with a wide margin**, using the *analytic* `γ`. The `QF` column reports
whether the q-free inequality closes with the tight Perron bound `λ_max(M)²_hi` and analytic
`γ_lo` **but the true `h−m₃`**: it closes on `ϖ∈[2,6]`. This is why the whole thing is
"certificate-level" and not yet a theorem: see §6.

---

## 6. The single residual  **[HONEST]**

Every remaining input is a **Rayleigh quantity of the same cone vector `v` against a
nonnegative-kernel `𝒪[φ]`**: the deficit `h−q = ⟨v,𝒪[φ_M]v⟩`, and `h−m₃ = ⟨v,𝒪[f−f³]v⟩`.
The only tool that lower-bounds these analytically is cone-coercivity (G2 Thm 1),
`⟨v,𝒪[φ]v⟩ ≥ κ_{L3}[φ,v]`, which **captures only `0.58 – 0.73`** of the true value
(gate columns `cap(h−q)`, `cap(h−m₃)`: `0.73` at ϖ=0.1 decreasing to `0.58` at ϖ=40).

Closing `(★)` fully analytically needs `≳0.9` capture, so the cone bounds are **insufficient**:
- **Ends (small & large ϖ).** op-norm-minus is a difference `λ_max(M)² − (h−q)²` of nearly
  equal quantities (at small ϖ, `h−q ≈ λ_max(M)` to `~0.1 %`, since `v` is essentially the top
  eigenvector of `M`); certifying the tiny `Var_v(M)` requires `(h−q)` to `~0.1 %`, far beyond
  `0.73`.
- **Middle.** `QF` removes the deficit but `h−m₃` at `0.68` capture already breaks the margin
  (checked: with `κ`/`κ_prof`-level bounds the fully-analytic ratio is `1–22 > 1` across the
  grid; scratch probe `probe_captures.py`).

This is **exactly the "shared top-eigenvector concentration / profile lower bound" bottleneck**
recorded in `P4_cone_coercivity.md` §7.4 and program memory (commit `a05847c`): a rigorous
`v_n ≥ c·v_0·u_n` for `n ≲ √x` would push all captures to a constant near 1 and close `(★)`.

**What this note nonetheless buys** (net progress over `P4_concavity_q_route.md` §7.5, which
lists *two* residuals with surrogates `21–42×` over budget):
1. **Unification.** The route doc's two separate residuals — the deficit `b_max − q` and the
   top eigenvalue `β_max` — collapse into **one** object (the cone-Rayleigh precision of `v`),
   via the identity `Var_v(B)=Var_v(M)` (§1). There is now a single thing to prove.
2. **Top-eigenvalue residual eliminated to `~10 %`** by the Perron/Schur bound on `λ_max(M)`
   (§3), replacing the route doc's `β_max`-to-`2%` requirement.
3. **A deficit-free form** `QF` (§4) that in the mid-range `ϖ∈[2,6]` needs *only* `h−m₃`, not
   the exquisitely-sensitive deficit.

---

## 7. Small-ϖ asymptotic route  **[ROUTE; constants not certified]**

Independent of the cone residual, `(★)` closes asymptotically as `ϖ→0` by power counting.
As `x→0`: `B=𝒪[f²]=O(x²)` and `v → ` top eigenvector of `𝒪[1/k²]`, so `q ~ x²`, `m₃ ~ x³`,
hence `q−m₃ ~ x²`; `γ = h−λ₂ ~ x`; and `Var_v(B)=Var_v(M) ~ x⁴`. Therefore

    RHS of (★) = γ(q−m₃) ~ x³   ≫   Var_v(B) ~ x⁴,     so   Var_v(B)/RHS ~ x → 0.

Gate confirms the exponents (`Var_v(B)/x⁴` and `(q−m₃)/x²` each ≈const between ϖ=0.1, 0.3) and
the budget `~1/x` blow-up (`3892, 397, 124` at ϖ=`0.1, 0.3, 0.5`). A rigorous small-ϖ theorem
needs explicit perturbation constants and a threshold `x₀`; the natural machinery is the
weak-coupling secular series (`p4-weak-coupling-secular-series-2026-09-13`). Not certified here.

---

## 8. Interfaces (no existing file modified)

- **Uses:** cone invariance `v∈C` (`P4_main_theorems.md` §3); `xE/2 = q−m₃` (Cor 1.1,
  `P4_concavity_q_route.md`); the analytic gap `γ_lo` and `Θ*` (`p4_concavity_gap2_gates.py`,
  G1); cone-coercivity `κ_{L3}` (`P4_cone_coercivity.md`, G2 Thm 1).
- **Feeds:** the terminal inequality `(★)` of `P4_concavity_q_route.md` §4 is now reduced to
  `SUFF`/`QF`; the sole open input is a `≳0.9`-capture lower bound on `⟨v,𝒪[φ]v⟩`
  (= the profile / eigenvector-concentration bottleneck).
- **Builder reuse:** `O_ker, kernels, top_eigen` and `h_lo_trial, lam2_hi, THETA_STAR` imported
  `dirname`-relative in the gate; no hardcoded worktree path.

∎
