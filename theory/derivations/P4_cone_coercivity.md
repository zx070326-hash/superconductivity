# P4 · Cone-coercivity lemma — a quantitative lower bound on E(x)

**状态 / status.** Delivered by agent `g2-cone-coercivity` on branch `g2-cone-coercivity`
(do **not** merge). The flagship analytic object of the concavity/q-monotonicity route
(`P4_concavity_q_route.md` §3, Theorem A) is here **upgraded from the qualitative
`E(x)>0` to a quantitative `E(x) ≥ κ(x) > 0`**, with three rigorous, independently
checkable bounds of increasing strength and an honest map of the residual obstruction.

Certificate: `theory/derivations/p4_cone_coercivity_probe.py`, **191 `ck` checks, ~5 s,
N=400/700**, ϖ ∈ [0.1, 40]. Every constant and every inequality below is reproduced there.

---

## 0. Notation and the object

On ℓ²(ℕ₀), `u_n := (2n+1)^{-1/2}` (strictly decreasing, `u_0=1`). For a kernel `φ:ℕ→ℝ`
the operator is (linear in φ; only the values `φ(k), k≥1`, ever enter)

    (𝒪[φ])_{nm} = u_n u_m [ φ(|n−m|)·𝟙_{n≠m} + φ(n+m+1) ]  −  δ_{nm}·(2/(2n+1))·Σ_{k=1}^{n} φ(k).

Lorentz kernel `f(k;x) = x/(k²+x)`, `x = ϖ² > 0` (convention `f(0)=1`, never used by 𝒪);
`A := 𝒪[f]`, `h(x) := λ₁(A) = H(x)`, `v = v_x` the top **unit** eigenvector.
By cone invariance (`P4_main_theorems.md` §3, proved) `v ∈ C := {w : w_n ≥ 0, w nonincreasing}`.
Derivative kernel `g(k) := ∂_x f² = 2x k²/(k²+x)³ ≥ 0`, `g(0)=0`, `g(1)=2x/(1+x)³ > 0`.
Target quantity `E(x) := ⟨v_x, 𝒪[g] v_x⟩`; by `P4_concavity_q_route.md` Cor 1.1, `xE/2 = q − m₃`.

Throughout write **`w := u∘v`, i.e. `w_n = u_n v_n`.** Since `u` and `v` are both nonnegative
and nonincreasing, **`w` is nonnegative and nonincreasing**, and `w_0 = u_0 v_0 = v_0`.
Note `‖v‖²=1` reads `Σ_n (2n+1) w_n² = 1`.

The subtlety (from `P4_concavity_q_route.md` §3, "重要的排除"): 𝒪[φ] is **not** PSD — its diagonal
carries the negative term `−(2/(2n+1))Σφ(k)`, and indeed `λ_min(𝒪[g]) < 0` numerically. So no
bound may go through `𝒪[φ] ⪰ 0`. All positivity below comes from the cone acting on the
**off-diagonal (L2/L3) structure**, made quantitative.

---

## 1. The master rung identity (pure algebra)

**Lemma 1.** For every `k ≥ 1` and every `v` with `w=u∘v ∈ ℓ²`,

    T_k := ⟨v, 𝒪[δ_k] v⟩ = 2·Σ_{n≥k} w_n (w_{n−k} − w_n)  +  Σ_{n=0}^{k−1} w_n w_{k−1−n}.   (1)

*Proof.* Insert `φ = δ_k` (so `φ(j)=𝟙_{j=k}`) into the definition:
- **L2** (`|n−m|=k`, `n≠m`): pairs `(n,n−k)` and `(n−k,n)` give `2 Σ_{n≥k} u_n u_{n−k} v_n v_{n−k} = 2 Σ_{n≥k} w_n w_{n−k}`.
- **L3** (`n+m+1=k`): `m=k−1−n≥0` gives `Σ_{n=0}^{k−1} u_n u_{k−1−n} v_n v_{k−1−n} = Σ_{n=0}^{k−1} w_n w_{k−1−n}`.
- **L1** diagonal: `Σ_{j=1}^{n} δ_k(j) = 𝟙_{k≤n}`, and `2/(2n+1)=2u_n²`, giving `−2 Σ_{n≥k} u_n² v_n² = −2 Σ_{n≥k} w_n²`.

Adding, `T_k = 2Σ_{n≥k} w_n w_{n−k} − 2Σ_{n≥k} w_n² + Σ_{n=0}^{k−1} w_n w_{k−1−n}`, and grouping the
first two sums as `2Σ_{n≥k} w_n(w_{n−k}−w_n)` gives (1). Absolute convergence of each sum follows
from `Σ|w_n w_{n−k}| ≤ ‖w‖²` (Cauchy–Schwarz). ∎

Verified to `≤ 4.4e-16` for `k ∈ {1,2,3,7,15}`, four values of x (probe block **(L1)**).

Since `𝒪[φ] = Σ_{k≥1} φ(k) 𝒪[δ_k]` (linearity; `φ(0)` never enters), Lemma 1 gives the exact
"rung expansion" for any kernel:

    ⟨v, 𝒪[φ] v⟩ = Σ_{k≥1} φ(k) T_k.                                                    (2)

---

## 2. The cone-coercivity theorem

**Theorem 1 (cone coercivity).** Let `v ∈ C` and let `φ:ℕ→[0,∞)` be any nonnegative kernel
(`φ(k) ≥ 0` for all `k≥1`). Then, with `w = u∘v`,

    ⟨v, 𝒪[φ] v⟩ ≥ κ_{L3}[φ,v] := Σ_{n,m≥0} φ(n+m+1) u_n u_m v_n v_m ≥ φ(1)·v_0².        (3)

*Proof.* Because `v ∈ C`, `w` is nonnegative and nonincreasing, hence for every `k≥1` the
telescoping sum in (1) satisfies `2Σ_{n≥k} w_n(w_{n−k}−w_n) ≥ 0` (each `w_{n−k} ≥ w_n ≥ 0`).
Therefore `T_k ≥ Σ_{n=0}^{k−1} w_n w_{k−1−n} ≥ 0`. Multiplying by `φ(k)≥0` and summing (2),

    ⟨v,𝒪[φ]v⟩ = Σ_k φ(k) T_k ≥ Σ_{k≥1} φ(k) Σ_{n=0}^{k−1} w_n w_{k−1−n}
              = Σ_{n,m≥0} φ(n+m+1) w_n w_m   (set k = n+m+1)
              = Σ_{n,m≥0} φ(n+m+1) u_n u_m v_n v_m = κ_{L3}[φ,v].

Finally, dropping every term but `(n,m)=(0,0)` (all terms are `≥0`) gives
`κ_{L3} ≥ φ(1) w_0² = φ(1) v_0²`. ∎

**Remarks.**
- Theorem 1 is exactly the **quantitative form of Corollary A′** (`P4_uniqueness_strategy.md` §2.A):
  A′ asserts `⟨v,𝒪[φ]v⟩ ≥ a_{L}²φ(1)Q_{LL}(δ₁) > 0` for *some* Abel coefficient; (3) sharpens the
  floor to `φ(1)v_0²`. The two are consistent: A′'s Lemma A gives `Q_{LL'}(δ₁) ≥ u_0² = 1` for **all**
  `L,L'` (its "group 3", the single term `n=0`), so with the Abel identity `Σ_L a_L = v_0` one has
  `Σ_{L,L'} a_L a_{L'} Q_{LL'}(δ₁) ≥ (Σ_L a_L)² = v_0²` — the same floor. Theorem 1's contribution is
  (i) the closed rung identity (1), and (ii) the far tighter intermediate `κ_{L3}` layer.
- It is **reusable and kernel-agnostic**: any `φ≥0`. Probe block **(T1)** verifies (3) for a *random*
  `v∈C` (not the eigenvector) and for `φ∈{g, f²}`, confirming the theorem does not use any spectral
  property of `v` beyond membership in the cone.
- `κ_{L3}` is genuinely the coercive content: numerically `κ_{L3}[g,v_x]/E(x) ∈ [0.582, 0.727]`
  **uniformly for x∈[0.01,1600]** — a constant-factor lower bound. The discarded telescoping part
  supplies the remaining ≤ 42 %.

---

## 3. Specialisation to E and the explicit floor

**Corollary 2 (E floor).** For every `x>0`,

    E(x) ≥ Σ_{n,m} g(n+m+1) u_n u_m (v_x)_n (v_x)_m ≥ g(1)·v_0(x)² = (2x/(1+x)³)·v_0(x)².   (4)

This is fully rigorous for all `x>0`. It is **tight for small/moderate x** (`g(1)v_0²/E → 0.70`
as `x→0`) but **loses a factor `x` at large x**: `g(1) = 2x/(1+x)³ ∼ 2/x²`, whereas `E ∼ 0.60/x`
(measured `E·x = 0.49, 0.55, 0.57, 0.59, 0.60, 0.61` at `x=4…6400`, probe/scaling), so
`g(1)v_0²/E ∼ (2v_0²/x²)/(0.6/x) ∼ 3v_0²/x → 0`. The floor keeps only the single rung `k=1`; the
mass of `g` sits near `k∼√x` (`g` peaks at `k=√x`, value `1/(4x)`), which (4) throws away.

To make (4) an explicit function of `x` we need a lower bound on the eigenvector head `v_0`.

---

## 4. Explicit lower bound on the eigenvector head v₀

**Proposition 3 (head bound via the test vector e₀).** With `f(1)=x/(1+x)`, `h=λ₁(A)`,
`λ₂ = λ₂(A)` the second eigenvalue,

    v_0(x)² ≥ (f(1) − λ₂) / (h − λ₂).                                                   (5)

*Proof.* Write `e_0 = v_0 v + r` with `r ⟂ v`, `‖r‖² = 1 − v_0²`. Then `⟨e_0,Ae_0⟩ = A_{00} = f(1)`
(the diagonal at `n=0`: L3 gives `u_0² f(1)=f(1)`, L1 is the empty sum). Expanding and using
`⟨v,Ar⟩ = h⟨v,r⟩ = 0`,

    f(1) = v_0² h + ⟨r,Ar⟩ ≤ v_0² h + λ₂‖r‖² = λ₂ + v_0²(h − λ₂),

since `r ∈ span{φ₂,φ₃,…}` where `A ⪯ λ₂`. Rearranging (and `h−λ₂ = gap > 0`) gives (5). ∎

The RHS of (5) is **decreasing in `λ₂`** (as `f(1)<h`), so replacing `λ₂` by any rigorous upper
bound `λ₂^hi < h` preserves the inequality:

    v_0² ≥ (f(1) − λ₂^hi)/(h − λ₂^hi),   whenever λ₂^hi < h.                            (5′)

This **reduces the head bound to the independently-established spectral gap** (no circularity: the
gap bounds of `2026-09-13_gap_lorentz_probe.md` and Lemma M do not use concavity). Two regimes:
`λ₂ ≤ 0.58 h` for `ϖ≤2`, `λ₂ ≤ Θ* = 0.6349` for `ϖ≥2`. Under these conservative bounds (5′) is
**positive and nonvacuous across the whole grid** (worst `0.062·v_0²` at `ϖ=2`, where the two
regimes meet; `0.20–0.50·v_0²` elsewhere; probe column `gapCr/v0^2`).

**The bound is intrinsically tight**, because the *true* `λ₂` is essentially zero: measured
`λ₂ = 4e-4` at `ϖ=1`, `1e-3` at `ϖ=2`, turning slightly negative for `ϖ≳20`. With the true `λ₂`,
(5) captures `0.64–0.90` of `v_0²` (probe column `v0lb/v0^2`), i.e. `v_0² ≳ f(1)/h`. Combining (4)
and (5′) yields the fully explicit, fully rigorous coercivity constant

    **κ_floor(x) = (2x/(1+x)³) · (f(1) − λ₂^hi)/(h − λ₂^hi).**                            (6)

---

## 5. Profile-improved κ with the correct order (up to a single log)

The factor-`x` loss of (4) is cured by keeping the whole `κ_{L3}` layer and lower-bounding the
eigenvector *profile*, not just its head.

**Proposition 4 (row-n profile bound).** For `n≥1`,

    v_n ≥ β_n·v_0,   β_n := u_n·(f(n)+f(n+1)) / (h − A_{nn}),                            (7)

where `A_{nn} = u_n²·(f(2n+1) − 2Σ_{j=1}^{n} f(j)) < h`; set `β_0 := 1`.

*Proof.* Row `n` of `Av=hv` reads `(h − A_{nn})v_n = Σ_{m≠n} A_{nm} v_m`. Every off-diagonal
`A_{nm} = u_n u_m(f(|n−m|)+f(n+m+1)) ≥ 0`, and every `v_m ≥ 0` (cone), so the RHS is `≥ A_{n0}v_0`.
Since `A_{n0} = u_n u_0(f(n)+f(n+1)) = u_n(f(n)+f(n+1))` and `h − A_{nn} = h − ⟨e_n,Ae_n⟩ > 0`
(as `h=λ₁ ≥ A_{nn}` with equality impossible for `n≥1`), divide to get (7). ∎

Feeding `w_n = u_n v_n ≥ (u_n β_n) v_0` into Theorem 1,

    **E(x) ≥ v_0²·P(x),   P(x) := Σ_{n,m≥0} g(n+m+1)·(u_n β_n)(u_m β_m).**               (8)

`P(x)` is explicit and finite (β decays like `u_n²`, β_n→0). Because (7) is decreasing in `h`,
replacing `h` by any rigorous upper bound `h^hi ≥ h` (e.g. the row-0 bound
`h ≤ f(1) + Σ_{m≥1}u_m(f(m)+f(m+1))`, or the enclosures in `P4_interchange_uniformity.md`) keeps
(8) a valid lower bound; combined with (5′) for `v_0²`, (8) is a fully explicit `κ_prof(x)`.

**Order.** For `n ≲ √x`, `f(n),f(n+1) ≈ 1` and `2u_n²Σ_{j≤n}f(j) ≈ 1`, so
`u_nβ_n ≈ 2u_n²/(h+1)` and `P ≈ 4Φ(x)/(h+1)²` with `Φ(x) := Σ_{n,m} g(n+m+1)u_n²u_m²`. Measured
(scaling probe): `Φ·x = 0.54, 0.77, 0.94, 1.12, 1.29, 1.46` at `x=4…6400`, i.e. `Φ ∼ 0.125·(ln x)/x`;
and `h ∼ 0.55·ln x`. Hence

    κ_prof = v_0²P ∼ v_0²·(0.5 ln x / x)/(0.31 ln²x) ∼ 1.6 v_0²/(x ln x),   while  E ∼ 0.60/x.

So `κ_prof` recovers `E` **up to a single factor `ln x`** — versus the factor-`x` loss of `κ_floor`.
Measured tightness `κ_prof/E`: `0.72 → 0.14` over `x∈[0.01,1600]` (probe column `kProf/E`), e.g.
`0.66` in the closure-critical window `ϖ∈[1.46,1.63]` (`x≈2.1–2.7`) and still `0.19` at `ϖ=20`,
against `0.49` and `0.004` for `κ_floor` there.

---

## 6. Certified tightness (probe, N=400)

| ϖ | x | E | κ_{L3}/E | κ_floor(true v₀)/E | v₀-bound/v₀² | κ_prof/E |
|---|---|---|---|---|---|---|
| 0.20 | 0.04 | 0.08501 | 0.726 | 0.702 | 0.901 | 0.724 |
| 1.00 | 1.00 | 0.32287 | 0.710 | 0.616 | 0.878 | 0.698 |
| 1.54 | 2.37 | 0.19011 | 0.695 | 0.490 | 0.853 | 0.661 |
| 2.00 | 4.00 | 0.12195 | 0.684 | 0.377 | 0.832 | 0.621 |
| 5.00 | 25.0 | 0.02180 | 0.640 | 0.076 | 0.745 | 0.404 |
| 20.0 | 400  | 0.00147 | 0.596 | 0.004 | 0.655 | 0.187 |
| 40.0 | 1600 | 0.000375| 0.582 | 0.001 | 0.636 | 0.138 |

`κ_{L3}` (Theorem 1, requires the profile `w`) is the constant-factor content (`≥0.58 E`);
`κ_prof` (explicit given `h^hi`) preserves the order up to `ln x`; `κ_floor` (closed form) is exact
for small `x`. N-consistency (N=400 vs 700) of `κ_{L3}` and `E`: relative drift `≤ 5.5e-7`.

---

## 7. Honest range and the residual obstruction

1. **Theorem 1 and Corollary 2** hold for **all `x>0`**, unconditionally. The floor `g(1)v_0²` is
   the honest "closed-form" `κ`; it is only order-correct for `x = O(1)`.

2. **Proposition 3 / (5′)** is rigorous **given a spectral-gap upper bound `λ₂^hi<f(1)`**. Such
   bounds exist independently (`0.58h`, `Θ*`); with them `κ_floor` (6) is fully explicit and
   nonvacuous on `[0.1,40]`. Its looseness (worst `0.062 v_0²` at `ϖ=2`) reflects the crudeness of
   the *imported* gap bound, not of (5): the true `λ₂≈0` makes (5) tight (`≥0.64 v_0²`).

3. **Proposition 4 / (8)** is rigorous **given `h^hi ≥ h`**; `κ_prof` loses only `∼ ln x` versus `E`.

4. **The open gap — a lower bound on the eigenvector profile.** Closing the residual `ln x`
   (turning `κ_prof ∼ 1/(x ln x)` into the true `∼ 1/x`) requires improving (7).
   A rigorous statement of the form **`v_n ≥ c·v_0·u_n` for `n ≲ √x`** (equivalently
   `w_n ≥ c v_0 u_n²`) would immediately give `κ_prof ≳ c²v_0²Φ(x) ∼ 1/x`, matching `E` to a
   constant. Its two sources of loss are precisely locatable:
   - **(a) keeping only the `m=0` rung.** The row-`n` equation has mass at `m≈n` (the L2 band), so
     `A_{n0}v_0` underestimates `Σ_{m≠n}A_{nm}v_m` by a per-component factor
     `v_n/(β_n v_0) ∈ [1.1, 2.8]` (measured, growing slowly with ϖ);
   - **(b) the `h` in the denominator.** In the regime that matters (`n≲√x`),
     `u_nβ_n ≈ 2u_n²/(h+1)` and `P ≈ 4Φ(x)/(h+1)²`; with `h ∼ 0.55 ln x` this injects the residual
     `ln x`. The numerically-observed profile is instead `v_n ≈ v_0·u_n` (i.e. `w_n ≈ v_0 u_n²`) for
     `n ≲ √x` — with **no** `h`-suppression — reflecting `A(x) → 2uuᵀ − I` as `x→∞`.

   This is exactly the **"shared top-eigenvector concentration" bottleneck** flagged in
   the program memory (commit `a05847c`, shared with `smooth-coulomb`): the top vector of `A`
   concentrates on the `u`-profile up to the crossover `n∼√x`, but a *lower* bound on that
   concentration is not yet proved. We record it as the single clean open problem.

5. **What is NOT claimed.** No closed-form `κ(x)` with the exact `1/x` large-`x` order is proved;
   only `κ_prof ∼ 1/(x ln x)` (order-correct up to a log) is. The constant-factor result
   `κ_{L3} ≥ 0.58 E` is **numerical** (Theorem 1 proves `κ_{L3} ≤ E` and `κ_{L3} ≥ g(1)v_0²`; the
   `0.58` ratio is measured, not proved).

---

## 8. Interfaces (no existing file modified)

- **Uses:** cone invariance `v∈C` (`P4_main_theorems.md` §3); the Abel/rung structure and Cor A′
  (`P4_uniqueness_strategy.md` §2.A) — Theorem 1 is its quantitative upgrade; the spectral-gap
  upper bounds `0.58h` / `Θ*` (`2026-09-13_gap_lorentz_probe.md`, Lemma M) for (5′)/(6).
- **Feeds:** the terminal inequality `Var_v(B) < γ·xE/2` of `P4_concavity_q_route.md` §4 now has an
  explicit lower bound on `xE/2 = q − m₃`, namely `xE/2 ≥ (x/2)·κ` for any of `κ∈{κ_floor, κ_prof}`
  (and `≥ (x/2)κ_{L3}` given the profile). For the closure-critical window `ϖ∈[1.46,1.63]`,
  `κ_prof ≥ 0.66 E`.
- **Builder reuse:** `O_ker, kernels, top_eigen` imported from `p4_concavity_q_route.py` via a
  `dirname`-relative path (probe header); no hardcoded worktree path.

∎
