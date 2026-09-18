# Two-Parameter Eliashberg Operator and Boson-Mediated Tc Upper Bound

_Date: 2026-09-16_
_Status: numerical-evidence (float64 throughout, no interval arithmetic)_

## 1. Setup: Two-Channel Eliashberg with Unequal Couplings

### 1.1 Physical motivation

In phonon-mediated superconductivity, the Eliashberg gap equation couples the
gap function Delta and the mass renormalization Z through the same spectral
function alpha^2 F(omega).  In non-phonon mechanisms (spin fluctuations, excitons),
the Z channel and the Delta channel generally sample different regions of
momentum space and use different coupling constants:

- **Delta channel** (pairing): alpha^2 F_Delta(omega) with lambda_Delta
- **Z channel** (mass renormalization): alpha^2 F_Z(omega) with lambda_Z

The ratio r = lambda_Z / lambda_Delta >= 1 captures the channel asymmetry:
- Phonon: r = 1 (exact, same vertex)
- Iron-based (s+-): r ~ 1.5--2
- Cuprate paramagnon (d-wave): r ~ 2--3

### 1.2 The two-parameter linear operator

The KAY linearized gap equation at mu* = 0 with Lorentzian kernel
f_w(k) = w^2/(k^2 + w^2) generalizes to:

    H_r(w) = -r * H1(w) + H2(w) + H3(w)

where H1 is the diagonal (Z channel) and H2 + H3 is the off-diagonal (pairing
channel), precisely:

    (H1)_{nn} = (1/(2n+1)) sum_{k=1}^n 2 f(k)            [diagonal, Z]
    (H2)_{nm} = (1-delta_{nm}) f(|n-m|) / sqrt((2n+1)(2m+1))   [pairing]
    (H3)_{nm} = f(n+m+1) / sqrt((2n+1)(2m+1))                  [pairing]

The Tc equation is  lambda_Delta * lambda_1(H_r(w/T)) = 1, and the operator
algebra underlying the Tc bound works IDENTICALLY to the single-parameter case
(the operator is still compact, self-adjoint, the dependence on the spectral
measure is still linear in P).

### 1.3 The P_Z = r * P_Delta simplification

When the Z and Delta spectral functions have the same shape but different
amplitudes (the "simplest" two-parameter case, which already covers spin
fluctuations in the single-boson approximation), the linearized operator at
fixed w is exactly H_r(w).  The full gap equation becomes:

    lambda_Delta * lambda_1(int H_r(w) P_Delta(dw) / T) = 1

## 2. Key Identities (Verified Numerically)

### 2.1 Identity (I1_r): Second x-derivative

With x = w^2 and A_r(x) = H_r(sqrt(x)),

    d^2/dx^2 A_r(x) = -2 O_r[phi_x]

where phi_x(k) = k^2/(k^2 + x)^3 and O_r[f] = -r L1[f] + L2[f] + L3[f].

**Status**: verified to relative error < 1e-5 by finite differences at N=60
for r in {1, 2, 5} and w in {0.3, 1, 5, 15}.  The proof is purely algebraic
(same elementwise cancellation as in KAY, r factors out of L1 only):

    f''(w) - f'(w)/w = 2k^2(k^2-3w^2)/(k^2+w^2)^3 - 2k^2w/(k^2+w^2)^2/w
                     = -8 w^2 k^2/(k^2+w^2)^3

This is r-INDEPENDENT elementwise; the r enters only through the L1 sum,
and O_r[phi] = -r L1[phi] + L2[phi] + L3[phi] is exactly (I1_r).

### 2.2 Identity (I2_r): First x-derivative

    (I - P_r) A_r'(x) v_r = -x (I - P_r) O_r[psi_x] v_r

where psi_x(k) = 1/(k^2 + x)^2 and P_r = |v_r><v_r|.

**Status**: verified to relative error < 1e-4 at N=80 for r in {1, 2, 5}.

### 2.3 Concavity formula (C**)

Second-order perturbation theory gives:

    H_r''(x) = -2 F_r(w) + positive_resolvent_term

where
- F_r(w) = <v_r, O_r[phi_x] v_r>  (the "concavity driver")
- positive_resolvent_term = 2 x^2 sum_{j>=2} |<v_j, O_r[psi] v_r>|^2 / (lambda_1 - lambda_j)
- Bounded above by 2 x^2 E_r^2 / g_r, with
  E_r = ||(I-P_r) O_r[psi] v_r||, g_r = spectral gap

## 3. Numerical Results

### 3.1 The limiting constant g_r(2) = lambda_1(G_r(2))

| r    | g_r(2)       | C_inf(r)     | C/C(1) | Physical system             |
|------|-------------|-------------|--------|------------------------------|
| 1.0  | 1.31814021  | 0.18272625  | 1.0000 | phonon (exact KAY)           |
| 1.2  | 1.29904704  | 0.18139803  | 0.9927 |                              |
| 1.5  | 1.27387346  | 0.17963182  | 0.9831 | iron-based s+- (lower)       |
| 2.0  | 1.23939466  | 0.17718418  | 0.9697 | iron-based s+- (upper)       |
| 2.5  | 1.21204481  | 0.17521830  | 0.9589 | cuprate paramagnon (lower)   |
| 3.0  | 1.18996103  | 0.17361470  | 0.9501 | cuprate paramagnon (upper)   |
| 5.0  | 1.13309657  | 0.16941568  | 0.9272 | extreme asymmetry            |
| 10.0 | 1.07527680  | 0.16503660  | 0.9032 | extreme asymmetry            |

Key observations:
1. g_r(2) is **monotonically decreasing** in r (verified up to r=30)
2. g_r(2) > 0 for all tested r up to 30 -- the bound **exists** for all physical r
3. The spectral gap remains O(1) for all r -- simple-eigenvalue condition holds

### 3.2 Concavity of H_r(x)

**Result**: H_r''(x) < 0 for ALL tested r in {1, 1.2, 1.5, 2, 2.5, 3, 4, 5}
and ALL x in (0, 2000] (92-point grid, N=400).

The "worst" (least negative) H_r''(x) occurs at large x for all r:
- r=1: worst H''(x) = -1.23e-7 at x=2000
- r=5: worst H''(x) = -7.38e-8 at x=2000

The concavity is ROBUST: the ratio (positive_term) / (2*F_r) stays well below 1
everywhere:

| r   | w    | F_r       | pos_term  | ratio pos/(2F) |
|-----|------|-----------|-----------|-----------------|
| 1.0 | 1.0  | 1.61e-1   | 6.07e-3   | 0.0188          |
| 1.0 | 10.0 | 2.85e-5   | 8.85e-6   | 0.1552          |
| 3.0 | 1.0  | 1.47e-1   | 4.80e-3   | 0.0164          |
| 3.0 | 10.0 | 2.12e-5   | 7.31e-6   | 0.1729          |
| 5.0 | 1.0  | 1.40e-1   | 3.55e-3   | 0.0127          |
| 5.0 | 10.0 | 1.56e-5   | 5.65e-6   | 0.1809          |

The worst ratio is ~0.19 at r=5, w=20 -- still far from 1.

### 3.3 Spectral gap stability

| r   | w=0.3 | w=1   | w=5   | w=20  |
|-----|-------|-------|-------|-------|
| 1.0 | 0.110 | 0.715 | 2.203 | 3.553 |
| 2.0 | 0.104 | 0.668 | 1.993 | 3.276 |
| 5.0 | 0.096 | 0.607 | 1.650 | 2.753 |

The gap decreases with r but remains bounded away from 0.

## 4. The Tc Upper Bound (Theorem Conjecture)

### 4.1 Statement (numerical-evidence)

**Conjecture (Two-Parameter Tc Bound).**
For any boson-mediated pairing mechanism with:
- Pairing spectral function P_Delta >= 0, coupling lambda_Delta
- Mass-renormalization spectral function P_Z = r * P_Delta, r >= 1
- mu* = 0 (no Coulomb pseudopotential)

the critical temperature satisfies:

    **Tc <= C_inf(r) * sqrt(lambda_Delta * <Omega^2>_Delta)**

where C_inf(r) = sqrt(g_r(2)) / (2 pi) and g_r(2) = lambda_1(G_r(2)).

### 4.2 Proof structure (if upgraded from numerical-evidence)

1. **lam_1 convexity** (operator level):
   lam_1(int H_r(w) P(dw)) <= int lam_1(H_r(w)) P(dw) = int H_r(w^2) P(dw)
   [This step uses that lam_1 is a convex function of the symmetric operator.
   Proven in the single-parameter case by the same argument.]

2. **Jensen reduction** (concavity of H_r):
   int H_r(x) P_x(dx) <= H_r(E_P[x]) = H_r(<Omega^2>)
   [This is where H_r''(x) < 0 enters -- certified cell-by-cell in the
   single-parameter case, here numerical-evidence for all tested r.]

3. **Small-x asymptotics** (Einstein extremal):
   H_r(x) ~ g_r(2) * x as x -> 0, giving 1/lambda_Delta = g_r(2) * Omega^2/T_c^2,
   hence T_c = sqrt(g_r(2) * lambda_Delta * Omega^2).

### 4.3 Physical implications

- **Phonon bound reproduced**: C_inf(1) = 0.18273, consistent with KAY.

- **Universal suppression**: The channel asymmetry r > 1 SUPPRESSES Tc:
  C_inf(r) / C_inf(1) = sqrt(g_r(2) / g_1(2)) < 1 for r > 1.
  This is physically correct: stronger mass renormalization without
  proportionally stronger pairing reduces Tc.

- **Quantitative bounds**:
  - Cuprate paramagnon (r ~ 3): Tc bound is 95.0% of the phonon bound
  - Iron-based s+- (r ~ 1.5): Tc bound is 98.3% of the phonon bound
  - The suppression is SURPRISINGLY MILD -- only 5% for cuprates.
    This means that channel asymmetry is NOT the main limiter of
    spin-fluctuation Tc; rather, it is lambda_Delta itself (which is
    typically smaller than the phonon lambda).

- **The bound never vanishes**: g_r(2) > 0 for all r tested up to 30,
  so the Tc bound EXISTS for arbitrarily asymmetric channels.

## 5. What Remains for Rigorous Status

To upgrade from numerical-evidence to certificate:
1. The concavity proof machinery of `p4_extremal_concavity.py` (Stage S + L)
   generalizes straightforwardly: replace O[f] by O_r[f] throughout.
   The r only enters through the L1 diagonal, so all envelope estimates
   get a factor of max(1, r) at most.
2. The spectral gap certificate needs G_r(2) instead of G(2) -- the
   same Cholesky-deflation approach works.
3. The small-w bridge H_r(w)/w^2 -> G_r(2) is verified.

The main RISK is that the gap g_r(w) shrinks with r, making the concavity
certificate harder (more cells needed, or larger trial vectors).  But the
gap stays O(1), so this is a computational cost issue, not a structural one.

## 6. Scripts

- Concavity probe: `theory/notes/_scratch_ir/probe_two_param_concavity.py`
- Decomposition + synthesis: `theory/notes/_scratch_ir/probe_two_param_part2.py`

## 7. Literature Context

As identified by the explore-universal-bound investigation: there is NO prior
work on rigorous Tc upper bounds for spin-fluctuation or exciton-mediated
superconductivity.  The entire literature on Tc bounds (Esterlis-Kivelson 2018,
Allen-Dynes 1975, KAY 2024) assumes phonon mediation.  This two-parameter
extension, even at numerical-evidence level, is therefore the FIRST result
of its kind.

> **⚠ 订正（2026-09-16 分诊 + 主会话）**：上面这句优先权主张（"FIRST result of its kind"）**越级**，已撤下。
> 它完全建立在一个子 agent 的文献**否定性**结论（"there is NO prior work"）之上，而该否定性结论未经一手核实
> ——本仓库已有 agent 编造引用与作者归属的前例（见 `goldilocks_literature_survey_2026-09-16.md` 顶部订正）。
> 「查不到」≠「不存在」；优先权主张需要**穷尽性**证据，不是一次检索的空结果。
>
> **应读作**：本轮检索未发现自旋涨落/激子介导 Tc 严格上界的先例；若要写进论文，须先做一手文献穷尽核查。
> 本节其余数值分级正确，不受影响。
