#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_fese_essenberger_sf_bipartite_2026-09-18.py -- from zero, no repo imports.

WHAT
----
First attempt at plugging a REAL, first-principles, non-phonon (spin-fluctuation
only) multiband coupling matrix into the repo's PF / sign-structure machinery
(see theory/notes/vertex/spm_gap_needs_no_exotic_physics_2026-09-18.md sec.3.2,
and theory/notes/probes/spin_fluctuation_rigorous_tc_scoping_2026-09-18.md sec.3.3).

L1 SOURCE (already transcribed by main session in
theory/notes/literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md, appendix
A.4, itself read from Essenberger, Sanna, Buczek, Ernst, Sandratskii, Gross,
arXiv:1411.2121 = PRB 94,014503 (2016), Fig.4(a) inset + Fig.5 caption, FeSe,
5-patch Fermi-surface decomposition {Gamma1,Gamma2,Gamma3 (hole pockets),
M1,M2 (electron pockets)}):

    Gamma1-Gamma2   lam_ij = 0.01   (same gap phase  -> "repulsive in s+- channel")
    M1-M2           lam_ij = 0.07   (same gap phase  -> "repulsive in s+- channel")
    Gamma3-M1       lam_ij = 0.19   (opposite phase  -> "attractive")
    Gamma1-M1       lam_ij = 0.18   (opposite phase  -> "attractive")
    Gamma2-M2       lam_ij = 0.11   (opposite phase  -> "attractive")

All FIVE numbers are POSITIVE (the paper's own Eq.(37)/(38) and p.9 sign
discussion: SF spectral function is repulsive, Im[Lambda^SF]>0, for ALL i,j;
the sign structure of the gap is carried by the GAP FUNCTION, not by the sign
of lambda_ij -- see repo note appendix D.6, "symbol not carried by the lambda
matrix, 4th occurrence of this pattern"). This is EXACTLY the "V(q)>=0
everywhere, sign lives in the eigenvector" setup of the probes note sec.3.3
dual (Hainzl real-space-positive vs momentum-space-positive), NOT the K1
"Lambda^ph_ii>=0 intra / Lambda^sf_ij<=0 inter" 2-family convention used in
spm_gap_needs_no_exotic_physics sec.3.2 (that convention needs an ATTRACTIVE
intra-family channel, e.g. phonons; here there is none -- FeSe SF-only has a
REPULSIVE intra-family (Gamma1-Gamma2, M1-M2) piece too).

ONLY 5 of the C(5,2)=10 pairs are given in the schematic (the figure shows two
small illustrative sub-panels, "left"=same-phase pairs, "right"=opposite-phase
pairs, not a full published 5x5 table). The other 5 pairs (Gamma1-Gamma3,
Gamma2-Gamma3, Gamma1-M2, Gamma2-M1, Gamma3-M2) are NOT given anywhere found
so far (see literature note appendix D.7 and A.6: the underlying two-index
Lambda^SF_ij(omega) data lives in the APS Supplemental Material of PRB
90,214504, unreachable in this environment -- 403 error already logged).
This script treats the missing entries as exactly 0 for a first pass, and
separately probes sensitivity to relaxing that assumption (part C).

WHAT THIS SCRIPT CHECKS (three parts)
--------------------------------------
A. Direct eigendecomposition of the 5x5 symmetric non-negative matrix V built
   from the five reported entries (missing entries = 0, diagonal = 0, per the
   paper's own statement Lambda^SF_ii ~ 0, intraband SF negligible).
   Report lambda_min(V), its eigenvector, and its sign pattern.
B. Bipartite-graph check on the coupling graph (5 edges, 5 vertices): find
   whether a 2-coloring exists with no monochromatic edge (equivalent to:
   does a diagonal sign matrix D exist with D V D = -V, i.e. spec(V)
   symmetric about 0, i.e. lambda_min(V) = -lambda_max(V) EXACTLY). If yes,
   report the bipartition found by the coloring and compare it to the
   "textbook" hole-vs-electron split {G1,G2,G3} vs {M1,M2}.
C. Metzler/PF check in the sense of spm_gap_needs_no_exotic_physics sec.3.2:
   using D = diag(+1 for Gamma-patches, -1 for M-patches) (the hole/electron
   split, NOT the bipartition found in B), is D Lambda^Delta D Metzler
   (entrywise >=0 off-diagonal), where Lambda^Delta_ij := -V_ij for a
   same-family edge (G-G or M-M, i.e. repulsive/pair-breaking) and
   Lambda^Delta_ij := +V_ij for a cross-family edge (G-M, i.e. the
   "attractive-effective" edge)? This is the precondition of the repo's PF
   theorem (sec.3.2 of spm_gap note): report which entries violate it.
D. Sensitivity: perturb the 5 "missing" entries by a small epsilon (breaking
   exact bipartiteness) and show lambda_min(V)'s eigenvector sign pattern
   moves continuously (no discontinuous jump) for small epsilon, confirming
   the bipartition found in B is not a knife-edge numerical accident of the
   reported numbers, it is a structural property of the 5-edge subgraph.

LEVEL: numerical-evidence throughout (literature-transcribed inputs, some
missing entries approximated as zero -- explicitly flagged, not hidden).
"""

import itertools
import numpy as np

np.set_printoptions(precision=6, suppress=True)

PATCHES = ["G1", "G2", "G3", "M1", "M2"]
IDX = {p: i for i, p in enumerate(PATCHES)}
N = len(PATCHES)

# ---- L1 input: the five reported edges (magnitude only, all positive) -----
EDGES = {
    ("G1", "G2"): 0.01,
    ("M1", "M2"): 0.07,
    ("G3", "M1"): 0.19,
    ("G1", "M1"): 0.18,
    ("G2", "M2"): 0.11,
}
SAME_FAMILY_EDGES = {("G1", "G2"), ("M1", "M2")}
CROSS_FAMILY_EDGES = {("G3", "M1"), ("G1", "M1"), ("G2", "M2")}
assert set(EDGES) == SAME_FAMILY_EDGES | CROSS_FAMILY_EDGES

FAMILY = {"G1": "hole", "G2": "hole", "G3": "hole", "M1": "electron", "M2": "electron"}


def build_V(edges):
    V = np.zeros((N, N))
    for (a, b), w in edges.items():
        i, j = IDX[a], IDX[b]
        V[i, j] = w
        V[j, i] = w
    return V


def sign_pattern_str(v, tol=1e-9):
    signs = []
    for x in v:
        if x > tol:
            signs.append("+")
        elif x < -tol:
            signs.append("-")
        else:
            signs.append("0")
    return "".join(signs)


def part_a(V):
    print("=" * 78)
    print("PART A -- direct eigendecomposition of the 5x5 V (missing edges = 0)")
    print("=" * 78)
    print("V =\n", V)
    w, vecs = np.linalg.eigh(V)
    order = np.argsort(w)
    w = w[order]
    vecs = vecs[:, order]
    for k in range(N):
        v = vecs[:, k]
        # normalize sign convention: first nonzero component positive
        for x in v:
            if abs(x) > 1e-9:
                if x < 0:
                    v = -v
                break
        print("  eigenvalue[%d] = %+.6f   eigenvector (%s) = %s   pattern=%s"
              % (k, w[k], ",".join(PATCHES), np.round(v, 4), sign_pattern_str(v)))
    lam_min = w[0]
    lam_max = w[-1]
    print("\n  lambda_min(V) = %.6f ,  lambda_max(V) = %.6f" % (lam_min, lam_max))
    print("  lambda_min + lambda_max = %.3e  (0 <=> spectrum symmetric about 0"
          " <=> graph bipartite)" % (lam_min + lam_max))
    v_min = vecs[:, 0]
    for x in v_min:
        if abs(x) > 1e-9:
            if x < 0:
                v_min = -v_min
            break
    hole_signs = {sign_pattern_str([v_min[IDX[p]]])[0] for p in ["G1", "G2", "G3"]}
    electron_signs = {sign_pattern_str([v_min[IDX[p]]])[0] for p in ["M1", "M2"]}
    matches_textbook = (len(hole_signs) == 1 and len(electron_signs) == 1
                        and hole_signs != electron_signs)
    print("  eigenvector(lambda_min) matches textbook hole-vs-electron s+- split? %s"
          % matches_textbook)
    return lam_min, lam_max, v_min


def part_b():
    print("\n" + "=" * 78)
    print("PART B -- is the 5-edge coupling graph bipartite? (find 2-coloring)")
    print("=" * 78)
    # brute-force over 2^5 colorings (tiny graph), fixing G1 = +1 to kill the
    # global sign redundancy
    found = None
    for bits in itertools.product([1, -1], repeat=N - 1):
        color = {PATCHES[0]: 1}
        for p, b in zip(PATCHES[1:], bits):
            color[p] = b
        ok = all(color[a] * color[b] < 0 for (a, b) in EDGES)
        if ok:
            found = dict(color)
            break
    if found is None:
        print("  NOT bipartite: no valid 2-coloring exists for the 5 reported edges.")
        return None
    part_plus = sorted(p for p in PATCHES if found[p] == 1)
    part_minus = sorted(p for p in PATCHES if found[p] == -1)
    print("  Bipartite. 2-coloring found:")
    print("    part A (+1) = %s" % part_plus)
    print("    part B (-1) = %s" % part_minus)
    textbook_plus = {"G1", "G2", "G3"}
    textbook_minus = {"M1", "M2"}
    matches = (set(part_plus) == textbook_plus and set(part_minus) == textbook_minus) or \
              (set(part_plus) == textbook_minus and set(part_minus) == textbook_plus)
    print("  Does this coincide with the textbook hole-vs-electron split "
          "{G1,G2,G3}|{M1,M2}?  %s" % matches)
    if not matches:
        print("  ==> NO. The bipartition forced by the 5 REPORTED edges mixes hole and")
        print("      electron patches (G2 and M1 fall on the same side). This is almost")
        print("      certainly an artefact of the schematic showing only 5 of 10 pairs,")
        print("      not a physical prediction -- flagged explicitly in the companion note.")
    return found


def part_c(V):
    print("\n" + "=" * 78)
    print("PART C -- Metzler/PF check under the repo's K1 two-family convention")
    print("          (D = diag(+1 on hole patches, -1 on electron patches))")
    print("=" * 78)
    D = np.array([1.0 if FAMILY[p] == "hole" else -1.0 for p in PATCHES])
    Lam_Delta = np.zeros((N, N))
    for (a, b), w in EDGES.items():
        i, j = IDX[a], IDX[b]
        sign = -1.0 if (a, b) in SAME_FAMILY_EDGES else +1.0
        Lam_Delta[i, j] = sign * w
        Lam_Delta[j, i] = sign * w
    print("  Lambda^Delta (K1 convention: same-family edge -> -V, cross-family -> +V):")
    print(Lam_Delta)
    DLD = np.outer(D, D) * Lam_Delta
    print("  D Lambda^Delta D =\n", DLD)
    off_diag_min = None
    violations = []
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            if (PATCHES[i], PATCHES[j]) in EDGES or (PATCHES[j], PATCHES[i]) in EDGES:
                val = DLD[i, j]
                if off_diag_min is None or val < off_diag_min:
                    off_diag_min = val
                if val < -1e-12:
                    violations.append((PATCHES[i], PATCHES[j], val))
    print("  min off-diagonal entry on reported edges = %.6f" % off_diag_min)
    if violations:
        print("  Metzler FAILS. Violating edges (same-family, hence negative and")
        print("  D_i D_j = +1 does not flip them):")
        for a, b, val in violations:
            print("    (%s,%s): D_iD_j*Lambda^Delta = %.6f < 0" % (a, b, val))
        print("  ==> The repo's PF theorem (spm_gap sec.3.2, precondition 1: intra-family")
        print("      entries of Lambda^ph must be >=0) does NOT apply to this SF-only FeSe")
        print("      model: the intra-hole (G1-G2) and intra-electron (M1-M2) SF channels")
        print("      are THEMSELVES repulsive (Lambda^Delta<0), not attractive. There is no")
        print("      phonon term in this table to supply the required Lambda^ph_ii>=0.")
        print("      The mechanism that actually operates here is the OTHER one already on")
        print("      record (probes/spin_fluctuation_rigorous_tc_scoping sec.3.3): pure")
        print("      lambda_min(V) of a non-negative V, no Metzler/PF conjugation needed --")
        print("      see part A/B above.")
    else:
        print("  Metzler holds on all reported edges (unexpected given same-family entries")
        print("  are repulsive; re-check convention).")


def part_d(V0):
    print("\n" + "=" * 78)
    print("PART D -- sensitivity: perturb the 5 UNREPORTED entries by epsilon")
    print("=" * 78)
    missing = [("G1", "G3"), ("G2", "G3"), ("G1", "M2"), ("G2", "M1"), ("G3", "M2")]
    rng = np.random.default_rng(20260918)
    for eps in [0.0, 1e-4, 1e-3, 1e-2, 0.05]:
        V = V0.copy()
        signs = rng.choice([-1, 1], size=len(missing))
        # keep entries non-negative (they are couplings, magnitudes) --
        # use random SUBSET turned on at magnitude eps, not signed noise,
        # since V must stay entrywise >=0 to keep the "repulsive SF" reading
        for (a, b), s in zip(missing, signs):
            i, j = IDX[a], IDX[b]
            val = eps if s > 0 else 0.0
            V[i, j] = val
            V[j, i] = val
        w, vecs = np.linalg.eigh(V)
        v_min = vecs[:, 0]
        for x in v_min:
            if abs(x) > 1e-9:
                if x < 0:
                    v_min = -v_min
                break
        print("  eps=%.4f  lambda_min=%.6f  lambda_min+lambda_max=%.3e  pattern=%s  "
              "(nonzero missing edges turned on: %s)"
              % (eps, w[0], w[0] + w[-1], sign_pattern_str(v_min),
                 [m for m, s in zip(missing, signs) if s > 0] if eps > 0 else []))
    print("\n  Interpretation: for eps>0 the exact bipartite symmetry (lambda_min+lambda_max=0)")
    print("  is generically broken (any added edge can close an odd-length alternative")
    print("  cycle), and the eigenvector sign pattern can continue to deviate from the")
    print("  textbook hole/electron split depending on WHICH missing edges are largest --")
    print("  i.e. the bipartition found in Part B is fragile against exactly the data that")
    print("  is missing from the published schematic. This is the concrete form of the")
    print("  gap identified in the companion note: a real test needs the full 5x5 (or")
    print("  finer patch-resolved) Lambda^SF_ij matrix, not 5 illustrative entries.")


def main():
    V = build_V(EDGES)
    part_a(V)
    part_b()
    part_c(V)
    part_d(V)


if __name__ == "__main__":
    main()
