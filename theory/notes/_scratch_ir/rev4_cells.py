"""rev4: (1) perturbation lemma (**) tested against float64 truth; (2) max Psi-bar recomputed
from the certified cell numbers printed by p4_conj2_small_w.py; (3) cell-uniform Lam_hi vs Lam at
the LEFT end (where Lam is largest); (4) Lemma A sanity on a 2-atom P."""
import numpy as np
from rev4_truth import build_fast, fG, fH, fK
N = 1500
G2 = build_fast(N, fG(2)); ev, V = np.linalg.eigh(G2); g2 = ev[-1]; v0 = V[:, -1]; v0 *= np.sign(v0[0])
gap = ev[-1] - ev[-2]
print(f"[1] lemma (**) vs truth, N={N}: g2={g2:.12f} gap={gap:.6f}")
for w in (0.1, 0.5, 1.0, 1.7, 2.4, 4.0, 6.0):
    K = build_fast(N, fK(w)); kap = v0 @ K @ v0; Kv = K @ v0
    eh = np.sqrt(max(0, Kv @ Kv - kap ** 2)); Lam = -np.linalg.eigvalsh(K)[0]
    s = np.linspace(0, 1, 100001)
    Psi = -(1 - s ** 2) * kap + s ** 2 * (Lam - gap / w ** 2) + 2 * s * eh
    ub = g2 + w * w * Psi.max()
    h_w2 = np.linalg.eigvalsh(build_fast(N, fH(w)))[-1] / w ** 2
    print(f"  w={w:4}: maxPsi={Psi.max():+.6f}  UB h/w^2={ub:.9f}  truth h/w^2={h_w2:.9f}  UB-truth={ub-h_w2:+.3e}  {'OK' if ub >= h_w2 else 'VIOLATION'}")

print("\n[2] max Psi-bar from the certified cell numbers (g_lo = 1.30443620929)")
g_lo = 1.30443620929
cells = [(0.0, 0.5, 0.9398847, 1.102995, 0.398982, 0.889745),
         (0.75, 0.975, 0.6097235, 0.7095244, 0.248381, 0.165901),
         (1.6286, 1.7386, 0.3058664, 0.2878522, 0.06252, 0.0186593),
         (2.5153, 2.5981, 0.1623498, 0.141258, 0.0212458, 0.00950121),
         (5.8899, 5.973, 0.03539221, 0.02888696, 0.00204713, 0.00358102),
         (5.973, 6.0, 0.035137823, 0.027671568, 0.0012350233, 0.00609272)]
for a, b, kl, Lh, eh, margin in cells:
    s = np.linspace(0, 1, 200001)
    Pb = -kl + s ** 2 * (kl + Lh - g_lo / b ** 2) + 2 * s * eh
    print(f"  [{a},{b}]: max Psibar = {Pb.max():+.6f} at s={s[Pb.argmax()]:.3f}; script margin {margin}; w^2(Lam+2eh)@b = {b*b*(Lh+2*eh):.5f}")

print("\n[3] cell-uniform Lam_hi / kappa_lo / ehat_hi vs truth at LEFT and RIGHT cell ends (N=1500)")
for a, b, kl, Lh, eh, margin in cells[1:]:
    for w in (a, b):
        K = build_fast(N, fK(w)); kap = v0 @ K @ v0; Kv = K @ v0
        e = np.sqrt(max(0, Kv @ Kv - kap ** 2)); Lam = -np.linalg.eigvalsh(K)[0]
        print(f"  w={w:7.4f}: kappa={kap:.7f} >= kl {kl} : {kap>=kl};  Lam={Lam:.7f} <= Lh {Lh} : {Lam<=Lh};  ehat={e:.7f} <= eh {eh} : {e<=eh}")

print("\n[4] Lemma A on a 2-atom P (N=800): lam1(sum p_i H(w_i)) <= sum p_i h(w_i)")
N = 800
for (w1, w2, p) in ((0.3, 2.4, 0.5), (1.0, 6.0, 0.2), (0.5, 0.7, 0.9)):
    H1 = build_fast(N, fH(w1)); H2 = build_fast(N, fH(w2))
    lhs = np.linalg.eigvalsh(p * H1 + (1 - p) * H2)[-1]
    rhs = p * np.linalg.eigvalsh(H1)[-1] + (1 - p) * np.linalg.eigvalsh(H2)[-1]
    print(f"  w=({w1},{w2}) p={p}: lam1(avg)={lhs:.6f} <= avg(h)={rhs:.6f} : {lhs<=rhs};  g2*<w^2>={g2*(p*w1**2+(1-p)*w2**2):.6f}")
