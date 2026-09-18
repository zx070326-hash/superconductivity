# (1) monotonicity of r(w)=h(w)/w^2 on a fine grid; (2) second-order coefficient at w->0:
#     A(w)=H(w)/w^2 = G(2) - w^2 G(4) + O(w^4)  =>  r(w) = g2 - w^2 <v0,G(4)v0> + ...
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "derivations"))
from p4_interchange_uniformity import H_mat, G_matrix
N = 600
G2 = G_matrix(N, 2.0); G4 = G_matrix(N, 4.0)
e, V = np.linalg.eigh(G2); g2 = e[-1]; v0 = V[:, -1]
c = float(v0 @ G4 @ v0)
print(f"g2={g2:.6f} (N={N}), lam2(G2)={e[-2]:.2e}, rho(G4)={np.linalg.eigvalsh(G4)[-1]:.5f}, <v0,G(4)v0>={c:.5f} (>0 needed)")
ws = np.concatenate([np.geomspace(0.01, 0.3, 30), np.linspace(0.3, 3.0, 55)])
r = []
for w in ws:
    Nw = 400 if w < 1 else 1200
    r.append(np.linalg.eigvalsh(H_mat(Nw, w))[-1] / w**2)
r = np.array(r)
print("r(w)=h/w^2 monotone decreasing on grid:", bool(np.all(np.diff(r) < 0)), " max r =", r.max(), "at w =", ws[np.argmax(r)])
for w, rr in list(zip(ws, r))[:8]:
    print(f"  w={w:.4f}  r={rr:.6f}  (g2 - r)/w^2 = {(g2-rr)/w**2:.4f}   [2nd-order coeff predicts {c:.4f}]")
