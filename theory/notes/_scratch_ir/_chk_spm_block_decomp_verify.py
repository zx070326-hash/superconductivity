"""
独立复核 team-lead 的块分解论证（未 import team-lead 的脚本，只 import 本文件同目录下
已修好的 _chk_spm_within_P_construction.py 取其 build_O/einstein，不复用其对外结论）。

验证：O = I_2 (x) A - sigma_x (x) B, A=T[Lam^ph]-L1[R], B=T[Lam^sf]
=> top eigenvalue of O 恰等于 max(lambda1(A+B), lambda1(A-B))
且当 B>=0（引理2, f in P）时 lambda1(A+B) >= lambda1(A-B)，顶本征向量必属反对称(s±)扇区。
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import importlib
mod = importlib.import_module('_chk_spm_within_P_construction')

def u(n):
    return (2*n+1)**-0.5

def Tmat(N, f):
    M = np.zeros((N, N))
    for n in range(N):
        for m in range(N):
            M[n, m] = u(n)*u(m)*(f(abs(n-m)) + f(n+m+1))
    return M

def L1mat(N, f):
    d = np.zeros(N)
    for n in range(N):
        s = f(0) + 2*sum(f(k) for k in range(1, n+1))
        d[n] = u(n)*u(n)*s
    return np.diag(d)

def check(N):
    P = lambda k: mod.einstein(k, w=1.0, amp=1.6)
    S = lambda k: mod.einstein(k, w=0.6, amp=2.6)
    def R(k):
        return P(k) + S(k)

    A = Tmat(N, P) - L1mat(N, R)
    B = Tmat(N, S)

    lam_Aplus = np.max(np.linalg.eigvalsh(A + B))
    lam_Aminus = np.max(np.linalg.eigvalsh(A - B))
    lam_Bmin = np.min(np.linalg.eigvalsh(B))

    def Lam_ph(i, j, k):
        return P(k) if i == j else 0.0
    def Lam_sf(i, j, k):
        return S(k) if i != j else 0.0
    O = mod.build_O(N, Lam_ph, Lam_sf, 2)
    w = np.linalg.eigvalsh(O)

    eigA_plus = np.linalg.eigvalsh(A + B)
    eigA_minus = np.linalg.eigvalsh(A - B)
    union_spec = np.sort(np.concatenate([eigA_plus, eigA_minus]))
    full_spec = np.sort(w)
    spec_match = np.max(np.abs(union_spec - full_spec)) < 1e-8

    print(f"N={N}: lam1(A+B)={lam_Aplus:.8f}  lam1(A-B)={lam_Aminus:.8f}  "
          f"lam_min(B)={lam_Bmin:.6e}  full_O_top2=({w[-1]:.8f},{w[-2]:.8f})")
    ok_match = abs(w[-1] - lam_Aplus) < 1e-8 and spec_match
    ok_Bpsd = lam_Bmin > -1e-10
    ok_order = lam_Aplus >= lam_Aminus - 1e-10
    print(f"   top-1(O)==lambda1(A+B): {abs(w[-1]-lam_Aplus)<1e-8}   全谱(O)==eig(A+B)∪eig(A-B): {spec_match}   "
          f"B半正定: {ok_Bpsd}   lam1(A+B)>=lam1(A-B): {ok_order}")
    return ok_match, ok_Bpsd, ok_order

if __name__ == "__main__":
    all_ok = True
    for N in [40, 60, 100, 200, 400]:
        r = check(N)
        all_ok = all_ok and all(r)
    print("ALL PASS" if all_ok else "FAIL")
