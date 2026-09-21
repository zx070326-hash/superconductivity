"""
验证 h^(d)(varpi) ~ C^(d) * ln(varpi) + O(1) 在 d=1,2,3 的系数 C^(d)。

扩展自 dim_construct_dXX.py（不改基线），复用其 tau_n^(d)/build_K/h_of。
严格性：整条继承 dim_construct_dXX.py 的 numerical-evidence 级物理建模——
tau_n^(1), tau_n^(2) 里的 KAPPA1, KAPPA2, LAMBDA 是未定标 O(1) 建模常数。
d=1 的 tau_n^(1) = KAPPA1 * Lambda^{1/2} * (2n+1)^{-3/2} 带有一个不随 n 变化、
仅依赖 Lambda 的整体标度因子 Lambda^{1/2}=31.6...，这会把 h^(1) 整体放大，
使拟合出的 C^(1) 强烈依赖 Lambda 的取值——这是建模假设的产物，不是普适常数。
下面报告里会显式标注这一点，不把 C^(1) 当作可与 C^(2),C^(3) 直接比较的「同一把尺子」。

M 收敛性（varpi=10 处独立测过，见会话记录）：
  d=3: M=120 -> 2.873516, M=150 应更准；d=1: M=150 vs 180 差 ~7e-6；d=2: M=150 vs 180 差 ~9e-6。
  取 M=150 对本窗口 varpi<=10 足够 4-5 位有效数字。
"""
import time
import sys
sys.path.insert(0, '.')
import dim_construct_dXX as dc
from mpmath import mp

mp.dps = 30
M = 150
WINDOW = [5, 6, 7, 8, 9, 10]


def compute_all():
    results = {1: [], 2: [], 3: []}
    t0 = time.time()
    for d in (1, 2, 3):
        for w in WINDOW:
            t1 = time.time()
            h = dc.h_of(w, d, M=M)
            results[d].append((w, float(h)))
            print(f"d={d} varpi={w:2d} h={mp.nstr(h, 14)}  ({time.time()-t1:.1f}s)", flush=True)
    print("TOTAL", round(time.time() - t0, 1), "s")
    return results


def fit_2param(pts):
    """h = C ln(w) + D，普通最小二乘（2 参数，全窗口线性拟合，会被 O(1/ln w) 项污染）。"""
    import math
    xs = [math.log(w) for w, _ in pts]
    ys = [h for _, h in pts]
    n = len(xs)
    sx = sum(xs); sy = sum(ys); sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, ys))
    C = (n*sxy - sx*sy) / (n*sxx - sx*sx)
    D = (sy - C*sx) / n
    # residuals
    resid = [y - (C*x + D) for x, y in zip(xs, ys)]
    rms = (sum(r*r for r in resid)/n) ** 0.5
    return C, D, rms


def fit_3param(pts):
    """h = C ln(w) + c0 + b1/ln(w)，最小二乘（3 参数），用来吸收 O(1/ln w) 修正，
    减小把子领先项误判进斜率 C 的系统偏差——这是本任务里代替「varpi->infty 直接算」的
    Richardson 式外推（同一功能形式的多点拟合去掉子领先项对 C 的污染）。"""
    import math
    # design matrix columns: [ln w, 1, 1/ln w]
    rows = []
    ys = []
    for w, h in pts:
        lw = math.log(w)
        rows.append([lw, 1.0, 1.0/lw])
        ys.append(h)
    # normal equations A^T A x = A^T y
    import itertools
    n = 3
    ATA = [[sum(rows[k][i]*rows[k][j] for k in range(len(rows))) for j in range(n)] for i in range(n)]
    ATy = [sum(rows[k][i]*ys[k] for k in range(len(rows))) for i in range(n)]
    # solve 3x3 via Cramer's rule
    def det3(m):
        return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
                - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
                + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    D0 = det3(ATA)
    sol = []
    for col in range(3):
        Mc = [row[:] for row in ATA]
        for r in range(3):
            Mc[r][col] = ATy[r]
        sol.append(det3(Mc)/D0)
    C, c0, b1 = sol
    resid = []
    for w, h in pts:
        lw = math.log(w)
        pred = C*lw + c0 + b1/lw
        resid.append(h - pred)
    rms = (sum(r*r for r in resid)/len(resid))**0.5
    return C, c0, b1, rms


if __name__ == '__main__':
    results = compute_all()
    print()
    print("=== 2-参数线性拟合 h=C ln(w)+D，窗口 varpi in [5,10] ===")
    for d in (1, 2, 3):
        C, D, rms = fit_2param(results[d])
        print(f"d={d}: C={C:.4f}  D={D:.4f}  rms={rms:.2e}")

    print()
    print("=== 3-参数拟合 h=C ln(w)+c0+b1/ln(w)（同一窗口内做子领先项修正） ===")
    for d in (1, 2, 3):
        C, c0, b1, rms = fit_3param(results[d])
        print(f"d={d}: C={C:.4f}  c0={c0:.4f}  b1={b1:.4f}  rms={rms:.2e}")

    print()
    print("=== 子窗口交叉检验（Richardson 式）：只用后4点 varpi in [7,10] 做 2 参数拟合 ===")
    for d in (1, 2, 3):
        sub = [p for p in results[d] if p[0] >= 7]
        C, D, rms = fit_2param(sub)
        print(f"d={d}: C(sub[7,10])={C:.4f}  rms={rms:.2e}")

    print()
    print("=== d=3 基线自检（已知解析值 C^(3)=1，因 h=ln w + c0 + b1/ln w + ...）===")
