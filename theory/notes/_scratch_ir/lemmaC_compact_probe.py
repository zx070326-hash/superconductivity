"""
Lemma C compact-interval probe (float64, numerical-evidence only).

Independent rebuild of the kernel-generic operator O[phi] (no import of any
gate script).  Measures, on w in (0, 8.5]:

  (1) the true (C**) decomposition  H'' = -2F + 2 x^2 <r,(h-A)^+ r>  and the
      Davis-Kahan-replaced positive term 2 x^2 E^2 / g, i.e. the true
      "margin" pos/F that any analytic proof must beat;
  (2) the q-route decomposition  q' = E_q + C,  |C|/E_q;
  (3) how well explicit trial vectors capture v: sin(angle), residual/gap,
      and the resulting DK-propagated (C**) ratio when the gap is replaced by
      the best *analytic* gap lower bounds available in the repo
      (0.4h on (0,2], h - Theta* for w>=2);
  (4) the relative size of the second-order term with A's exact reduced
      resolvent (to see how much DK loses).

Run:  OMP_NUM_THREADS=1 python lemmaC_compact_probe.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

THETA_STAR = 0.6349287917881079


def O_of(phi, N):
    """(O[phi])_{nm} = u_n u_m [phi(|n-m|) 1_{n!=m} + phi(n+m+1)] - delta_nm (2/(2n+1)) sum_{k<=n} phi(k)."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    K = np.arange(0, 2 * N + 2, dtype=float)
    pk = phi(K)
    pk[0] = 0.0
    M = (pk[D] + pk[S]) * np.outer(u, u)
    cum = np.cumsum(pk[:N])          # sum_{k=1}^n phi(k)  (pk[0]=0)
    M[np.diag_indices(N)] -= 2.0 * cum / (2.0 * n + 1.0)
    return M


def top(M):
    w, V = np.linalg.eigh(M)
    v = V[:, -1]
    if v[0] < 0:
        v = -v
    return w[-1], v, w[-2], w, V


def probe(w, N=600):
    x = w * w
    f = lambda k: x / (k * k + x)
    A = O_of(f, N)
    h, v, l2, ev, EV = top(A)
    g = h - l2
    # (C**) pieces
    phi = lambda k: k * k / (k * k + x) ** 3
    psi = lambda k: 1.0 / (k * k + x) ** 2
    F = v @ O_of(phi, N) @ v
    r = O_of(psi, N) @ v
    r = r - v * (v @ r)
    E = np.linalg.norm(r)
    # exact reduced resolvent term 2 x^2 <r,(h-A)^+ r>
    c = EV.T @ r
    c[-1] = 0.0
    exact_pos = 2 * x * x * np.sum(c[:-1] ** 2 / (h - ev[:-1]))
    dk_pos = 2 * x * x * E * E / g
    Hpp_true = -2 * F + exact_pos
    # q-route
    f2 = lambda k: f(k) ** 2
    gk = lambda k: 2 * x * k * k / (k * k + x) ** 3
    Eq = v @ O_of(gk, N) @ v
    B = O_of(f2, N)
    Bv = B @ v
    q = v @ Bv
    Bp = Bv - q * v
    cb = EV.T @ Bp
    cb[-1] = 0.0
    Cq = -(2.0 / x) * np.sum(cb[:-1] ** 2 / (h - ev[:-1]))
    return dict(w=w, x=x, h=h, l2=l2, g=g, F=F, E=E, exact_pos=exact_pos,
                dk_pos=dk_pos, Hpp=Hpp_true, ratio_true=exact_pos / (2 * F),
                ratio_dk=dk_pos / (2 * F), Eq=Eq, Cq=Cq, CoverE=-Cq / Eq,
                v=v, A=A, N=N)


def trial_dk(pr, t, glo):
    """DK propagation for an explicit trial vector t with analytic gap glo.
    Returns sin bound, and the (C**) ratio bound  x^2 E_hi^2/(g_lo F_lo)."""
    A, v, x, N = pr["A"], pr["v"], pr["x"], pr["N"]
    t = t / np.linalg.norm(t)
    rho = t @ A @ t
    res = np.linalg.norm(A @ t - rho * t)
    s_true = np.sqrt(max(0.0, 1 - (t @ v) ** 2))
    s_bnd = res / glo if glo > 0 else np.inf
    phi = lambda k: k * k / (k * k + x) ** 3
    psi = lambda k: 1.0 / (k * k + x) ** 2
    Ophi, Opsi = O_of(phi, N), O_of(psi, N)
    Ft = t @ Ophi @ t
    e_phi = np.linalg.norm(Ophi @ t - Ft * t)
    Kt = Opsi @ t
    qpsi = t @ Kt
    e_psi = np.linalg.norm(Kt - qpsi * t)
    nphi = np.abs(np.linalg.eigvalsh(Ophi)).max()
    npsi = np.abs(np.linalg.eigvalsh(Opsi)).max()
    s = min(s_bnd, 1.0)
    F_lo = (1 - s * s) * Ft - 2 * s * e_phi - s * s * nphi
    E_hi = e_psi + s * (abs(qpsi) + npsi)
    ratio = x * x * E_hi ** 2 / (glo * F_lo) if (F_lo > 0 and glo > 0) else np.inf
    return dict(rho=rho, res=res, s_true=s_true, s_bnd=s_bnd, F_lo=F_lo, E_hi=E_hi, ratio=ratio,
                Ft=Ft, e_psi=e_psi)


if __name__ == "__main__":
    ws = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 8.5, 10.0]
    print("== (1)(2)(4): true structure of the two routes  (N=600) ==")
    print(f"{'w':>5} {'h':>8} {'g':>8} {'g/h':>6} {'F':>10} {'E':>10} {'pos_exact/2F':>13} {'pos_DK/2F':>10} {'DK/exact':>9} {'|C|/E_q':>8} {'H''':>12}")
    rows = {}
    for w in ws:
        p = probe(w)
        rows[w] = p
        print(f"{w:5.2f} {p['h']:8.4f} {p['g']:8.4f} {p['g']/p['h']:6.3f} {p['F']:10.4e} {p['E']:10.4e} "
              f"{p['ratio_true']:13.4f} {p['ratio_dk']:10.4f} {p['dk_pos']/p['exact_pos']:9.3f} {p['CoverE']:8.4f} {p['Hpp']:12.4e}")

    print("\n== (3): explicit trial vectors + analytic gap bounds (0.4h on (0,2], h-Theta* for w>=2; both need h_lo, here true h) ==")
    print("trial families: T1 lorentz tau_n=p_n G((n+1/2)/w); T2 u_n/(1+(n/w)^2); T3 e0; T4 G2-like t_n=u_n/(n+1)^a best a")
    print(f"{'w':>5} {'g_an':>7} | {'T1 s_true':>9} {'s_bnd':>7} {'ratio':>8} | {'T2 s_true':>9} {'s_bnd':>7} {'ratio':>8} | {'e0 s_true':>9} {'s_bnd':>7} {'ratio':>8}")
    for w in ws:
        p = rows[w]
        h = p["h"]
        g_an = max(0.4 * h, h - THETA_STAR) if w <= 2 else max(h - THETA_STAR, 0.4 * h if w <= 2 else 0)
        if w > 2:
            g_an = h - THETA_STAR
        N = p["N"]
        n = np.arange(N)
        u = 1 / np.sqrt(2 * n + 1)
        t1 = u / (1 + ((n + 0.5) / w) ** 2)
        t2 = u / (1 + (n / w) ** 2)
        t3 = np.zeros(N); t3[0] = 1.0
        out = []
        for t in (t1, t2, t3):
            d = trial_dk(p, t, g_an)
            out.append(d)
        print(f"{w:5.2f} {g_an:7.3f} | " + " | ".join(f"{d['s_true']:9.4f} {d['s_bnd']:7.3f} {d['ratio']:8.3f}" for d in out))
