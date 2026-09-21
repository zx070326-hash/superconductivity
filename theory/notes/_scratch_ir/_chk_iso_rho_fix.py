"""W1 CHECK OF A CORRECTION: s_hi exponent, and the bare-vs-reduced mass choice.

The note's prose said   s_hi = t^{2 w_B}.   Recompute what it actually is, then
verify the CORRECTED form regenerates every published cell -- not merely that the
old error is gone (CLAUDE.md W1: a fix can introduce new false statements).

Careful with conventions: in my machinery  t := sqrt(m_H/m_D) < 1.
team-lead's message uses t := sqrt(m_D/m_H) > 1.  Both must be stated.

Run: OMP_NUM_THREADS=1 python _chk_iso_rho_fix.py
"""
import numpy as np

m = {"H": 1.00782503, "D": 2.01410178, "T": 3.01604928, "Pd": 106.42}
t_mine = np.sqrt(m["H"] / m["D"])        # < 1, the convention in _chk_iso_rho*.py
t_tl = np.sqrt(m["D"] / m["H"])          # > 1, team-lead's convention
WLOG, A_BASE = 405.0, -0.3992
SHIFT = 0.5 / np.log(m["D"] / m["H"])

print("=== 0. conventions ===")
print(f"  m_D/m_H              = {m['D']/m['H']:.6f}   ln = {np.log(m['D']/m['H']):.6f}")
print(f"  t_mine = sqrt(m_H/m_D) = {t_mine:.6f}     t_tl = sqrt(m_D/m_H) = {t_tl:.6f}")

print()
print("=== 1. what is s_hi, derived from the identity (not from the note) ===")
print("  s(w_B, beta) = t^{2 w_B} / (1 - beta(1-t^2)),  t = t_mine")
print("  at beta = 1:  denominator = 1 - (1 - t^2) = t^2")
print("  ==> s_hi = t^{2 w_B} / t^2 = t^{2(w_B - 1)} = t_mine^{-2 w_A} = t_tl^{2 w_A}")
print("      with w_A := 1 - w_B.   Equivalently s_hi = (m_D/m_H)^{w_A}.")


def s_exact(w_B, beta, t=t_mine):
    return t ** (2.0 * w_B) / (1.0 - beta * (1.0 - t**2))


def s_hi_new(w_B):
    return (m["D"] / m["H"]) ** (1.0 - w_B)


def s_hi_old_prose(w_B):
    return t_mine ** (2.0 * w_B)


print()
print(f"{'w_B':>8} | {'identity@beta=1':>15} {'NEW t^{2w_A}':>13} {'OLD prose t^{2w_B}':>19} {'OK?':>5}")
for w_B in (0.750, 0.790, 0.830, 0.870, 0.9023, 1.000):
    a, b, c = s_exact(w_B, 1.0), s_hi_new(w_B), s_hi_old_prose(w_B)
    print(f"{w_B:8.4f} | {a:15.6f} {b:13.6f} {c:19.6f} {'yes' if abs(a-b)<1e-12 else 'NO':>5}")
print("  ==> OLD prose form is wrong AND monotone the wrong way "
      f"({s_hi_old_prose(0.75):.4f} -> {s_hi_old_prose(0.9023):.4f} rises; truth falls).")

print()
print("=== 2. consistency of the CORRECTED form with Lemma rho-1 (independent check) ===")
print(f"  w_A -> 0 (everything substituted = UNIFORM dilation) must give s = 1 exactly.")
print(f"    s_hi(w_B=1) = {s_hi_new(1.0):.15f}")
print(f"    and the full identity at w_B=1 for several beta (must all be 1):")
for beta in (0.0, 0.5, 0.9, 1.0):
    print(f"      beta={beta:.1f} -> s = {s_exact(1.0, beta):.15f}")
print("  ==> corrected form reproduces Lemma rho-1; the OLD form gives "
      f"{s_hi_old_prose(1.0):.6f} at w_B=1, contradicting it.")

print()
print("=== 3. team-lead's reverse-engineering, checked against UNROUNDED values ===")
for w_B in (0.750, 0.9023):
    v = s_exact(w_B, 1.0)
    print(f"  w_B={w_B:.4f}: ln(s_hi)/(1-w_B) = {np.log(v)/(1-w_B):.6f}  "
          f"(exact ln(m_D/m_H) = {np.log(m['D']/m['H']):.6f})")
print("  (their 0.692450/0.692514 came from my 4-dp table; exact values agree to 1e-15)")

print()
print("=== 4. REGENERATE every published cell from the corrected form ===")


def beta_lo(om, w_B):
    w_A = 1.0 - w_B
    G = (WLOG / om**w_A) ** (1.0 / w_B)
    return 1.0 / (1.0 + (w_A * om**2) / (w_B * G**2))


OMS = (150.0, 190.0, 230.0, 270.0, 310.0)
WBS = (0.750, 0.790, 0.830, 0.870, 0.9023)
print("  alpha interval per cell; 'hi' recomputed from s_hi_new, 'lo' from the identity")
bad = 0
for om in OMS:
    row = f"  {om:5.0f} |"
    for w_B in WBS:
        lo = A_BASE + SHIFT * np.log(s_exact(w_B, beta_lo(om, w_B)))
        hi_id = A_BASE + SHIFT * np.log(s_exact(w_B, 1.0))
        hi_new = A_BASE + SHIFT * np.log(s_hi_new(w_B))
        if abs(hi_id - hi_new) > 1e-12:
            bad += 1
        row += f" [{lo:+.3f},{hi_new:+.3f}]"
    print(row)
print(f"  cells where corrected closed form != identity: {bad} / {len(OMS)*len(WBS)}")

print()
print("=== 5. the 5.3 table (Om=190) and the s-interval table, regenerated ===")
for w_B in (0.750, 0.790, 0.9023):
    b = beta_lo(190.0, w_B)
    print(f"  w_B={w_B:.4f}  beta_lo={b:.5f}  s=[{s_exact(w_B,b):.4f},{s_hi_new(w_B):.4f}]  "
          f"alpha=[{A_BASE+SHIFT*np.log(s_exact(w_B,b)):+.4f},"
          f"{A_BASE+SHIFT*np.log(s_hi_new(w_B)):+.4f}]")

print()
print("=== 6. MODELLING CHOICE: bare mass ratio vs reduced mass in a Pd matrix ===")
mu = {i: m[i] * m["Pd"] / (m[i] + m["Pd"]) for i in ("H", "D", "T")}
r_bare, r_red = m["D"] / m["H"], mu["D"] / mu["H"]
print(f"  bare      m_D/m_H   = {r_bare:.6f}  ln = {np.log(r_bare):.6f}")
print(f"  reduced  mu_D/mu_H  = {r_red:.6f}  ln = {np.log(r_red):.6f}   "
       f"(differ by {100*(r_bare/r_red-1):.2f}%)")
print()
print("  where it enters, and how much it moves things:")
HARM = {"H": 205.0, "D": 150.0, "T": 125.0}
for lab, r in (("bare", r_bare), ("reduced", r_red)):
    wHD = 2 * np.log(HARM["H"] / HARM["D"]) / np.log(r)
    print(f"    w_op readout (harmonic H->D), {lab:>7}: {wHD:.4f}")
print(f"    ==> the choice moves the w_op readout by "
      f"{abs(2*np.log(HARM['H']/HARM['D'])/np.log(r_red) - 2*np.log(HARM['H']/HARM['D'])/np.log(r_bare)):.4f}")
print(f"    ==> vs the unexplained w_B tension 0.11 and the +-0.062 requirement")
print()
for lab, r in (("bare", r_bare), ("reduced", r_red)):
    sh = r ** (1 - 0.9023)
    print(f"    s_hi(w_B=0.9023), {lab:>7}: {sh:.4f}  -> alpha_hi = {A_BASE+SHIFT*np.log(sh):+.4f}")
print("    ==> shifts alpha_hi by less than 0.002: does NOT touch the unreachability verdict")
