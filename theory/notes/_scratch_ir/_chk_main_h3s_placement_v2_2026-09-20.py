# Corrected H3S/D3S placement on the AD75 trend, using COMPLETE self-consistent
# (lambda, w_log, wbar2) triples integrated from ONE alpha2F each (Mishra Figshare).
import math
CINF, K_PER_MEV = 0.1827262477, 11.604518
A,B,C = -1.0691, 0.9972, 0.7198          # AD75 regression; band +-7.94%
BAND = 0.0794

def row(name, lam, wlog, wb2, tc):
    ceil = CINF*math.sqrt(lam)*wb2*K_PER_MEV
    shape = wlog/wb2
    fp = math.exp(A + B*math.log(lam) + C*math.log(shape))
    fo = tc/ceil
    r  = fo/fp - 1.0
    # Tc at which the residual hits the lower band edge
    tc_edge = tc*(1-BAND)/(1+r)
    print(f"  {name:22s} lam={lam:.4f} wlog={wlog:7.3f} wb2={wb2:7.3f} shape={shape:.4f} "
          f"ceil={ceil:7.2f}K f_obs={fo:.4f} f_pred={fp:.4f} resid={r*100:+6.2f}% "
          f"{'PASS' if abs(r)<BAND else '**FAIL**'}   Tc@edge={tc_edge:6.2f}K (margin {tc-tc_edge:+5.2f}K)")
    return r

# moments recomputed by me from the raw .a2f (_chk_main_h3s_a2f_moments_2026-09-20.py)
H = dict(lam=1.720958, wlog=109.995163, wb2=149.262922)
Dd= dict(lam=1.783761, wlog= 85.535326, wb2=105.430104)

print("=== CORRECTED: complete consistent triples, Mishra alpha2F @200 GPa ===")
row("H3S Tc=184.04(200GPa)", H['lam'],H['wlog'],H['wb2'],184.04)
row("D3S Tc=150  (P UNVERIF)",Dd['lam'],Dd['wlog'],Dd['wb2'],150.0)

print("\n=== what the earlier rows had wrong ===")
row("agent H3S Tc=203@155GPa", H['lam'],H['wlog'],H['wb2'],203.0)
row("agent D3S wlog=1/sqrt2 ", Dd['lam'],109.995163/math.sqrt(2),Dd['wb2'],150.0)
row("my  D3S wlog=Errea ratio",Dd['lam'],109.995163*(73.3/92.86),Dd['wb2'],150.0)

print("\n=== D3S Tc sensitivity (pressure pairing never checked) ===")
for tc in (150,148,146,144,142,140):
    row(f"D3S Tc={tc}", Dd['lam'],Dd['wlog'],Dd['wb2'],float(tc))

print("\n=== isotope-scaling diagnostics from the two a2F pairs ===")
print(f"  anharmonic: lam_D/lam_H = {1.783761/1.720958:.6f}   (uniform mass scaling requires exactly 1)")
print(f"  harmonic  : lam_D/lam_H = {2.293767/2.292857:.6f}   <- invariant to 0.04%, as Hopfield requires")
print(f"  anh wbar2 ratio = {105.430104/149.262922:.6f}   w_log ratio = {85.535326/109.995163:.6f}   1/sqrt2={2**-0.5:.6f}")
print(f"  har wbar2 ratio = { 92.973347/129.314988:.6f}   w_log ratio = {79.561668/107.165325:.6f}")
