# Main session 2026-09-29: archive the Richardson post-processing that the 09-27 note quoted but never saved.
# Reads Tc(N=40,60,90) triples from _out_hofmann_largeN_1_2026-09-27.txt (produced by _chk_sonnet_hofmann_largeN_1_*.py),
# extrapolates in 1/N^2 (2-pt on N=60,90; 3-pt with 1/N^4 term) and reports r = Tc/(C_inf sqrt(lam)).
import re, os, numpy as np
d = os.path.dirname(os.path.abspath(__file__))
C = 0.1827262477
Ns = np.array([40., 60., 90.])
for line in open(os.path.join(d, "_out_hofmann_largeN_1_2026-09-27.txt")):
    m = re.match(r"lam=\s*([\d.]+)\s+Tc\(N=40,60,90\)=\[(.*?)\]", line)
    if not m: continue
    lam = float(m.group(1)); T = np.array([float(x) for x in re.findall(r"float64\(([\d.eE+-]+)\)", m.group(2))])
    t2 = np.linalg.solve(np.c_[np.ones(2), Ns[1:]**-2.], T[1:])[0]
    t3 = np.linalg.solve(np.c_[np.ones(3), Ns**-2., Ns**-4.], T)[0]
    print(f"lam={lam:5.1f}  Tc_2pt={t2:.8f}  Tc_3pt={t3:.8f}  r_2pt={t2/(C*np.sqrt(lam)):.10f}  r_3pt={t3/(C*np.sqrt(lam)):.10f}")
