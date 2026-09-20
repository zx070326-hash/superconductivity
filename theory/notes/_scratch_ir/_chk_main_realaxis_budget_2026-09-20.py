# Budget for the real-axis correction c := Delta_0 / Delta(i w_0).
# Consumption point: the verdict R > 1 (lambda-axis signal beats mu* interference).
# Pure algebra + Carbotte closed form. No Eliashberg solve.
import math

KB = 0.08617333262  # meV/K
WLOG = 106.9        # meV (Einstein => w_log = w_E)

# --- measured, self-consistent solver.tc() branch ---
# source: _chk_main_gapratio_tcsrc_2026-09-20.py
TC = {'mu+':170.807,'mu-':209.333,'lam+':213.481,'lam-':160.660}
OME= {'mu+':4.5278 ,'mu-':4.7346 ,'lam+':4.8231 ,'lam-':4.4078 }

def carbotte(tc):
    t = KB*tc/WLOG
    return 3.53*(1.0 + 12.5*t*t*math.log(1.0/(2.0*t)))

OC = {k: carbotte(v) for k,v in TC.items()}

half = lambda d,a,b: 0.5*(d[a]-d[b])
D_lam_ME = half(OME,'lam+','lam-');  D_mu_ME = half(OME,'mu+','mu-')
D_lam_C  = half(OC ,'lam+','lam-');  D_mu_C  = half(OC ,'mu+','mu-')
O0_ME = 0.25*sum(OME.values());      O0_C = 0.25*sum(OC.values())

print("recompute Eliashberg half-signals (must match note)")
print("  dO_lam=%+.5f  dO_mu=%+.5f   R=%.4f   O0=%.4f"%(D_lam_ME,D_mu_ME,abs(D_lam_ME/D_mu_ME),O0_ME))
print("recompute Carbotte on the SAME self-consistent Tc values")
print("  dO_lam=%+.5f  dO_mu=%+.5f   R=%.4f   O0=%.4f"%(D_lam_C,D_mu_C,abs(D_lam_C/D_mu_C),O0_C))
for k in TC: print("    %-5s Tc=%8.3f  O_ME=%.4f  O_Carb=%.4f  c=%.4f"%(k,TC[k],OME[k],OC[k],OC[k]/OME[k]))

# --- budget: O_real = c*O  =>  S_p = c0*D_p + O0*dc_p ;  r_p := dc_p/c0
# R_real(r_lam,r_mu) = |D_lam + O0*r_lam| / |D_mu + O0*r_mu|
def R_real(rl,rm): return abs(D_lam_ME+O0_ME*rl)/abs(D_mu_ME+O0_ME*rm)

print("\n=== budget 1: full cancellation of the lambda-axis signal ===")
rl_kill = -D_lam_ME/O0_ME
print("  need r_lam = %+.5f  (c changes %+.2f%% across the half-swing Dlam/lam=0.15)"%(rl_kill,100*rl_kill))
print("  i.e. dln c/dln lambda = %+.4f"%(rl_kill/0.15))

print("\n=== budget 2: single-channel flip thresholds R_real = 1 ===")
for s in (+1,-1):
    rl = (s*abs(D_mu_ME)-D_lam_ME)/O0_ME
    print("  r_mu=0 : r_lam = %+.5f  (%+.2f%%)"%(rl,100*rl))
for s in (+1,-1):
    rm = (s*abs(D_lam_ME)-D_mu_ME)/O0_ME
    print("  r_lam=0: r_mu  = %+.5f  (%+.2f%%)"%(rm,100*rm))

print("\n=== empirical r from Carbotte (the only third-party real-axis formula) ===")
rl_emp = D_lam_C/O0_C - D_lam_ME/O0_ME
rm_emp = D_mu_C /O0_C - D_mu_ME /O0_ME
print("  r_lam_emp = %+.5f (%+.2f%%)   r_mu_emp = %+.5f (%+.2f%%)"%(rl_emp,100*rl_emp,rm_emp,100*rm_emp))
print("  R_real at empirical r = %.4f   (cross-check vs R_Carb = %.4f)"%(R_real(rl_emp,rm_emp),abs(D_lam_C/D_mu_C)))

print("\n=== how much stronger than Carbotte before the verdict flips ===")
lo,hi=0.0,50.0
for _ in range(200):
    mid=0.5*(lo+hi)
    if R_real(mid*rl_emp,mid*rm_emp)>1.0: lo=mid
    else: hi=mid
print("  scale factor s such that R_real(s*r_emp) = 1 :  s = %.3f"%(0.5*(lo+hi)))
print("  single-channel margins:  lambda %.2fx   mu* %.2fx"%(
      abs((abs(D_mu_ME)-D_lam_ME)/O0_ME/rl_emp), abs((-abs(D_lam_ME)-D_mu_ME)/O0_ME/rm_emp)))

print("\n=== sensitivity of the budget to the O0 proxy (central O not measured) ===")
for f in (0.90,0.95,1.00,1.05,1.10):
    O0=O0_ME*f
    s_lo,s_hi=0.0,50.0
    for _ in range(200):
        m=0.5*(s_lo+s_hi)
        v=abs(D_lam_ME+O0*m*rl_emp)/abs(D_mu_ME+O0*m*rm_emp)
        if v>1.0: s_lo=m
        else: s_hi=m
    print("  O0=%.4f (%.0f%%)  ->  s = %.3f"%(O0,100*f,0.5*(s_lo+s_hi)))
