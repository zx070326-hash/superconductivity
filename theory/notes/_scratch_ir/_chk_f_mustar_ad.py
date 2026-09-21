"""Is our Coulomb deficit D really 2-3x 'the' weak-coupling prediction, or was the
comparison made against the wrong baseline?  Main session, independent of the agent.

Three baselines, all closed-form (no matrices):
  BCS/Morel-Anderson : Tc ~ exp(-1/(lam - mus))              -> what the agent used
  Allen-Dynes 1975   : Tc ~ exp(-1.04(1+lam)/(lam - mus(1+0.62 lam)))
  AD with our cutoff : same, but mus first run from our r=10 cutoff down to the
                       phonon scale via 1/mus(r) = 1/mus(10) - ln(r/10), r=1.
D := 1 - Tc(mus)/Tc(0).
"""
import math

def D_bcs(lam, mus):
    return 1.0 - math.exp(-1.0/(lam-mus) + 1.0/lam) if lam > mus else float('nan')

def D_ad(lam, mus):
    d = lam - mus*(1.0 + 0.62*lam)
    if d <= 0: return float('nan')
    return 1.0 - math.exp(-1.04*(1+lam)/d + 1.04*(1+lam)/lam)

def run_mus(mus10, r):
    return 1.0/(1.0/mus10 - math.log(r/10.0))

# (lam, mus, D_obs%) straight off the agent's table, which I am treating as DATA
obs = [(0.5,0.05,36.453),(0.5,0.10,60.081),(0.5,0.13,69.972),(0.5,0.15,75.230),
       (0.7,0.05,24.810),(0.7,0.10,43.088),(0.7,0.13,51.656),(0.7,0.20,66.548),
       (1.0,0.05,17.036),(1.0,0.10,30.375),(1.0,0.13,36.884),(1.0,0.20,49.296),
       (1.5,0.05,12.231),(1.5,0.10,21.885),(1.5,0.13,26.889),(1.5,0.20,36.334),
       (2.0,0.05,10.369),(2.0,0.10,18.603),(2.0,0.13,22.844),(2.0,0.20,30.987),
       (3.0,0.05, 8.837),(3.0,0.10,15.938),(3.0,0.13,19.488),(3.0,0.20,26.614),
       (5.0,0.05, 7.712),(5.0,0.10,14.000),(5.0,0.13,17.327),(5.0,0.20,23.727)]

print(f"{'lam':>5}{'mus':>7}{'D_obs%':>9}{'D_BCS%':>9}{'D_AD%':>9}"
      f"{'obs/BCS':>9}{'obs/AD':>9}   mus@r=1 -> D_AD'%   obs/AD'")
rb=[]; ra=[]; rap=[]
for lam, mus, d in obs:
    b, a = 100*D_bcs(lam,mus), 100*D_ad(lam,mus)
    m1 = run_mus(mus,1.0); ap = 100*D_ad(lam,m1)
    rb.append(d/b); ra.append(d/a); rap.append(d/ap)
    print(f"{lam:5.1f}{mus:7.3f}{d:9.2f}{b:9.2f}{a:9.2f}{d/b:9.2f}{d/a:9.2f}"
          f"     {m1:7.4f} {ap:9.2f}{d/ap:9.2f}")
import statistics as st
for name, rr in (("obs/BCS",rb),("obs/AD",ra),("obs/AD(mus run to r=1)",rap)):
    print(f"\n{name:26s} mean {st.mean(rr):.3f}  sd {st.pstdev(rr):.3f}"
          f"  range [{min(rr):.3f}, {max(rr):.3f}]")

# --- isolate the mechanism: is it the (1+0.62 lam) enhancement doing the work? ---
def D_ad_noenh(lam, mus):          # AD prefactor/exponent but mus NOT enhanced
    d = lam - mus
    return 1.0 - math.exp(-1.04*(1+lam)/d + 1.04*(1+lam)/lam)
print("\n--- mechanism isolation: AD exponent WITHOUT the (1+0.62 lam) enhancement ---")
r0=[]
for lam, mus, d in obs:
    v = 100*D_ad_noenh(lam,mus); r0.append(d/v)
print(f"obs/AD_no_enhancement     mean {st.mean(r0):.3f}  sd {st.pstdev(r0):.3f}"
      f"  range [{min(r0):.3f}, {max(r0):.3f}]")
print("\nper-lambda ratio at mus=0.13 (drift diagnostic):")
print(f"{'lam':>5}{'obs/AD':>9}{'obs/AD_noenh':>14}")
for lam, mus, d in obs:
    if abs(mus-0.13) < 1e-9:
        print(f"{lam:5.1f}{d/(100*D_ad(lam,mus)):9.3f}{d/(100*D_ad_noenh(lam,mus)):14.3f}")
