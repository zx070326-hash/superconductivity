# rev7: wbar_A with the eps1/eps2 formulas replaced by the ones the written chain (L3)+(L4) actually yields:
#   eps1_rig = q [ pi^2/12 + ln(1+Vw) + V/(V-D) ],   eps2_rig = 7z(3)D^4/576 + q^2 [ 4 ln(1+Vw) + 2 V^2/(V-D)^2 ]
import sys, numpy as np
exec(open(r"D:\physics\superconductivity\theory\notes\_scratch_ir\rev7_asm_a.py").read())
import p4_cii_m2 as m2
def patched(w, rig1=True, rig2=True):
    D = 1/w; q = D*D/8
    e1o, e2o = m2.eps1, m2.eps2
    def e1(nm, ww):
        V = VF[nm]; DD = 1/ww; qq = DD*DD/8
        return qq*(np.pi**2/12 + np.log(1+V*ww) + V/(V-DD)) if rig1 else e1o(nm, ww)
    def e2(nm, ww):
        V = VF[nm]; DD = 1/ww; qq = DD*DD/8
        return Z3H*DD**4/576 + qq*qq*(4*np.log(1+V*ww) + 2*V*V/(V-DD)**2) if rig2 else e2o(nm, ww)
    m2.eps1, m2.eps2 = e1, e2
    try:
        r = m2.bounds2(w)["ratio"]
    finally:
        m2.eps1, m2.eps2 = e1o, e2o
    return r
for lab, kw in (("original", dict(rig1=False, rig2=False)), ("eps1 rigorous chain", dict(rig1=True, rig2=False)), ("eps2 rigorous chain", dict(rig1=False, rig2=True)), ("both", dict(rig1=True, rig2=True))):
    f = lambda w, kw=kw: patched(w, **kw)
    print(f"  {lab:22s}: wbar_A = {bis(f):.4f}   ratio@8.5 = {f(8.5):.5f}  @10 = {f(10.0):.5f}  @45 = {f(45.0):.6f}")
