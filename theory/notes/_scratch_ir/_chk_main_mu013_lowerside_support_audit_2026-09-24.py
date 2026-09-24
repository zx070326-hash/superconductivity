"""同一支撑链检验施加到**已入库** mu*=0.13 表的下侧端点（gateM-cells 日志原值）。"""
import os, importlib.util; os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
_d = os.path.join("theory","notes","_scratch_ir")
def _load(n):
    s = importlib.util.spec_from_file_location(n, os.path.join(_d, n+".py"))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
KT = _load("_chk_main_kato_temple_two_atom_prototype_2026-09-23")
ROWS = [(3.2953, 0.668045, 6, 1.011682278), (6.0, 0.456360, 4, 1.013784411),
        (10.0,   0.339065, 3, 1.010339659)]
print(f"{'lam':>8}{'w_lo':>10}{'Nc':>4} | " + "".join(f"{'P_'+str(n):>15}" for n in (40,160,600))
      + f"{'mu_lo':>15}{'判':>22}")
bad = 0
for lam, w, nc, ml in ROWS:
    assert int(np.floor(10*w)) == nc, (lam, w, nc, int(np.floor(10*w)))
    vs = [KT.lam_k([lam],[w],0.13,nc,n) for n in (40,160,600)]
    ok = ml <= vs[0] + 1e-12; bad += not ok
    v = "OK (<=P_40)" if ok else ("** > P_600 fail-open **" if ml > vs[-1] else "** > P_40 **")
    print(f"{lam:8}{w:10.6f}{nc:4d} | " + "".join(f"{x:15.9f}" for x in vs)
          + f"{ml:15.9f}{v:>22}  (mu_lo-P_600 {ml-vs[-1]:+.2e})")
print(f"\n==> mu*=0.13 已入库下侧端点：{3-bad}/3 通过支撑链检验")
