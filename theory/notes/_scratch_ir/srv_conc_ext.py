# Server wrapper: run the concavity certificate beyond wbar=45 with Na scaled ~14*wbar (cap 1000),
# monkeypatching na_band only (no edit to the committed script).  Evidence for how far C-I can go.
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_extremal_concavity as C
def na_band(w):
    if w <= 5.0: return 70
    if w <= 12.0: return 140
    if w <= 20.0: return 280
    return min(1000, int(math.ceil(14.0 * float(w) / 20.0) * 20))
C.na_band = na_band
print("[wrapper] na_band patched: Na=ceil(14 w/20)*20 for w>20 (cap 1000)", flush=True)
sys.argv = ["p4_extremal_concavity.py", "EXT", sys.argv[1] if len(sys.argv) > 1 else "120"]
C.main()
