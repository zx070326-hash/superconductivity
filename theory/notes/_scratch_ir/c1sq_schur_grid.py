import sys
sys.path.insert(0, 'theory/notes/_scratch_ir')
from c1sq_schur_bound import run_one_schur

grid = [0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
for w in grid:
    run_one_schur(w, N)
    print()
