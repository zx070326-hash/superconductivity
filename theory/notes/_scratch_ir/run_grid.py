import sys
from c1sq_lower_bound import run_one_dk

grid = [0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
for w in grid:
    run_one_dk(w, N)
    print()
