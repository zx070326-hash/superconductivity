import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
np.random.seed(0)

# ---- DEFECT 1: does Rayleigh-on-point-eigenvectors enclose the spectrum?
# Take a symmetric matrix, perturb its eigenvectors the way float64 eigh would,
# and compare R(v_k) to the true mu_k.
N=60
B=np.random.randn(N,N); A=(B+B.T)/2
ev,V=np.linalg.eigh(A)
R=np.array([ (V[:,k]@A@V[:,k])/(V[:,k]@V[:,k]) for k in range(N)])
print("exact eigenvectors: max|R-mu| =", np.abs(R-ev).max())
# now with a realistic eigenvector error eps
for eps in (1e-8,1e-5):
    W=V+eps*np.random.randn(N,N); W/=np.linalg.norm(W,axis=0)
    R=np.array([ (W[:,k]@A@W[:,k]) for k in range(N)])
    d=R-ev
    print(f"  eps={eps:.0e}: R-mu  top={d[-1]:+.3e}  min={d.min():+.3e}  max={d.max():+.3e}"
          f"   (top must be <=0 always; interior sign is free)")

# ---- DEFECT 2: L1 off-by-one. Compare cert's loop to the baseline formula.
f=np.array([1.0/(1+k) for k in range(2*20+1)])
n=20
base = np.sum(2*f[1:n+1])/(2*n+1)      # KAY: sum_{k=1}^{n}
cert = np.sum(2*f[1:n])/(2*n+1)        # cert script: range(1,n)
print(f"\nL1_20  KAY/baseline={base:.10f}   cert={cert:.10f}   reldiff={abs(cert-base)/base:.3e}")
