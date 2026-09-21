import numpy as np
rng = np.random.default_rng(11)

def O(f, N=140):
    n = np.arange(N); u = (2*n+1)**-0.5
    fv = np.array([f(k) for k in range(2*N+2)])
    L1 = np.diag(np.concatenate([[0.0], np.cumsum(2*fv[1:N])]) / (2*n+1))
    L2 = fv[np.abs(n[:,None]-n[None,:])] * np.outer(u,u); np.fill_diagonal(L2, 0.0)
    L3 = fv[n[:,None]+n[None,:]+1] * np.outer(u,u)
    return -L1 + L2 + L3

# sanity: g(2) on f(k)=<w^2>/k^2 with <w^2>=1  -> lam1 should approach 1.3181...
f2 = lambda k: 0.0 if k == 0 else 1.0/k**2
print("sanity  lam1(O[1/k^2]) at N=140 :", np.linalg.eigvalsh(O(f2)).max(),
      " (repo g(2)=1.318140212132, N=3200)")

# Claim: for a sign-INDEFINITE pairing kernel Lambda_ij, the (G') chain survives
# verbatim with M2 replaced ELEMENTWISE by |M2|.
#   (G'-1) N_ij := ||X_ij||  is >= 0 BY CONSTRUCTION, whatever the sign of Lambda_ij
#   O[-|f|] = -O[|f|]  =>  ||O[Lambda_ij]|| = ||O[|Lambda_ij|]||
#   (G'-3) then needs Lemma (G) on |Lambda_ij| only.
print("\n--- O is elementwise-linear in f, so O[-f] = -O[f] (norms equal) ---")
A = O(f2); B = O(lambda k: -f2(k))
print("   max|O[-f] + O[f]| =", np.abs(A+B).max(), "   ||.|| equal?",
      abs(np.linalg.norm(A,2)-np.linalg.norm(B,2)) < 1e-12)

print("\n--- full chain test with sign-indefinite blocks (s+- style) ---")
g2, kappa = 1.318140212132, 1.0
bad1 = bad3 = 0
for trial in range(400):
    B_ = rng.integers(2,5)
    w2 = rng.uniform(0.5, 8.0, (B_,B_)); w2 = (w2+w2.T)/2          # <w^2>_ij  > 0
    c  = rng.uniform(0.2, 2.0, (B_,B_)); c = (c+c.T)/2             # magnitudes
    sgn = rng.choice([-1.0,1.0], (B_,B_)); sgn = np.triu(sgn); sgn = sgn+np.triu(sgn,1).T
    Lam = c*sgn                                                    # SIGN-INDEFINITE
    M2  = Lam*w2                                                   # M2_ij = int 2w a2F_ij  (signed)
    Nn = 90
    blocks = [[O(lambda k,i=i,j=j: 0.0 if k==0 else Lam[i,j]*w2[i,j]/k**2, Nn)
               for j in range(B_)] for i in range(B_)]
    X = np.block(blocks)
    Nmat = np.zeros((B_,B_))
    for i in range(B_):
        for j in range(B_):
            Nmat[i,j] = np.linalg.norm(blocks[i][j], 2)            # >= 0 by construction
    if np.linalg.eigvalsh(X).max() > np.linalg.eigvalsh(Nmat).max() + 1e-9: bad1 += 1
    ceil = g2*max(abs(np.linalg.eigvals(np.abs(M2))))/kappa        # rho(|M2|)
    if np.linalg.eigvalsh(Nmat).max() > ceil + 1e-9: bad3 += 1
print("  (G'-1) lam1(X) <= lam_max(N)               violations:", bad1, "/400")
print("  (G'-2,3) lam_max(N) <= g(2)*rho(|M2|)/kap  violations:", bad3, "/400")

print("\n--- counter-check: using rho(M2) with SIGNS kept can UNDERSHOOT ---")
worse = 0
for _ in range(20000):
    B_=3; A_=rng.uniform(-3,3,(B_,B_)); M=(A_+A_.T)/2
    if max(abs(np.linalg.eigvals(M))) > max(abs(np.linalg.eigvals(np.abs(M))))+1e-12: worse+=1
print("  rho(M2) > rho(|M2|)  in", worse, "/20000   (<=0 expected: rho is monotone on |.|)")
