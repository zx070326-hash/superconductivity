import numpy as np
rng = np.random.default_rng(20260918)

# Convention (A): Coulomb term written  sum_j mu_ij N_j Delta_j , mu SYMMETRIC.
#   kernel = mu @ diag(N)   -> symmetrise: D^{1/2} mu D^{1/2}, entries sqrt(N_i) mu_ij sqrt(N_j)
# Convention (B): Coulomb term written  sum_j mu_ij Delta_j , mu_ij = N_j K_ij, K SYMMETRIC.
#   kernel = K @ diag(N)    -> symmetrise: D^{1/2} K D^{1/2}, entries sqrt(N_i/N_j) mu_ij
#
# CLAIM to test: the criterion  mu_ii >= 0 and  mu_11 mu_22 >= mu_12 * mu_21
# gives the SAME verdict as PSD-of-the-correct-symmetrisation, under BOTH conventions.

def sym_A(mu, N):
    s = np.sqrt(N); return s[:,None]*mu*s[None,:]
def sym_B(mu, N):
    s = np.sqrt(N); return (s[:,None]/s[None,:])*mu

bad = 0
for trial in range(200000):
    B = 2
    N = rng.uniform(0.05, 5.0, B)
    if trial % 2 == 0:
        # build a convention-(A) matrix: mu symmetric, arbitrary sign
        M = rng.normal(0, 1, (B,B)); mu = (M+M.T)/2
        S = sym_A(mu, N)
    else:
        # build a convention-(B) matrix: K symmetric, mu_ij = N_j K_ij
        M = rng.normal(0, 1, (B,B)); K = (M+M.T)/2
        mu = K * N[None,:]
        S = sym_B(mu, N)
    psd = np.linalg.eigvalsh((S+S.T)/2).min() >= -1e-12
    # is S actually symmetric?
    assert abs(S[0,1]-S[1,0]) < 1e-12, (trial, S)
    crit = (mu[0,0] >= 0) and (mu[1,1] >= 0) and (mu[0,0]*mu[1,1] >= mu[0,1]*mu[1,0] - 1e-12)
    if psd != crit:
        bad += 1
        if bad <= 3: print("MISMATCH", trial, "psd", psd, "crit", crit, "\n mu=",mu,"\n N=",N)
print("2x2 both conventions: mismatches =", bad, "/ 200000")

# counter-check: the WRONG criterion mu_11 mu_22 >= mu_12**2 under convention (B)
bad2 = 0
for trial in range(200000):
    N = rng.uniform(0.05, 5.0, 2)
    M = rng.normal(0,1,(2,2)); K=(M+M.T)/2; mu = K*N[None,:]
    S = sym_B(mu,N); psd = np.linalg.eigvalsh((S+S.T)/2).min() >= -1e-12
    crit_wrong = (mu[0,0]>=0) and (mu[1,1]>=0) and (mu[0,0]*mu[1,1] >= mu[0,1]**2 - 1e-12)
    if psd != crit_wrong: bad2 += 1
print("convention (B) with mu_12**2 (WRONG form): mismatches =", bad2, "/ 200000")

# and: my old formula sym_A applied to a convention-(B) matrix -> not even symmetric
asym = []
for _ in range(5000):
    N = rng.uniform(0.05,5.0,2); M=rng.normal(0,1,(2,2)); K=(M+M.T)/2; mu=K*N[None,:]
    W = sym_A(mu,N); asym.append(abs(W[0,1]-W[1,0]))
print("sym_A on a convention-(B) matrix, max |W12-W21| =", max(asym))

# B>2: general statement is det/leading-minors, check PSD equivalence of K vs D^1/2 K D^1/2
worst = 0.0
for _ in range(20000):
    B = rng.integers(2,6)
    N = rng.uniform(0.05,5.0,B)
    M = rng.normal(0,1,(B,B)); K=(M+M.T)/2
    s=np.sqrt(N); S=s[:,None]*K*s[None,:]
    a = np.linalg.eigvalsh(K).min() >= 0
    b = np.linalg.eigvalsh(S).min() >= 0
    if a!=b: worst += 1
print("B=2..5 congruence preserves PSD verdict, disagreements =", int(worst), "/ 20000")
