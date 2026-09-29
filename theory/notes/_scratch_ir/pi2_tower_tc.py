from mpmath import mp, mpf
mp.dps=50
b1 = mpf('0.3584070531484623077978476')
kap2 = mpf('0.09873042407289524101111089465604911963035')
C  = mpf('1.463510026021423479440976')   # gamma+ln4-1/2
c0 = C-1                                   # gamma+ln4-3/2
b2 = 2*kap2 - b1*C
tc2_direct  = b1*c0 + b1**2/2 + b2         # w2's physical Tc lambda^2 coeff
tc2_reduced = 2*kap2 - b1 + b1**2/2        # claimed identity (gamma,ln4 cancel)
print("b2            =", b2)
print("Tc2 (direct)  =", tc2_direct)
print("Tc2 (reduced) =", tc2_reduced)
print("identity diff =", tc2_direct - tc2_reduced)
print("w2 -0.0967184 diff =", tc2_reduced - mpf('-0.0967184'))
