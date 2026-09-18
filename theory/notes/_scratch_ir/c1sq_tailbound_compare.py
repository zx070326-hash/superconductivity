import numpy as np

def f(k, w): return w*w/(w*w+k*k)

def tail_true(J, w, JMAX=2000000):
    j = np.arange(J, JMAX)
    return f(j, w).sum()

def tail_arctanfree(J, w):
    y = (J-1)/w
    return w*(np.pi/2 - y/(1+y*y))

def tail_rat(J, w):
    if J <= 1:
        return 2.0*w*w
    return w*w/(J-1)

for w in (0.5, 1, 5, 30, 100):
    for J in (1, 2, 10, int(w)+1, 2*int(w)+5, 501):
        t = tail_true(J, w)
        af = tail_arctanfree(J, w)
        rt = tail_rat(J, w)
        print("w=%6.1f J=%5d true=%12.6f arctanfree=%12.6f (x%6.2f) rat=%12.6f (x%8.2f)" %
              (w, J, t, af, af/t, rt, rt/max(t, 1e-12)))
