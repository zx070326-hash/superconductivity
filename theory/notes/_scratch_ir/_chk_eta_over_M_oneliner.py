# 天花板 + 亏损预算 ==> 零自由参数的 Tc 公式  Tc[K] = K0 * sqrt(eta/M [meV^2])
# 严格性：C_infty 的陈述 analytic（g(2) 数值 rational-exact，12格/58有理量）；
#         f (填充因子) 是 numerical-evidence（4 行，float64，alpha2F 来自外部）。
import math
G2_LO, G2_HI = 1.3181093293757127, 1.3181615784372800   # CERT-B (rational-exact)
C_LO, C_HI = math.sqrt(G2_LO)/(2*math.pi), math.sqrt(G2_HI)/(2*math.pi)
MEV_K = 11.604518                                        # 1 meV in K
print(f"C_inf in [{C_LO:.10f}, {C_HI:.10f}]")

rows = [("H3S-200 anh", 38341.4, 203.0), ("H3S-200 har", 38347.4, 203.0),
        ("D3S-200 anh", 19830.0, 150.0), ("H3S-160 anh", 35659.3, 203.0)]
print(f"{'row':>12} {'eta/M':>9} {'ceil[K]':>8} {'Tc_exp':>7} {'f':>6} {'K0':>7}")
K0s = []
for name, em, tc in rows:
    ceil = C_HI*math.sqrt(em)*MEV_K          # 取 C_HI = 天花板上估，方向：上界要上估 ✓
    f = tc/ceil
    K0 = tc/math.sqrt(em)
    K0s.append(K0)
    print(f"{name:>12} {em:9.1f} {ceil:8.2f} {tc:7.1f} {f:6.3f} {K0:7.4f}")
K0m = sum(K0s)/len(K0s)
print(f"\nK0 range [{min(K0s):.4f}, {max(K0s):.4f}]  mean {K0m:.4f}  spread {100*(max(K0s)-min(K0s))/K0m:.1f}%")
print(f"ceiling constant (f=1): {C_HI*MEV_K:.4f}   => f_mean = {K0m/(C_HI*MEV_K):.4f}")
print("\n回代检验 Tc = K0_mean*sqrt(eta/M):")
for name, em, tc in rows:
    p = K0m*math.sqrt(em)
    print(f"  {name:>12}  pred {p:7.2f} K   exp {tc:6.1f}   resid {100*(p-tc)/tc:+6.2f}%")
