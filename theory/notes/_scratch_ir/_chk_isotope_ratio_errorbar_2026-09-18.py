import random, math
for _ in range(5):
    fp_H, fp_D, fo_H, fo_D = [random.uniform(0.1, 2) for _ in range(4)]
    r_H = fp_H / fo_H - 1
    r_D = fp_D / fo_D - 1
    lhs = (1 + r_H) / (1 + r_D)
    ratio_pred = fp_D / fp_H
    ratio_obs = fo_D / fo_H
    rhs = ratio_obs / ratio_pred
    print(lhs, rhs, lhs - rhs)

# error propagation with independently-extracted D3S double-marker scatter
D3S_pairs = [(200.00, 145.00, 145.99), (216.98, 140.04, 142.02),
             (242.00, 132.04, 137.99), (256.00, 126.03, 130.00)]
rel_diffs = []
for P, a, b in D3S_pairs:
    m = 0.5 * (a + b)
    rel_diffs.append(abs(a - b) / m)
    print(f"P={P:6.1f}  {a:.2f}/{b:.2f}  rel_diff={100*rel_diffs[-1]:.2f}%")
mean_rd = sum(rel_diffs) / len(rel_diffs)
median_rd = sorted(rel_diffs)[1:3]
median_rd = sum(median_rd) / 2
print(f"mean rel diff = {100*mean_rd:.2f}%   median = {100*median_rd:.2f}%")

# H3S duplicate at 182.97 GPa: 175.99 / 185.03
h3s_pair_reldiff = abs(185.03 - 175.99) / (0.5*(185.03+175.99))
print(f"H3S duplicate @182.97 GPa: 175.99/185.03  rel_diff={100*h3s_pair_reldiff:.2f}%")

# per-point sigma from pair differences: Var(diff)=2*sigma^2 -> sigma = diff/sqrt(2)
sigma_D3S = [rd/math.sqrt(2) for rd in rel_diffs]
sigma_H3S = h3s_pair_reldiff/math.sqrt(2)
print("sigma per D3S point (approx):", [f"{100*s:.2f}%" for s in sigma_D3S])
print(f"sigma per H3S point (from its own duplicate): {100*sigma_H3S:.2f}%")

# representative single-point sigma (use median of D3S sigmas as the working estimate)
sigma_rep = sorted(sigma_D3S)[1:3]
sigma_rep = sum(sigma_rep)/2
print(f"representative single-point sigma (median of D3S) = {100*sigma_rep:.2f}%")

# f_obs(D3S) is the MEAN of two points at 200 GPa -> sigma_mean = sigma_rep/sqrt(2)
sigma_D_mean = sigma_rep/math.sqrt(2)
# f_obs(H3S) is a SINGLE point (no duplicate at 199.97 GPa) -> use sigma_rep directly (no averaging)
sigma_H_single = sigma_rep
sigma_ratio = math.sqrt(sigma_D_mean**2 + sigma_H_single**2)
print(f"propagated sigma on ratio_obs (D mean vs H single) = {100*sigma_ratio:.2f}%")

# conservative alternative: both single-point (no averaging benefit assumed)
sigma_ratio_conservative = math.sqrt(sigma_rep**2 + sigma_rep**2)
print(f"conservative sigma (both single-point, no averaging credit) = {100*sigma_ratio_conservative:.2f}%")

isotope_resid = -1.24
print(f"isotope-ratio residual = {isotope_resid}%  vs propagated sigma range [{100*sigma_ratio:.2f}%, {100*sigma_ratio_conservative:.2f}%]")
print(f"|residual|/sigma in [{abs(isotope_resid)/(100*sigma_ratio_conservative):.2f}, {abs(isotope_resid)/(100*sigma_ratio):.2f}]")

# --- ROUND 3 (team-lead's H3S-dominant convention, independently reproduced) ---
# SD_D3S combined via RMS (quadrature mean) of the four per-pair SDs, not arithmetic mean/median.
sd_D3S_rms = math.sqrt(sum(s**2 for s in sigma_D3S) / len(sigma_D3S))
SE_D3S = sd_D3S_rms / math.sqrt(2)          # D3S@200 anchor is a 2-point average
SE_H3S = sigma_H3S                           # H3S@200 anchor is a single point, no /sqrt(2) credit
sigma_ratio_v3 = math.sqrt(SE_H3S**2 + SE_D3S**2)
print(f"\n[round3] SD_D3S (RMS of 4 pairs) = {100*sd_D3S_rms:.2f}%   SE_D3S = {100*SE_D3S:.2f}%")
print(f"[round3] SD_H3S (single pair)     = {100*SE_H3S:.2f}%   SE_H3S = {100*SE_H3S:.2f}%")
print(f"[round3] sigma(ratio_obs) = sqrt(SE_H3S^2+SE_D3S^2) = {100*sigma_ratio_v3:.2f}%   (H3S-dominated)")
print(f"[round3] |residual|/sigma = {abs(isotope_resid)/(100*sigma_ratio_v3):.2f}")

# common-offset significance: offset is a geometric mean -> log-space sigma = sigma_ratio/2
sigma_log = sigma_ratio_v3 / 2
common_offset = 7.44
print(f"[round3] sigma_log (common-offset) = {100*sigma_log:.2f}%   "
      f"+{common_offset}% / sigma_log = {common_offset/(100*sigma_log):.2f} sigma")

