# -*- coding: utf-8 -*-
"""任务 ① 收尾：KT 里唯一剩下的 fail-open 件——谱隙 γ 的认证代价。

方向审计：`gamma = λ₂(P_N 𝔎_ph P_N)` 是真 λ₂ 的**下估** ⟹ μ−γ 偏大 ⟹ KT 偏紧
⟹ 这一处是 fail-open，必须换成 λ₂ 的**上界**。问题：那个上界要多准？
本脚本把 γ 朝**判据变难**的方向退化（γ → γ + s(μ−γ)，s∈[0,1)，即假装 λ₂ 被高估到
离 μ 只剩 (1−s) 的间距），扫 δ=KT/λ₁−1 在哪里开始吃掉预算（消费点：δ=1e−3 ⟹ λ_ceil −0.26%）。
方向：退化使 μ−γ 变小 ⟹ KT 变大 ⟹ 若退化后仍 δ≪1e−3，则任何**更好**的 λ₂ 上界也够
（单边界方向规则的标准应用）。  numerical-evidence。
"""
import os; os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util, numpy as np
p = "theory/notes/_scratch_ir/_chk_main_kato_temple_two_atom_prototype_2026-09-23.py"
s_ = importlib.util.spec_from_file_location("kt", p); m_ = importlib.util.module_from_spec(s_)
s_.loader.exec_module(m_)

r, Nref = 10.0, 1500
print("=" * 100)
print("GATE G  γ 退化敏感度： γ_used = γ + s(μ−γ)  （s=0 即点算值；s→1 即谱隙被吃光）")
print(f"  {'case':<24}{'μ':>12}{'γ(pt)':>12}{'μ−γ':>10}" + "".join(f"{f'δ(s={s})':>13}" for s in (0.0,0.5,0.9,0.99)))
for lam, mus in ((2.0,0.13),(3.2953,0.13),(6.0,0.20)):
    wE = m_.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
    for tag,(x,rho) in (("Einstein",(0.0,1.0)),("2atom x=.4 ρ=3",(0.4,3.0)),("2atom x=.1 ρ=20",(0.1,20.0))):
        L_,W_ = (np.array([lam]),np.array([wE])) if tag=="Einstein" else m_.members(lam,wE,x,rho)
        ref = m_.lam_k(L_,W_,mus,Nc,Nref); g0 = m_.lam_k(L_,W_,0.0,0,Nref,k=2)
        _,mu,_,_ = m_.kato_temple(L_,W_,mus,Nc,120,g0)
        row = f"  {tag+f' λ={lam:.3g}':<24}{mu:>12.7f}{g0:>12.7f}{mu-g0:>10.5f}"
        for s in (0.0,0.5,0.9,0.99):
            kt,_,_,_ = m_.kato_temple(L_,W_,mus,Nc,120, g0 + s*(mu-g0))
            row += f"{kt/ref-1.0:>13.2e}"
        print(row)
print()
print("  读法：s=0.9 意味着 λ₂ 的上界只要把谱隙估准到 10% 就够；若该列仍 ≪1e−3，")
print("        则 γ 的认证是**宽松活**，不是瓶颈。")
