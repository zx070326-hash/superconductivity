"""
独立核验 team-lead 第四轮追加的判据（阻挫判据的精确划界）：

Dias & Marques, arXiv:1102.0496, "Frustrated multiband superconductivity" 用
「带间排斥总数的奇偶」判阻挫；本笔记 §3.3 用的是「符号图是否平衡」（每个环负边数为偶）。
两者不是一回事，只在 n=3（Dias-Marques 真正求解的情形）上等价。

本脚本从零独立写：对完全图 K_n（n=3,4,5），穷举全部 2^{C(n,2)} 种边符号赋值
（全部带间耦合都非零——不用零耦合，零耦合是 CLAUDE.md 说的原子型退化集，物理上不一般），
对每种赋值：
  (a) 直接判定「符号图平衡」= 是否存在顶点符号 s_i∈{±1} 使得对每条边 sign(i,j)=s_i*s_j
      （对完全图穷举 2^{n-1} 种顶点符号赋值即可，s_1 固定为 +1 不失一般性）；
  (b) 数负边条数 m；
  (c) 检验「m 奇 ⟹ 不平衡（阻挫）」是否对每个赋值成立，数反例；
  (d) 检验「m 偶 ⟹ 平衡（不阻挫）」是否对每个赋值成立，数反例。

再独立验证 analytic 恒等式：平衡态下负边数 = |S|(n-|S|)（S = 顶点符号为 -1 的集合），
从而：n 奇 ⟹ |S| 与 n-|S| 至少一个是偶数 ⟹ 平衡时负边数必偶 ⟹「m奇⟹阻挫」对 n 奇成立；
n 偶 ⟹ |S|=n/2 时 |S|(n-|S|) 可以是奇数（n/2 奇）⟹ 存在奇数负边的平衡态 ⟹ 该蕴含式失效。

Level: analytic（有限精确整数枚举，无浮点、无随机）。
"""
import itertools


def edges(n):
    return list(itertools.combinations(range(n), 2))


def is_balanced(n, edge_sign):
    """edge_sign: dict (i,j)->+1/-1 for i<j. 穷举顶点符号(s_1固定+1)看是否存在
    s_i 使 s_i*s_j == edge_sign[(i,j)] 对所有边成立。"""
    idx = list(range(1, n))
    for bits in itertools.product([1, -1], repeat=n - 1):
        s = [1] + list(bits)
        ok = True
        for (i, j) in edge_sign:
            if s[i] * s[j] != edge_sign[(i, j)]:
                ok = False
                break
        if ok:
            return True, s
    return False, None


def neg_count(edge_sign):
    return sum(1 for v in edge_sign.values() if v == -1)


def scan(n):
    E = edges(n)
    m_total = len(E)
    odd_implies_frustrated_counterex = []
    even_implies_balanced_counterex = []
    balanced_neg_counts = set()
    n_configs = 0
    for bits in itertools.product([1, -1], repeat=m_total):
        edge_sign = dict(zip(E, bits))
        n_configs += 1
        bal, s = is_balanced(n, edge_sign)
        m = neg_count(edge_sign)
        if bal:
            balanced_neg_counts.add(m)
        if (m % 2 == 1) and bal:
            # 负边数是奇数，但图仍平衡（不阻挫）——违反"奇 => 阻挫"
            odd_implies_frustrated_counterex.append((edge_sign, m, s))
        if (m % 2 == 0) and (not bal):
            # 负边数是偶数，但图不平衡（阻挫）——违反"偶 => 不阻挫"
            even_implies_balanced_counterex.append((edge_sign, m))
    return {
        "n": n,
        "n_edges": m_total,
        "n_configs": n_configs,
        "odd_counterex": odd_implies_frustrated_counterex,
        "even_counterex": even_implies_balanced_counterex,
        "balanced_neg_counts": sorted(balanced_neg_counts),
    }


def analytic_predicted_neg_counts(n):
    """平衡态下负边数 = |S|(n-|S|)，S 取遍 0..n。"""
    return sorted(set(s * (n - s) for s in range(0, n + 1)))


def main():
    for n in (3, 4, 5):
        res = scan(n)
        pred = analytic_predicted_neg_counts(n)
        print(f"=== K_{n}: 边数={res['n_edges']}，全部符号赋值={res['n_configs']} ===")
        print(f"  「负边数奇 ⟹ 阻挫(不平衡)」反例个数 = {len(res['odd_counterex'])}")
        print(f"  「负边数偶 ⟹ 不阻挫(平衡)」反例个数 = {len(res['even_counterex'])}")
        print(f"  平衡态下出现过的负边数集合（穷举实测） = {res['balanced_neg_counts']}")
        print(f"  |S|(n-|S|) 解析预测集合                = {pred}")
        assert res["balanced_neg_counts"] == pred, "穷举与解析公式 |S|(n-|S|) 不符"
        # 打印一个「负边数奇但平衡」的具体反例（若存在）
        if res["odd_counterex"]:
            es, m, s = res["odd_counterex"][0]
            print(f"  显式反例（奇数负边仍平衡）：负边={[e for e,v in es.items() if v==-1]}"
                  f"，负边数={m}，平衡顶点符号={s}")
        print()

    print("=== 解析划界（自己重推，两行）===")
    for n in (3, 4, 5, 6, 7):
        # n 奇 => |S| 与 n-|S| 至少一个偶 => |S|(n-|S|) 恒偶 => "奇=>阻挫" 应该对 n 奇成立
        # n 偶 => 取 |S|=n/2，若 n/2 为奇数，|S|(n-|S|)=(n/2)^2 也可能是奇数 (奇^2=奇)
        odd_n = (n % 2 == 1)
        exists_odd_balanced = any((s * (n - s)) % 2 == 1 for s in range(0, n + 1))
        print(f"  n={n} ({'奇' if odd_n else '偶'}): "
              f"存在奇数负边的平衡态？{exists_odd_balanced}  "
              f"⟹ 「负边数奇⟹阻挫」{'应成立' if not exists_odd_balanced else '应失效'}")
        if odd_n:
            assert not exists_odd_balanced, "n 奇时不应存在奇数负边的平衡态"
        # 对 n>=4，"偶=>不阻挫"总该失效（除非负边数集合里不含非平衡偶数，这里直接引用 scan 结果验证）

    print("\n=== 与 scan() 结果交叉核对（n=3,4,5 的划界预测应与穷举一致）===")
    for n in (3, 4, 5):
        res = scan(n)
        odd_ok = (len(res["odd_counterex"]) == 0)
        expect_odd_ok = (n % 2 == 1)
        print(f"  n={n}: 「奇⟹阻挫」穷举成立？{odd_ok}  解析预测应成立？{expect_odd_ok}")
        assert odd_ok == expect_odd_ok, f"n={n} 处解析划界与穷举不符"
        even_ok = (len(res["even_counterex"]) == 0)
        print(f"  n={n}: 「偶⟹不阻挫」穷举成立？{even_ok}（n=3 才可能侥幸成立，n>=4 应失效）")


if __name__ == "__main__":
    main()
