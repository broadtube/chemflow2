r"""検証: 入口ガスの炭素析出限界温度 — 有識者の値（465 ℃ / 525 ℃）と突き合わせる。

`references/carbon_activity.html` §8.2 に記録した伝聞の検証。

────────────────────────────────────────────────────────────────────────
検証の対象
────────────────────────────────────────────────────────────────────────
実務者から次の 2 点が伝えられている（2026-08-28 記録・出典は口頭）:

    改質器 A の条件      入口ガスの a_C が 約 465 ℃ で 1 → 入口温度は 450 ℃ 程度
    旧 chemflow pattern1  同 約 525 ℃ で 1（入口 H2 9.7 vol% / CO2 36.5 vol%）→ 同 450 ℃

pattern1 の組成は `example_carbon_limit.py` の運転点と**完全に一致**する
（CH4 23.40 / CO2 36.50 / H2O 30.40 / H2 9.70 vol%）ので、ケースの同定は済んでいる。

────────────────────────────────────────────────────────────────────────
分かったこと
────────────────────────────────────────────────────────────────────────
1. **入口では CH4 分解経路しか効かない。** 入口ガスに CO が無いので Boudouard と
   CO 還元は厳密に 0。千代田の特許の定義（Boudouard 1 本）だと入口の a_C は 0 に
   なってしまうので、**有識者の計算は CH4 分解経路を使っている**ことが確定する。

2. **黒鉛基準だと 2 ケースとも 27 ℃ 低く出る**（438.0 対 465、497.5 対 525）。
   この 2 ケースは p_CH4/p_H2² が 6.757 対 2.208 と **3 倍違う**のに、ずれが揃う。

3. **炭素の基準を ΔG_dev ≈ 3.0〜3.2 kJ/mol ずらすと両方が同時に合う。**
   これは §5.3 の**ウィスカーカーボン補正**そのもののサイズ（Alstrup の式で
   d ≈ 200〜300 nm 相当）。補正量の一貫性で比べると:

       ΔG_dev 一定   3243 / 3092 J/mol   差 5%
       ΔT 一定       +27.0 / +27.5 ℃     差 2%
       K_M 倍率一定  0.5896 / 0.6275     差 6%

   **K_M を一定倍する形（＝平衡定数の出典違い）よりも、炭素基準をずらす形の方が
   一貫する。** ただし 2 点では ΔG 一定と ΔT 一定を判別できない（どちらも ±4 ℃）。

⚠️ **仮説の域を出ない。** 有識者側がウィスカー基準を使っているという確認は取れて
   いない。千代田の特許 5 件を全文検索しても「ウィスカー」「黒鉛」「グラファイト」
   「粒子径」は **0 件**で、炭素の基準についての記載自体が無い。
   判別には**温度の離れた 3 点目**（例えば 600 ℃ 付近）が要る。

実行: PYTHONPATH=. python3 -u examples/example_inlet_carbon_limit.py
要 Cantera: pip install chemflow2[gibbs]
"""

import math

from scipy.optimize import brentq

from chemflow2.core.carbon_activity import (
    R_GAS,
    carbon_activity,
    equilibrium_constants,
    whisker_dG_dev,
)
from chemflow2.core.pressure import parse_pressure

NV = 22.414   # ノルマルモル体積 [L/mol]

#: (名前, 供給 [mol/h], 圧力, 有識者が言う a_C = 1 の温度 [℃])
CASES = {
    "改質器 A": (dict(CH4=23.20 / NV, CO2=52.00 / NV, H2O=32.70 / NV, H2=5.10 / NV),
                "1.41MPaG", 465.0),
    "pattern1": (dict(CH4=2.137257, CO2=3.333756, H2O=2.776607, H2=0.885957),
                 "1.04MPaG", 525.0),
}


def inlet_partials(feed, P):
    """入口ガスの分圧 [atm]。CO は 0（原料に CO が無い）。"""
    total = sum(feed.values())
    P_atm = parse_pressure(P) / 101325.0
    p = {k: v / total * P_atm for k, v in feed.items()}
    p["CO"] = 0.0
    return p


def crossing(p_atm, dG_dev=0.0):
    """a_C = 1 になる温度 [℃]。"""
    return brentq(lambda T: carbon_activity(p_atm, T, dG_dev=dG_dev) - 1.0,
                  200.0, 900.0, xtol=1e-4)


def main():
    print("=== 入口ガスの炭素析出限界温度（a_C = 1）===\n")
    prepared = {}
    for name, (feed, P, expert) in CASES.items():
        p = inlet_partials(feed, P)
        total = sum(feed.values())
        x = {k: v / total * 100 for k, v in feed.items()}
        prepared[name] = (p, P, expert)
        print(f"--- {name} / {P} ---")
        print("  入口組成 [vol%]: " + "  ".join(f"{k} {x[k]:5.2f}"
                                             for k in ("CH4", "CO2", "H2O", "H2")))
        print(f"  p_CH4 = {p['CH4']:.3f} atm、p_H2 = {p['H2']:.4f} atm"
              f"  → p_CH4/p_H2² = {p['CH4'] / p['H2'] ** 2:.3f}")
        r = carbon_activity(p, expert, per_route=True)
        print(f"  {expert:.0f} ℃ での経路別 a_C: "
              f"CH4分解 {r['CH4_cracking']:.4f} / Boudouard {r['Boudouard']:.1f}"
              f" / CO還元 {r['CO_reduction']:.1f}")
        print("  → 入口は CO が無いので **CH4 分解しか効かない**\n")

    print(f"{'炭素の基準':34s}{'ΔG_dev':>9s}"
          + "".join(f"{n:>12s}" for n in CASES) + f"{'A の差':>9s}{'p1 の差':>9s}")
    rows = [("黒鉛（既定）", 0.0),
            ("ウィスカー 黒鉛極限 d→∞", whisker_dG_dev()),
            ("ウィスカー d = 300 nm", whisker_dG_dev(300)),
            ("ウィスカー d = 200 nm", whisker_dG_dev(200)),
            ("ウィスカー d = 100 nm", whisker_dG_dev(100))]
    for label, dG in rows:
        ts = {n: crossing(prepared[n][0], dG) for n in CASES}
        dA = ts["改質器 A"] - CASES["改質器 A"][2]
        dP = ts["pattern1"] - CASES["pattern1"][2]
        print(f"{label:34s}{dG:9.0f}"
              + "".join(f"{ts[n]:12.1f}" for n in CASES)
              + f"{dA:+9.1f}{dP:+9.1f}")

    print("\n--- 有識者の値に合わせるのに必要な補正量（3 つの形で比較）---")
    print(f"{'ケース':12s}{'有識者[℃]':>10s}{'本実装[℃]':>11s}"
          f"{'ΔG_dev[J/mol]':>15s}{'K_M 倍率':>11s}{'ΔT[℃]':>9s}")
    forms = {}
    for name, (p, P, expert) in prepared.items():
        a = carbon_activity(p, expert)
        dG = R_GAS * (expert + 273.15) * math.log(a)
        base = crossing(p)
        forms[name] = (dG, 1.0 / a, expert - base)
        print(f"{name:12s}{expert:10.0f}{base:11.1f}{dG:15.0f}{1 / a:11.4f}"
              f"{expert - base:+9.1f}")
    keys = list(forms)
    for i, label in enumerate(("ΔG_dev 一定", "K_M 倍率一定", "ΔT 一定")):
        v = [forms[k][i] for k in keys]
        print(f"  {label:14s} {v[0]:10.4f} 対 {v[1]:10.4f}"
              f"   → 差 {abs(v[0] - v[1]) / max(abs(v[0]), abs(v[1])) * 100:.1f}%")
    print("  → **炭素基準をずらす形（ΔG_dev）と温度シフトが一貫し、K_M 倍率は劣る**。")
    print("     ΔG_dev ≈ 3.0〜3.2 kJ/mol は Alstrup の式で d ≈ 200〜300 nm 相当。")

    print("\n--- 判別に必要なもの: 限界温度が離れた 3 点目 ---")
    print("  2 点では ΔG_dev 一定と ΔT 一定を判別できない。両者が分かれる条件を探す。")
    print("  入口の a_C は a_C = K_M(T)·Q（Q = p_CH4/p_H2²）なので、**限界温度は Q だけで決まる**。")
    print("  黒鉛基準の限界温度 T0 を横に振って、2 つの仮説の予測を並べる:")
    print(f"  {'T0（黒鉛基準）':>14s}{'Q = p_CH4/p_H2²':>18s}"
          f"{'ΔG_dev=3053':>14s}{'ΔT=+27.2':>12s}{'差':>8s}")
    for T0 in (300.0, 400.0, 465.0, 500.0, 600.0, 700.0, 800.0):
        Q = 1.0 / equilibrium_constants(T0)[1]
        p = {"CH4": Q, "H2": 1.0, "CO": 0.0, "CO2": 1.0, "H2O": 1.0}
        t_dg = crossing(p, 3053.0)
        t_dt = T0 + 27.2
        print(f"  {T0:14.0f}{Q:18.4f}{t_dg:14.1f}{t_dt:12.1f}{t_dg - t_dt:+8.1f}")
    print("  → 500 ℃ 前後では両仮説の差が 1 ℃ 未満で判別できない。")
    print("     **限界温度が 300 ℃ 台か 700 ℃ 以上になる条件**の値を 1 つもらえれば決着する。")
    print("     もっと確実なのは、有識者に**炭素の基準（黒鉛かウィスカーか）を直接聞く**こと。")


if __name__ == "__main__":
    main()
