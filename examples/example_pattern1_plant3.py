"""サンプル: 前工程 pattern1（改質器）と後工程 plant3（MA プラント）をつなぐ。

plant3 の新鮮供給 Steam1 は pattern1 の DryGas そのものなので、改質条件を振れば
下流の挙動まで一気通貫で追える。**改質器へ戻る流れが無い一方向結合**なので、
逐次に解けば厳密（巨大な連立を一度に解く必要はなく、桁違いに速い）。

    [pattern1]  RG_feed(固定) + CO2_feed(未知) + H2O_feed(固定)
                    → Mixed → Gibbs(T,P) → 凝縮(25℃) → DryGas
    [plant3]    DryGas → M1 → R1(合成+脱水) → R2(カルボニル化) → … → 酢酸メチル

制約: **生成ガスの H2/CO = 1.3**。example_pattern1.py は Mixed の組成を固定して
Feed を逆算する形だったが、ここでは逆に **H2/CO を指定して CO2 供給量を逆算する**
（CO2 供給が未知、比の制約が 1 本増えて自由度は 20/20 で閉じる）。

背景: plant3 では **CO2 も CH4 も反応網の外にあり、パージ率の逆数（約20倍）まで
循環系に濃縮する**。既定条件（H2/CO=1.384）では DryGas の CO2 16.09% が反応器入口
51.8 mol% まで濃縮し、CO を 12.7 mol% まで薄めていた。そこで DryGas 中の
**CO2 + CH4（= 不活性計）** を減らす条件を探す。

実行: PYTHONPATH=.:../reaction_rate/src python3 examples/example_pattern1_plant3.py
要 Cantera（改質器）と reaction_rate（速度論反応器）。
"""

import os
import sys

from scipy.optimize import brentq

import examples.example_plant3 as p3
from chemflow2 import (
    GibbsReactor,
    Mixer,
    Problem,
    Separator,
    Stream,
    StreamCondition,
    stream_table,
)

# --- 改質器の固定条件 ---
SPECIES = ["CO2", "CH4", "H2O", "CO", "H2"]     # Gibbs 平衡種
DISSOLVED = ["H2", "CO", "CO2", "CH4"]          # 液水に溶けるガス
RG_FEED = {"H2": 0.885957, "CH4": 2.137257}     # 原料ガス（pattern1 の逆算値・固定）
T_COND = 25.0

#: 生成ガスの H2/CO。**現状値を維持する**（= 元の pattern1 の DryGas の比）。
#: MA 合成の総括量論 3CO + 4H2 → CH3COOCH3 + H2O が要求する 1.333 に近い。
H2_CO_TARGET = 4.911989 / 3.549828              # = 1.3837259

# --- 比較する改質条件（いずれも H2/CO は上の値に固定）---
# A は元の pattern1 の条件そのものなので、DryGas は Steam1 を再現するはず（検算になる）。
# B は制約下での最良点。制約は 2 つ:
#   (1) 改質温度は **900 ℃ 以下**
#   (2) 改質器出口の **炭素活量 a_C ≤ 1**（examples/example_carbon_limit.py 参照）
# 900 ℃・1.04 MPaG で H2/CO を維持したまま水蒸気を絞ると:
#     H2O/CH4  不活性(CO2+CH4)   a_C
#       0.50      10.97 mol%    1.467  ✗ 析出域
#       0.75       9.82 mol%    1.040  ✗ 析出域（わずかに超える）
#       1.00       9.56 mol%    0.789  OK  ← 採用
#       2.78      14.47 mol%    0.242  OK（現状水準）
# → 水蒸気を絞るほど不活性は減るが a_C が上がる。両立する下限が H2O/CH4 = 1.0。
# ⚠️ a_C は上表のとおり事前に確認した値であって、まだソルバの制約式にはなっていない。
CASES = {
    "A_base":       {"T": 850.0, "P": "1.04MPaG", "h2o": 2.776607},
    "B_900C_steam1": {"T": 900.0, "P": "1.04MPaG", "h2o": 1.0},
}

#: 検算用: A_base の DryGas はこの値（= plant3 の従来 Steam1）になるはず。
STEAM1_REFERENCE = {"H2": 4.911989, "CO": 3.549828, "CO2": 1.675074,
                    "CH4": 0.243186, "H2O": 0.028800}


def solve_reformer(*, T, P, h2o, h2_co=H2_CO_TARGET, verbose=True):
    """H2/CO を指定して CO2 供給量を逆算する形で pattern1 を解き、DryGas を返す。

    自由度 20/20:
        変数   CO2_feed 1 + Mixed 4 + ReactOut 5 + DryGas 5 + Condensate 5
        方程式 Mixer 4 + Gibbs 5 + Separator 5 + 飽和 1 + Henry 4 + 比 1
    """
    hot = StreamCondition(T=T, P=P, phase="gas")
    cold_g = StreamCondition(T=T_COND, P=P, phase="gas")
    cold_l = StreamCondition(T=T_COND, P=P, phase="liquid")

    RG   = Stream(["H2", "CH4"], name="1. RG_feed", order=1, condition=hot, flows=RG_FEED)
    CO2f = Stream(["CO2"],       name="2. CO2_feed", order=2, condition=hot)   # ← 未知
    H2Of = Stream(["H2O"],       name="3. H2O_feed", order=3, condition=hot,
                  flows={"H2O": h2o})
    Mixed = Stream(["H2", "CO2", "CH4", "H2O"], name="4. Mixed", order=4, condition=hot)
    ReactOut   = Stream(SPECIES, name="5. ReactOut", order=5, condition=hot)
    DryGas     = Stream(SPECIES, name="6. DryGas", order=6, condition=cold_g)
    Condensate = Stream(SPECIES, name="7. Condensate", order=7, condition=cold_l)

    problem = Problem(
        streams=[RG, CO2f, H2Of, Mixed, ReactOut, DryGas, Condensate],
        units=[Mixer([RG, CO2f, H2Of], Mixed, name="M0"),
               GibbsReactor(inlet=Mixed, outlet=ReactOut, species=SPECIES,
                            T=T, P=P, name="G1"),
               Separator(ReactOut, [DryGas, Condensate], name="Condenser")],
        name=f"Reformer {T:.0f}C {P} H2O={h2o}")
    problem.constrain_saturation(DryGas, "H2O", T=T_COND, P=P)
    problem.constrain_henry(DryGas, Condensate, DISSOLVED, T=T_COND, P=P)
    # 生成ガスの H2/CO を指定（残り 1 本）
    problem.constrain(DryGas.flow_expr("H2"), h2_co * DryGas.flow_expr("CO"),
                      name=f"H2/CO={h2_co}")

    # CO2 供給の初期推定を入れておくと Gibbs が素直に収束する
    CO2f.molar_flows[:] = 3.0
    sol = problem.solve()
    if not sol.success:
        raise RuntimeError(f"改質器が収束せず: {sol}")
    if verbose:
        print(f"  改質器: T={T:.0f}℃ P={P} H2O={h2o:.4f} → "
              f"CO2供給 {CO2f.flow_of('CO2'):.4f} mol/h  (自由度 {problem.degrees_of_freedom()})")
    return problem, DryGas, CO2f


def drygas_summary(DryGas) -> dict:
    tot = float(DryGas.total_flow.eval())
    return {
        "total": tot,
        "H2/CO": DryGas.flow_of("H2") / DryGas.flow_of("CO"),
        "CO": DryGas.flow_of("CO"),
        "CO2%": DryGas.flow_of("CO2") / tot * 100,
        "CH4%": DryGas.flow_of("CH4") / tot * 100,
        "CO%": DryGas.flow_of("CO") / tot * 100,
    }


def as_plant_feed(DryGas) -> dict[str, float]:
    """DryGas（5成分）を plant3 の新鮮供給の形（無い成分は 0）に直す。"""
    return {f: DryGas.flow_of(f) for f in DryGas.formulas}


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None

    # --- 改質側だけ先に比較（Gibbs のみなので一瞬）---
    print(f"=== 前工程 pattern1: H2/CO={H2_CO_TARGET} 制約下での改質条件比較 ===")
    print(f"{'case':>10s} {'T[℃]':>6s} {'H2O':>7s} {'CO2feed':>8s} {'DryGas':>8s} "
          f"{'H2/CO':>6s} {'CO2%':>7s} {'CH4%':>7s} {'不活性計':>8s} {'CO%':>7s} {'CO':>7s}")
    feeds, summaries = {}, {}
    for tag, cond in CASES.items():
        _, DryGas, CO2f = solve_reformer(T=cond["T"], P=cond["P"], h2o=cond["h2o"],
                                         verbose=False)
        s = drygas_summary(DryGas)
        feeds[tag], summaries[tag] = as_plant_feed(DryGas), s
        print(f"{tag:>10s} {cond['T']:6.0f} {cond['h2o']:7.4f} {CO2f.flow_of('CO2'):8.4f} "
              f"{s['total']:8.4f} {s['H2/CO']:6.4f} {s['CO2%']:7.2f} {s['CH4%']:7.2f} "
              f"{s['CO2%']+s['CH4%']:8.2f} {s['CO%']:7.2f} {s['CO']:7.4f}")

    # A_base は元の pattern1 条件そのものなので、DryGas は Steam1 と一致するはず
    worst = max(abs(feeds["A_base"][f] - v) for f, v in STEAM1_REFERENCE.items())
    print(f"\n検算: A_base の DryGas と plant3 の従来 Steam1 の最大差 = {worst:.2e} mol/h "
          f"→ {'一致' if worst < 1e-5 else '不一致'}")

    # --- 各 DryGas を plant3 に流す（速度論なので 1 ケース 10 分前後）---
    results = {}
    for tag, feed in feeds.items():
        if only and tag != only:
            continue
        print(f"\n=== 後工程 plant3: 供給 = {tag} の DryGas ===")
        problem, streams, v_tot = p3.solve_with_sv(feed=feed)
        (Steam1, ReactorIn, HybridOut, ReactorOut, CondGas, CondLiq, Purge,
         RecycleGas, MA, Sep1Liq, RecycleMeOH, Water) = streams
        n_in = float(ReactorIn.total_flow.eval())
        results[tag] = {
            "V_cat": v_tot * 1e6,
            "ReactorIn": n_in,
            "CO2_in%": ReactorIn.flow_of("CO2") / n_in * 100,
            "CH4_in%": ReactorIn.flow_of("CH4") / n_in * 100,
            "CO_in%": ReactorIn.flow_of("CO") / n_in * 100,
            "MA": MA.flow_of("CH3COOCH3"),
            "CO_purge": Purge.flow_of("CO"),
            "CO_feed": Steam1.flow_of("CO"),
        }
        print(stream_table(problem.streams, basis=["mol"]))

        out = os.path.join(os.path.dirname(__file__), "output")
        os.makedirs(out, exist_ok=True)
        from chemflow2 import to_csv
        to_csv(problem.streams, os.path.join(out, f"plant3_{tag}_streams.csv"),
               basis=["mol", "mole_frac"])

    # --- まとめ ---
    if len(results) > 1:
        print("\n=== 前工程条件が後工程に効く度合い ===")
        print(f"{'指標':>26s}" + "".join(f"{t:>14s}" for t in results))
        rows = [("DryGas 中 CO2 [mol%]", lambda t: summaries[t]["CO2%"]),
                ("DryGas 中 CH4 [mol%]", lambda t: summaries[t]["CH4%"]),
                ("反応器入口 CO2 [mol%]", lambda t: results[t]["CO2_in%"]),
                ("反応器入口 CH4 [mol%]", lambda t: results[t]["CH4_in%"]),
                ("反応器入口 CO [mol%]", lambda t: results[t]["CO_in%"]),
                ("反応器入口 [mol/h]", lambda t: results[t]["ReactorIn"]),
                ("全触媒体積 [mL]", lambda t: results[t]["V_cat"]),
                ("酢酸メチル [mol/h]", lambda t: results[t]["MA"]),
                ("CO パージ損失 [mol/h]", lambda t: results[t]["CO_purge"]),
                ("新鮮CO基準 炭素収率 [%]",
                 lambda t: results[t]["MA"] * 3 / results[t]["CO_feed"] * 100)]
        for label, fn in rows:
            print(f"{label:>26s}" + "".join(f"{fn(t):14.4f}" for t in results))


if __name__ == "__main__":
    main()
