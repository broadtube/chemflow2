"""サンプル: 前工程 pattern1（改質器）と後工程 plant3（MA プラント）をつなぐ。

plant3 の新鮮供給 Steam1 は pattern1 の DryGas そのものなので、改質条件を振れば
下流の挙動まで一気通貫で追える。**改質器へ戻る流れが無い一方向結合**なので、
逐次に解けば厳密（巨大な連立を一度に解く必要はなく、桁違いに速い）。

    [pattern1]  RG_feed(固定) + CO2_feed(未知) + H2O_feed(固定)
                    → Mixed → Gibbs(T,P) → 凝縮(25℃) → DryGas
    [plant3]    DryGas → M1 → R1(合成+脱水) → R2(カルボニル化) → … → 酢酸メチル

制約: **生成ガスの H2/CO は現状値 1.3837 を維持**。example_pattern1.py は Mixed の組成を固定して
Feed を逆算する形だったが、ここでは逆に **H2/CO を指定して CO2 供給量を逆算する**
（CO2 供給が未知、比の制約が 1 本増えて自由度は 20/20 で閉じる）。

背景: plant3 では **CO2 も CH4 も反応網の外にあり、パージ率の逆数（約20倍）まで
循環系に濃縮する**。既定条件では DryGas の CO2 16.09% が反応器入口
51.8 mol% まで濃縮し、CO を 12.7 mol% まで薄めていた。そこで DryGas 中の
**CO2 + CH4（= 不活性計）** を減らす条件を探す。

実行: PYTHONPATH=.:../reaction_rate/src python3 examples/example_pattern1_plant3.py
要 Cantera（改質器）と reaction_rate（速度論反応器）。
"""

import os
import sys

from scipy.optimize import brentq

import examples.example_plant3 as p3
from chemflow2.core.carbon_activity import carbon_activity_of
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
#
# 900 ℃・1.04 MPaG で H2/CO を維持したまま水蒸気供給を振ると（CH4 供給 2.1373 mol/h）:
#   H2O[mol/h]  H2O/CH4   CO2供給   不活性(CO2+CH4)    a_C
#      0.50      0.234    1.611      10.97 mol%     1.467  ✗ 析出域
#      0.75      0.351    1.827       9.82 mol%     1.040  ✗ 析出域
#      0.783     0.366    1.841       9.80 mol%     1.000  ← 炭素析出限界ちょうど
#      0.95      0.445    1.983       9.565 mol%    0.831  OK
#      1.00      0.468    2.020       9.561 mol%    0.789  OK  ← 採用（不活性の最小）
#      1.05      0.491    2.057       9.575 mol%    0.751  OK
#      2.78      1.299    3.163      14.47 mol%     0.242  OK（現状水準）
#
# 水蒸気を絞ると CO2 供給も減るので残存 CO2 は減るが、CH4 の未転化が増える。
# その和（不活性計）は **H2O ≈ 1.0 mol/h で極小**になり、そこでの a_C は 0.789。
# **つまりこの温度では炭素活量の制約は効いていない**（無制約の最適点が既に炭素フリー）。
# 制約が効き始めるのは H2O < 0.783 で、そこまで絞ると不活性はかえって増える。
CASES = {
    "A_base":     {"T": 850.0, "P": "1.04MPaG", "h2o": 2.776607},
    "B_900C_min": {"T": 900.0, "P": "1.04MPaG", "h2o": 1.0},
}

#: 出口ストリームの初期推定。Gibbs + 炭素活量の制約を同時に解くとき、既定の
#: 「全成分 1.0」からでは収束しないことがあるため、妥当な合成ガス組成を与える。
SEED_OUTLET = {"CO": 3.5, "H2": 4.9, "CO2": 1.7, "CH4": 0.25, "H2O": 1.0}

#: 検算用: A_base の DryGas はこの値（= plant3 の従来 Steam1）になるはず。
STEAM1_REFERENCE = {"H2": 4.911989, "CO": 3.549828, "CO2": 1.675074,
                    "CH4": 0.243186, "H2O": 0.028800}


def solve_reformer(*, T, P, h2o=None, h2_co=H2_CO_TARGET, carbon_limit=None,
                   verbose=True):
    """H2/CO を指定して CO2 供給量を逆算する形で pattern1 を解き、DryGas を返す。

    h2o を数値で渡すと水蒸気供給を固定する。自由度 20/20:
        変数   CO2_feed 1 + Mixed 4 + ReactOut 5 + DryGas 5 + Condensate 5
        方程式 Mixer 4 + Gibbs 5 + Separator 5 + 飽和 1 + Henry 4 + 比 1

    h2o=None・carbon_limit=1.0 とすると、**水蒸気供給も未知にして炭素活量が
    その値になる点（＝炭素析出限界）を直接解く**。変数・方程式とも 1 本増えて 21/21。
    a_C ≤ 1 は不等式なので等式系には載らないが、限界そのものを解けば
    「炭素析出しない最小の水蒸気量」が求まる。
    """
    if (h2o is None) == (carbon_limit is None):
        raise ValueError("h2o か carbon_limit のどちらか一方を指定してください")
    hot = StreamCondition(T=T, P=P, phase="gas")
    cold_g = StreamCondition(T=T_COND, P=P, phase="gas")
    cold_l = StreamCondition(T=T_COND, P=P, phase="liquid")

    RG   = Stream(["H2", "CH4"], name="1. RG_feed", order=1, condition=hot, flows=RG_FEED)
    CO2f = Stream(["CO2"],       name="2. CO2_feed", order=2, condition=hot)   # ← 未知
    H2Of = (Stream(["H2O"], name="3. H2O_feed", order=3, condition=hot, flows={"H2O": h2o})
            if h2o is not None else
            Stream(["H2O"], name="3. H2O_feed", order=3, condition=hot))       # ← 未知
    Mixed = Stream(["H2", "CO2", "CH4", "H2O"], name="4. Mixed", order=4, condition=hot)
    ReactOut   = Stream(SPECIES, name="5. ReactOut", order=5, condition=hot,
                        guess=SEED_OUTLET)
    DryGas     = Stream(SPECIES, name="6. DryGas", order=6, condition=cold_g,
                        guess=SEED_OUTLET)
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
    # 炭素析出限界を解く場合はさらに 1 本（水蒸気供給が未知になったぶん）
    if carbon_limit is not None:
        problem.constrain_carbon_activity(ReactOut, carbon_limit, T=T, P=P)

    # 初期推定を入れておくと Gibbs が素直に収束する
    CO2f.molar_flows[:] = 3.0
    if h2o is None:
        H2Of.molar_flows[:] = 2.0
    sol = problem.solve()
    if not sol.success:
        raise RuntimeError(f"改質器が収束せず: {sol}")
    if verbose:
        print(f"  改質器: T={T:.0f}℃ P={P} H2O={H2Of.flow_of('H2O'):.4f} → "
              f"CO2供給 {CO2f.flow_of('CO2'):.4f} mol/h  (自由度 {problem.degrees_of_freedom()})")
    return problem, DryGas, CO2f, H2Of


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
        _, DryGas, CO2f, _ = solve_reformer(T=cond["T"], P=cond["P"], h2o=cond["h2o"],
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

    # --- 炭素析出限界そのものを解く（a_C = 1 を課し、水蒸気供給を未知にする）---
    print("\n=== 炭素析出限界（a_C = 1）を直接解く ===")
    print(f"{'T[℃]':>6s}{'H2O[mol/h]':>11s}{'H2O/CH4':>9s}{'CO2供給':>9s}"
          f"{'不活性計':>9s}{'a_C':>9s}")
    for T in (850.0, 900.0):
        pr, dry, co2f, h2of = solve_reformer(T=T, P="1.04MPaG", carbon_limit=1.0,
                                             verbose=False)
        react_out = next(s for s in pr.streams if s.name == "5. ReactOut")
        s = drygas_summary(dry)
        h2o = h2of.flow_of("H2O")
        print(f"{T:6.0f}{h2o:11.4f}{h2o / RG_FEED['CH4']:9.4f}"
              f"{co2f.flow_of('CO2'):9.4f}{s['CO2%'] + s['CH4%']:9.3f}"
              f"{carbon_activity_of(react_out, T=T, P='1.04MPaG'):9.6f}")
    print("→ この限界より水蒸気を絞ると析出域に入る。ただし B の運転点（H2O=1.0）は\n"
          "  限界より水蒸気が多い側にあり、**炭素活量の制約は効いていない**。")

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
