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

背景: plant3 では CO2 も CH4 もパージ率の逆数（約20倍）まで循環系に濃縮する。
既定条件では DryGas の CO2 16.09% が反応器入口 51.8 mol% まで濃縮し、CO を
12.7 mol% まで薄めていた。そこで DryGas 中の **CO2 + CH4（= 循環蓄積成分）** を
減らす条件を探した。

────────────────────────────────────────────────────────────────────────
✅ 結論（2026-08-20 改訂）: **不活性最小化は正しかった。当初「裏目」としたのは
   指標の誤りだった。C/D 条件のほうが炭素利用率で大きく優る。**
────────────────────────────────────────────────────────────────────────

  指標                    A_base    B_900C  C_03MPaG  D_01MPaG
  改質器 CO2 供給 [mol/h]   3.3338    2.3980    1.8670    1.8458   ← **投入炭素が違う**
  改質 CH4 転化率 [%]      88.62     88.62     91.19     95.24
  改質 CO2 転化率 [%]      49.75     69.08     90.69     94.91   ← A は半分捨てている
  改質器出口 a_C            0.394     0.495     0.900     0.900
  DryGas 不活性計 [mol%]   18.43     10.39      3.97      2.09
  反応器入口 CO2 [mol%]    51.80     35.51     21.03     19.39
  反応器入口 [mol/h]       68.06     52.76     47.56     48.38
  全触媒体積 [mL]         305.11    236.52    213.23    216.91   ← −30%
  酢酸メチル [mol/h]       1.0416    1.0036    0.9653    0.9966   ← 微減
  新鮮CO基準 炭素収率 [%]   88.03     84.81     79.51     78.94   ← ❌ 不適切な指標
  **投入炭素基準 収率 [%]   57.12     66.39     72.32     75.06**  ← ✅ **D が最良**

**当初の誤り 1 — 指標:** 「新鮮 CO 基準」は分母が DryGas 中の CO だけで、**改質器に
投入した CO2 と CH4 を勘定に入れていない**。A_base は CO2 を 3.334 mol/h も食って
MA 1.04 を作っているのに、その消費が指標に現れなかった。分母を投入炭素
（CH4 2.137 + 改質器 CO2 供給）に取り直すと順位が完全に逆転し、**D_01MPaG が最良、
A_base が最下位**になる。改質器の段階で既に CO2 転化率 49.75% 対 94.91% と差がある。

**当初の誤り 2 — 機構の説明:** 「CO2 は RWGS を通じて CO を供給する」と書いたが、
物質収支と合っていなかった。**CO2 は全ケースで循環系の中で正味"生成"側**
（供給 < パージ）で、動いているのは順方向の水性ガスシフト CO + H2O → CO2 + H2、
つまり **CO を食う**反応。ループ内の CO2 は WGS の生成物側阻害として CO を守って
いた、というのが正しい読み。詳細は example_co2_removal.py の炭素収支を参照。

**それでも MA が微減するのは事実**（1.0416 → 0.9653〜0.9966）。ただし触媒量 −30%、
投入炭素 −27%（5.471 → 3.983 mol-C/h）と引き換えなので、**総合では C/D が優る**。

**次の一手（実施済み）:** 供給 CO2 を絞るのとは別に、**循環ガスから CO2 を除去して
改質器へ戻す**構成を評価した（example_co2_removal.py）。改質器の CO2 供給は
H2/CO 制約で総量が決まるので、回収 CO2 を戻しても改質器の解は不変。その前提で
**D_01MPaG + 循環除去 η=0.9 が 86.43%** と全構成で最良になった。
ただし NaOH は不可逆なので、回収 CO2 を戻すにはアミンや PSA が要る。

実行: PYTHONPATH=.:../reaction_rate/src python3 examples/example_pattern1_plant3.py
要 Cantera（改質器）と reaction_rate（速度論反応器）。1 ケース約 10 分。
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
    to_excel,
)

# --- 改質器の固定条件 ---
SPECIES = ["CO2", "CH4", "H2O", "CO", "H2"]     # Gibbs 平衡種
DISSOLVED = ["H2", "CO", "CO2", "CH4"]          # 液水に溶けるガス
RG_FEED = {"H2": 0.885957, "CH4": 2.137257}     # 原料ガス（pattern1 の逆算値・固定）
T_COND = 25.0

#: 生成ガスの H2/CO。**現状値を維持する**（= 元の pattern1 の DryGas の比）。
#: MA 合成の総括量論 3CO + 4H2 → CH3COOCH3 + H2O が要求する 1.333 に近い。
H2_CO_TARGET = 4.911989 / 3.549828              # = 1.3837259

#: CH4 転化率の下限＝**現状（元の pattern1・850 ℃）の値**。これ以上を維持する。
#: A_base（現状条件）を解くと 0.886207988601 になる。この値自体が基準なので
#: A_base は下限ちょうどに乗る。丸め差で「制約違反」と誤判定しないよう、
#: 判定には下の許容差を使う（1e-6 = 0.0001 %ポイント。物理的に無意味な差）。
CH4_CONVERSION_MIN = 0.886207988601
CH4_CONVERSION_TOL = 1e-6

#: 炭素活量の設計値。千代田の特許 WO2012140994A1 表 5 では**実施例 1〜4 のすべてで
#: 出口カーボン活性を 0.90 に設計**しているので、それに倣う（a_C ≤ 1 を余裕をもって満たす）。
CARBON_ACTIVITY_DESIGN = 0.90

# --- 比較する改質条件 ---
# 制約は 4 つ:
#   (1) 生成ガス H2/CO = 1.3837（現状維持）      → CO2 供給が決まる
#   (2) CH4 転化率 ≥ 88.6208%（現状以上）
#   (3) 改質温度 ≤ 900 ℃
#   (4) 改質器出口の炭素活量 a_C ≤ 1
# 目的は DryGas 中の **CO2 + CH4（plant3 で循環蓄積する成分）を最小化**すること。
# 圧力は設計変数として動かしてよい（ただし下流 plant3 は 5 MPaG なので昇圧コストは増える）。
#
# 転化率を上げるほど水蒸気が要り、H2/CO を戻すため CO2 も要るので不活性は増える。
# よって最適は転化率＝下限ちょうど。ただし**低圧側では下限が自動的に満たされる**ようになり、
# そこからは a_C が下限を決める。その場合は千代田の実務値 0.90 を採る。
#
#   ケース    T      P          決め方      H2O    CO2   不活性  CH4転化率   a_C
#   A_base   850  1.04MPaG   転化率       2.777  3.334  18.43%  88.62%   0.394  ← 現状
#   B_900C   900  1.04MPaG   転化率       1.545  2.398  10.39%  88.62%   0.495  ← 圧力据え置き
#   C_03MPaG 900  0.3MPaG    a_C=0.90     0.562  1.867   3.97%  91.19%   0.900
#   D_01MPaG 900  0.1MPaG    a_C=0.90     0.449  1.846   2.09%  95.24%   0.900  ← 不活性最小
#
# A は元の pattern1 の条件そのものなので DryGas は Steam1 を再現する（検算になる）。
CASES = {
    "A_base":   {"T": 850.0, "P": "1.04MPaG", "h2o": 2.776607},
    "B_900C":   {"T": 900.0, "P": "1.04MPaG", "ch4_conversion": CH4_CONVERSION_MIN},
    "C_03MPaG": {"T": 900.0, "P": "0.3MPaG",  "carbon_limit": CARBON_ACTIVITY_DESIGN},
    "D_01MPaG": {"T": 900.0, "P": "0.1MPaG",  "carbon_limit": CARBON_ACTIVITY_DESIGN},
}

#: 出口ストリームの初期推定。Gibbs + 炭素活量の制約を同時に解くとき、既定の
#: 「全成分 1.0」からでは収束しないことがあるため、妥当な合成ガス組成を与える。
SEED_OUTLET = {"CO": 3.5, "H2": 4.9, "CO2": 1.7, "CH4": 0.25, "H2O": 1.0}

#: 検算用: A_base の DryGas はこの値（= plant3 の従来 Steam1）になるはず。
STEAM1_REFERENCE = {"H2": 4.911989, "CO": 3.549828, "CO2": 1.675074,
                    "CH4": 0.243186, "H2O": 0.028800}


def solve_reformer(*, T, P, h2o=None, h2_co=H2_CO_TARGET, carbon_limit=None,
                   ch4_conversion=None, co2_removal=None, verbose=True):
    """H2/CO を指定して CO2 供給量を逆算する形で pattern1 を解き、DryGas を返す。

    co2_removal に除去率 η を渡すと、**改質器の下流・plant3 へ渡す手前**に CO2 除去塔を
    置く（= 除去位置 (a)）。CO2 は改質器の中では CO 源として働かせたうえで、下流の循環系
    には持ち込まない、という狙い。改質器への CO2 供給を絞る C/D ケースと違い、
    ドライリフォーミングの取り分は捨てない。

        DryGas ─[T-100]─→ PlantFeed（plant3 へ）
                  └────→ CO2Vent

    戻り値は (problem, DryGas, CO2f, H2Of, PlantFeed)。除去塔が無いときは
    PlantFeed is DryGas。

    h2o を数値で渡すと水蒸気供給を固定する。自由度 20/20:
        変数   CO2_feed 1 + Mixed 4 + ReactOut 5 + DryGas 5 + Condensate 5
        方程式 Mixer 4 + Gibbs 5 + Separator 5 + 飽和 1 + Henry 4 + 比 1

    h2o=None・carbon_limit=1.0 とすると、**水蒸気供給も未知にして炭素活量が
    その値になる点（＝炭素析出限界）を直接解く**。変数・方程式とも 1 本増えて 21/21。
    a_C ≤ 1 は不等式なので等式系には載らないが、限界そのものを解けば
    「炭素析出しない最小の水蒸気量」が求まる。
    """
    # 水蒸気供給を決める式は 1 本だけ: 固定 / 炭素析出限界 / CH4 転化率 のいずれか
    specs = [h2o is not None, carbon_limit is not None, ch4_conversion is not None]
    if sum(specs) != 1:
        raise ValueError("h2o / carbon_limit / ch4_conversion のいずれか 1 つを指定してください")
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

    # --- CO2 除去塔（オプション・除去位置 (a)）---
    scrub_streams: list[Stream] = []
    scrub_units: list = []
    PlantFeed = DryGas
    if co2_removal is not None:
        PlantFeed = Stream(SPECIES, name="8. PlantFeed", order=8, condition=cold_g,
                           guess=SEED_OUTLET)
        CO2Vent   = Stream(["CO2"], name="9. CO2Vent", order=9, condition=cold_g)
        scrub_streams = [PlantFeed, CO2Vent]
        scrub_units = [Separator(DryGas, [PlantFeed, CO2Vent], name="T-100 CO2除去塔")]

    problem = Problem(
        streams=[RG, CO2f, H2Of, Mixed, ReactOut, DryGas, Condensate] + scrub_streams,
        units=[Mixer([RG, CO2f, H2Of], Mixed, name="M0"),
               GibbsReactor(inlet=Mixed, outlet=ReactOut, species=SPECIES,
                            T=T, P=P, name="G1"),
               Separator(ReactOut, [DryGas, Condensate], name="Condenser")] + scrub_units,
        name=f"Reformer {T:.0f}C {P} H2O={h2o}")
    if co2_removal is not None:
        problem.constrain_recovery(DryGas, scrub_streams[1], {"CO2": float(co2_removal)},
                                   name=f"CO2除去率={co2_removal}")
    problem.constrain_saturation(DryGas, "H2O", T=T_COND, P=P)
    problem.constrain_henry(DryGas, Condensate, DISSOLVED, T=T_COND, P=P)
    # 生成ガスの H2/CO を指定（残り 1 本）
    problem.constrain(DryGas.flow_expr("H2"), h2_co * DryGas.flow_expr("CO"),
                      name=f"H2/CO={h2_co}")
    # 水蒸気供給が未知のとき、それを決める式を 1 本足す
    if carbon_limit is not None:
        problem.constrain_carbon_activity(ReactOut, carbon_limit, T=T, P=P)
    elif ch4_conversion is not None:
        # CH4 転化率 X ⇔ 出口 CH4 = (1 − X) × 供給 CH4（供給は RG_feed で固定）
        problem.constrain(ReactOut.flow_expr("CH4"),
                          (1.0 - ch4_conversion) * RG_FEED["CH4"],
                          name=f"CH4 conversion={ch4_conversion}")

    # 初期推定を入れておくと Gibbs が素直に収束する
    CO2f.molar_flows[:] = 3.0
    if h2o is None:
        # CH4 転化率を高く要求するほど水蒸気も CO2 も桁が上がるので初期推定を合わせる
        scale = 1.0 if ch4_conversion is None else max(1.0, 30.0 * ch4_conversion ** 12)
        H2Of.molar_flows[:] = 2.0 * scale
        CO2f.molar_flows[:] = 3.0 * scale
    sol = problem.solve()
    if not sol.success:
        raise RuntimeError(f"改質器が収束せず: {sol}")
    if verbose:
        print(f"  改質器: T={T:.0f}℃ P={P} H2O={H2Of.flow_of('H2O'):.4f} → "
              f"CO2供給 {CO2f.flow_of('CO2'):.4f} mol/h  (自由度 {problem.degrees_of_freedom()})")
    return problem, DryGas, CO2f, H2Of, PlantFeed


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
    print(f"=== 前工程 pattern1: H2/CO={H2_CO_TARGET:.4f} 維持 + CH4 転化率 ≥ "
          f"{CH4_CONVERSION_MIN*100:.4f}% + a_C ≤ 1 の下での改質条件比較 ===")
    print(f"{'case':>10s} {'T[℃]':>6s} {'P':>9s} {'H2O':>7s} {'CO2feed':>8s} {'DryGas':>8s} "
          f"{'H2/CO':>7s} {'CO2%':>7s} {'CH4%':>7s} {'不活性計':>8s} {'CO':>7s} "
          f"{'CH4転化':>8s} {'a_C':>7s}")
    feeds, summaries = {}, {}
    for tag, cond in CASES.items():
        problem, DryGas, CO2f, H2Of, _ = solve_reformer(verbose=False, **cond)
        react_out = next(s for s in problem.streams if s.name == "5. ReactOut")
        s = drygas_summary(DryGas)
        conv = (RG_FEED["CH4"] - react_out.flow_of("CH4")) / RG_FEED["CH4"]
        a_c = carbon_activity_of(react_out, T=cond["T"], P=cond["P"])
        feeds[tag], summaries[tag] = as_plant_feed(DryGas), dict(s, conv=conv, aC=a_c)
        flag = ("" if (conv >= CH4_CONVERSION_MIN - CH4_CONVERSION_TOL and a_c <= 1.0)
                else "  ⚠制約違反")
        print(f"{tag:>10s} {cond['T']:6.0f} {cond['P']:>9s} {H2Of.flow_of('H2O'):7.4f} "
              f"{CO2f.flow_of('CO2'):8.4f} {s['total']:8.4f} {s['H2/CO']:7.4f} "
              f"{s['CO2%']:7.2f} {s['CH4%']:7.2f} {s['CO2%']+s['CH4%']:8.2f} {s['CO']:7.4f} "
              f"{conv*100:7.3f}% {a_c:7.4f}{flag}")

    # A_base は元の pattern1 条件そのものなので、DryGas は Steam1 と一致するはず
    worst = max(abs(feeds["A_base"][f] - v) for f, v in STEAM1_REFERENCE.items())
    print(f"\n検算: A_base の DryGas と plant3 の従来 Steam1 の最大差 = {worst:.2e} mol/h "
          f"→ {'一致' if worst < 1e-5 else '不一致'}")

    # --- 炭素析出限界そのものを解く（a_C = 1 を課し、水蒸気供給を未知にする）---
    print("\n=== 炭素析出限界（a_C = 1）を直接解く ===")
    print(f"{'T[℃]':>6s}{'H2O[mol/h]':>11s}{'H2O/CH4':>9s}{'CO2供給':>9s}"
          f"{'不活性計':>9s}{'a_C':>9s}")
    for T in (850.0, 900.0):
        pr, dry, co2f, h2of, _ = solve_reformer(T=T, P="1.04MPaG", carbon_limit=1.0,
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
        xlsx = os.path.join(out, f"plant3_{tag}_streams.xlsx")
        to_excel(problem.streams, xlsx, sheet=tag,
                 basis=["mol", "mole_frac", "mass", "normal_volume"])
        print(f"Excel 出力: {xlsx}")

    # --- まとめ ---
    if len(results) > 1:
        print("\n=== 前工程条件が後工程に効く度合い ===")
        print(f"{'指標':>26s}" + "".join(f"{t:>14s}" for t in results))
        rows = [("DryGas 中 CO2 [mol%]", lambda t: summaries[t]["CO2%"]),
                ("DryGas 中 CH4 [mol%]", lambda t: summaries[t]["CH4%"]),
                ("DryGas 不活性計 [mol%]",
                 lambda t: summaries[t]["CO2%"] + summaries[t]["CH4%"]),
                ("改質 CH4 転化率 [%]", lambda t: summaries[t]["conv"] * 100),
                ("改質器出口 a_C", lambda t: summaries[t]["aC"]),
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
