r"""サンプル: 改質器 A〜D — 供給を NL/h で与え直した新しい前提での比較。

example_pattern1_plant3.py とは**別の前提**で組み直したもの。既存ファイルには触れず、
出力も examples/output/reformer/ に分けてある。

────────────────────────────────────────────────────────────────────────
前提（A〜D 共通）
────────────────────────────────────────────────────────────────────────
    供給合計        113 NL/h（= 5.041492 mol/h、22.414 L/mol 換算）
    平衡種          CO2, CH4, H2O, CO, H2 の 5 種（Cantera が Gibbs 平衡を解く）
    凝縮温度        25 ℃（Antoine 飽和 + Henry 溶解で DryGas と凝縮水に分ける）
    液に溶けるガス  H2, CO, CO2, CH4

**A だけ特殊**で、供給 4 成分をすべて実測値で固定する。したがって H2/CO は制約では
なく**結果**として出る（= 1.067）。B〜D は H2/CO = 4/3 を課して供給を逆算する。

────────────────────────────────────────────────────────────────────────
酢酸メチルの量論と H2/CO の狙い
────────────────────────────────────────────────────────────────────────
    3CO + 4H2 → CH3COOCH3 + H2O     H2/CO = 1.3333
    4CO + 3H2 → CH3COOCH3 + CO2     H2/CO = 0.7500

**この 2 式の差は水性ガスシフト 1 回分**（前者 + CO+H2O→CO2+H2 = 後者）。したがって
適正な H2/CO は 0.75〜1.333 の範囲にあり、**系内で順方向 WGS がどれだけ進むかで決まる**。
A の 1.067 はちょうど中間で、WGS が進む前提の設計。B〜D の 1.333 は WGS が進まない
前提（上式のみ）に対応する。

────────────────────────────────────────────────────────────────────────
B〜D の制約と、効く制約の切り替わり
────────────────────────────────────────────────────────────────────────
    (1) 供給合計 113 NL/h
    (2) 原料ガス中の CH4:H2 = 78.9:21.1
    (3) 出口 H2/CO = 4/3
    (4) CH4・CO2 転化率が A を上回る       CH4 87.93% / CO2 41.58%
    (5) a_C ≤ 0.90
    (6) 温度 ≤ 875 ℃、圧力 0.1〜1.41 MPaG

未知は CH4 / CO2 / H2O / H2 の 4 つで、(1)(2)(3) が 3 本。残り 1 本を (4) か (5) の
**どちらで埋めるかは圧力から一意に決まる**（binding_mode()。設計判断ではない）。
**B〜D は同列で、温度と圧力だけが違う 3 点**:

    ≥ 0.417 MPaG   CH4 転化率 88% を課す      a_C は 0.41〜0.90 で余裕
    ≤ 0.417 MPaG   a_C = 0.90 を課す          CH4 転化率は 88% を自然に超える

⚠ ここを取り違えると答えが出ない。**a_C = 0.90 は「炭素析出しない最小の酸化剤量」**を
意味するので、高圧側でこれを課すと水蒸気が最小になり CH4 転化率が 74〜79% にしか
届かず A に負ける。a_C は不等式制約なので、効いていない圧力域では課してはいけない。

交差点 0.417 MPaG は brentq で求めた（両制約が同時に成立する唯一の圧力）。

────────────────────────────────────────────────────────────────────────
結果
────────────────────────────────────────────────────────────────────────
                        A         B         C         D
    温度 [℃]           835       875       875       875
    圧力 [MPaG]        1.41      1.41     0.417       0.1
    供給 CH4 [NL/h]   23.20     28.87     44.95     47.54
    供給 CO2 [NL/h]   52.00     41.25     39.80     40.60
    供給 H2O [NL/h]   32.70     35.15     16.23     12.15
    供給 H2  [NL/h]    5.10      7.72     12.02     12.71
    H2O/CH4           1.409     1.217     0.361     0.256
    O/C               1.818     1.678     1.131     1.059
    CH4 転化率 [%]    87.93     88.00     88.01     94.66
    CO2 転化率 [%]    41.58     52.11     83.94     92.58
    a_C              0.4314    0.4100    0.9000    0.9000
    改質 C 効率 [%]   55.77     66.83     86.09     93.70   ← DryGas CO ÷ 投入 C
    DryGas H2/CO     1.0671    1.3333    1.3333    1.3333
    DryGas [NL/h]    120.13    132.84    183.15    201.41
    DryGas CO2 [mol%] 25.29     14.87      3.49      1.49

**低圧ほど改質器として優れる。** 改質反応はモル数が増えるので低圧が有利。D は投入炭素
の 93.7% を CO に変換し、下流に持ち込む CO2 は 1.49 mol% しかない。

⚠ 下流では逆に働く可能性を予想したが、**実測では外れた**（下記）。

────────────────────────────────────────────────────────────────────────
下流まで通した結果（example_reformer_plant3.py / _plant4.py）
────────────────────────────────────────────────────────────────────────
                        A         B         C         D
    plant3 投入炭素基準 [%]  43.60     57.78     69.16     73.82   ← D が最良
    plant4 投入炭素基準 [%]  44.77     58.92     61.49     61.27
    plant3 酢酸メチル      0.4876    0.6026    0.8717    0.9675   mol/h
    plant4 酢酸メチル      0.5007    0.6144    0.7750    0.8030   mol/h

**改質器の優位がそのまま下流まで通った。** A → D で投入炭素基準の収率 1.7 倍、
酢酸メチル 2.0 倍。

事前に「D は DryGas の CO2 が 1.5 mol% しかないので、ループ CO2 が減って順方向 WGS が
CO を食い不利になる」と予想したが外れた:

    反応器入口 CO2 [mol%]  67.42 → 22.58   狙いどおり減る
    CO2 正味生成 [mol/h]    0.180 → 0.404   WGS は確かに進んだ
    酢酸メチル              0.488 → 0.968   それでも 2 倍

WGS による CO 損失は確かに増えたが、**改質器が供給する CO の絶対量が 41.9 → 82.6 NL/h
と倍増する効果が圧倒的に上回った**。過去の検討（example_co2_removal.py）は改質器の
供給を固定したまま下流で CO2 を抜いていたので損失だけが見えたが、今回は改質器側で
CO を増やしているので構図が違う。

**plant3 と plant4 が C で逆転する。** A・B は plant4 が勝ち、C・D は plant3 が勝つ:

    DME パージ [mol/h]   plant3: 0.004 → 0.111   plant4: 0.003 → 0.374

plant4 は未反応 DME をあえて循環させる設計なので、ループ CO2 が減ると DME が濃縮して
パージで失われる。D では plant4 の DME 損失が plant3 の 3.4 倍。R2 入口の CO 濃度も
9.45 → 5.87 mol% と下がり、カルボニル化が追いつかず DME が溜まる構図。過去の知見
（plant4 は「CO2 希釈剤を失うこと」全般に弱い）と一致する。

⚠ **指標は投入炭素基準で見ること。** 新鮮 CO 基準だと B が最良に見えてしまう:

    新鮮CO基準 [%]    A 78.18   B 86.47   C 80.33   D 78.78   ← 誤読
    投入炭素基準 [%]  A 43.60   B 57.78   C 69.16   D 73.82   ← 正しい

分母に改質器で消費した CH4・CO2 を含めないと、D の「原料を有効に使っている」という
最大の長所が消える。

────────────────────────────────────────────────────────────────────────
実行方法
────────────────────────────────────────────────────────────────────────
改質器だけなら Cantera のみで数秒。

    PYTHONPATH=. python3 -u examples/example_reformer.py

下流まで通す場合は reaction_rate も要る（1 ケース 5〜20 分）。

    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_reformer_plant3.py
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_reformer_plant4.py

出力は examples/output/reformer/ 以下（既存の検討と混ざらないよう分けてある）。
"""

import os

from chemflow2 import (
    GibbsReactor,
    Mixer,
    Problem,
    Separator,
    Stream,
    StreamCondition,
    export_mermaid,
    stream_table,
    to_excel,
)
from chemflow2.core.carbon_activity import carbon_activity_of
from chemflow2.core.pressure import parse_pressure

#: ノルマルモル体積 [L/mol]（0 ℃・1 atm）。NL/h ↔ mol/h の換算に使う。
NV = 22.414

OUT = os.path.join(os.path.dirname(__file__), "output", "reformer")

SPECIES = ["CO2", "CH4", "H2O", "CO", "H2"]   # Gibbs 平衡種
DISSOLVED = ["H2", "CO", "CO2", "CH4"]        # 液水に溶けるガス
FEED_C = ["CH4", "CO2", "H2O", "H2"]
T_COND = 25.0

#: 供給合計 [NL/h] — A〜D すべてで固定（スループットを揃えて比較するため）
TOTAL_NL = 113.0
#: 原料ガス中の CH4:H2（B〜D。A は実測値なので 82.0:18.0 とずれる）
CH4_FRAC, H2_FRAC = 0.789, 0.211
#: B〜D が狙う出口 H2/CO。3CO + 4H2 → MA + H2O の量論そのもの。
H2_CO_TARGET = 4.0 / 3.0
#: 炭素活量の設計値（千代田 WO2012140994A1 表 5 の実施例に倣う）
CARBON_ACTIVITY_DESIGN = 0.90

#: A の供給 [NL/h] — 実測値。4 成分すべて固定するので H2/CO は結果として出る。
A_FEED_NL = {"CH4": 23.2, "CO2": 52.0, "H2O": 32.7, "H2": 5.1}

#: A が達成する転化率。B〜D はこれを上回ることが要件。
A_CH4_CONV, A_CO2_CONV = 87.93, 41.58

#: 効く制約が切り替わる圧力 [MPaG]。brentq で求めた値。
#: これ以上なら CH4 転化率が、これ以下なら a_C が水蒸気量を決める。
P_CROSSOVER = 0.417

#: B〜D に課す CH4 転化率（A の 87.93% をわずかに上回る）
CH4_CONVERSION_MIN = 0.88

#: ケース定義。**B〜D は温度と圧力だけを指定する**。水蒸気量を決める制約は
#: 圧力から一意に決まる従属量なので、手で書かない（binding_mode 参照）。
#: A だけが本当に特殊で、供給 4 成分をすべて固定するため H2/CO が結果になる。
CASES = {
    "A": dict(T=835.0, P="1.41MPaG", feed_nl=A_FEED_NL,
              note="実測供給を全固定。H2/CO は結果（1.067）"),
    "B": dict(T=875.0, P="1.41MPaG", note="A と同圧"),
    "C": dict(T=875.0, P=f"{P_CROSSOVER}MPaG", note="交差点。両制約が同時に成立"),
    "D": dict(T=875.0, P="0.1MPaG", note="指定範囲の下限。転化率が最大"),
}


def binding_mode(P):
    """圧力から、水蒸気量を決める 1 本の制約を選ぶ。

    **これは設計判断ではなく従属量**。B〜D は「圧力だけが違う 3 点」であり、
    どちらの制約が効くかは圧力で決まる:

        P >  0.417 MPaG   CH4 転化率が効く（a_C は自動的に 0.90 未満に収まる）
        P <= 0.417 MPaG   a_C が効く（CH4 転化率は自動的に 88% を超える）

    ⚠ 取り違えると答えが出ない。a_C = 0.90 は「炭素析出しない最小の酸化剤量」を
    意味するので、高圧側でこれを課すと水蒸気が最小になり CH4 転化率が 74〜79% に
    しか届かず A に負ける。a_C は不等式制約なので、効いていない圧力域では
    課してはいけない。
    """
    return "aC" if parse_pressure(P) <= parse_pressure(f"{P_CROSSOVER}MPaG") \
        else CH4_CONVERSION_MIN

#: 低圧側は初期推定が悪いと縮退解（CH4→0）に落ちる。高圧から段階的にウォームスタート
#: して降りる必要がある。C・D を解く前に通す圧力の列。
CONTINUATION = ["1.41MPaG", "1.0MPaG", "0.7MPaG", "0.5MPaG",
                f"{P_CROSSOVER}MPaG", "0.3MPaG", "0.2MPaG", "0.15MPaG", "0.1MPaG"]

#: 出力の basis（4 種。plant3/plant4 側と揃えてある）
BASIS = ["mol", "mole_frac", "mass", "normal_volume"]


def build(T, P, *, mode, feed_nl=None, seed=None):
    """改質器のフローシートを組む。

    自由度:
        A（mode="fixed"）  変数 19 / 方程式 19
            供給が既知なので Mixed 4 + ReactOut 5 + DryGas 5 + Condensate 5 = 19、
            Mixer 4 + Gibbs 5 + Separator 5 + 飽和 1 + Henry 4 = 19。
        B〜D               変数 23 / 方程式 23
            供給 4 が未知になり、合計 1 + 比 1 + H2/CO 1 + 水蒸気を決める 1 が増える。
    """
    hot = StreamCondition(T=T, P=P, phase="gas")
    cold_g = StreamCondition(T=T_COND, P=P, phase="gas")
    cold_l = StreamCondition(T=T_COND, P=P, phase="liquid")
    g = (seed or {}).get

    Feed = Stream(FEED_C, name="1. Feed", order=1, condition=hot,
                  flows={k: v / NV for k, v in feed_nl.items()} if feed_nl else None,
                  guess=None if feed_nl else g("1. Feed", {"CH4": 1.9, "CO2": 1.8,
                                                           "H2O": 0.8, "H2": 0.5}))
    Mixed = Stream(FEED_C, name="2. Mixed", order=2, condition=hot, guess=g("2. Mixed"))
    ReactOut = Stream(SPECIES, name="3. ReactOut", order=3, condition=hot,
                      guess=g("3. ReactOut", {"CO": 2.5, "H2": 3.3, "CO2": 0.8,
                                              "CH4": 0.2, "H2O": 0.6}))
    DryGas = Stream(SPECIES, name="4. DryGas", order=4, condition=cold_g,
                    guess=g("4. DryGas", {"CO": 2.5, "H2": 3.3, "CO2": 0.8,
                                          "CH4": 0.2, "H2O": 0.02}))
    Condensate = Stream(SPECIES, name="5. Condensate", order=5, condition=cold_l,
                        guess=g("5. Condensate"))

    problem = Problem(
        streams=[Feed, Mixed, ReactOut, DryGas, Condensate],
        units=[Mixer([Feed], Mixed, name="M0"),
               GibbsReactor(inlet=Mixed, outlet=ReactOut, species=SPECIES,
                            T=T, P=P, name="G1"),
               Separator(ReactOut, [DryGas, Condensate], name="Condenser")],
        name=f"Reformer {T:.0f}C {P}")
    problem.constrain_saturation(DryGas, "H2O", T=T_COND, P=P)
    problem.constrain_henry(DryGas, Condensate, DISSOLVED, T=T_COND, P=P)

    if feed_nl is None:
        problem.constrain(Feed.total_flow, TOTAL_NL / NV, name=f"合計 {TOTAL_NL} NL/h")
        problem.constrain(Feed.flow_expr("CH4") * H2_FRAC,
                          Feed.flow_expr("H2") * CH4_FRAC, name="CH4:H2=78.9:21.1")
        problem.constrain(DryGas.flow_expr("H2"), H2_CO_TARGET * DryGas.flow_expr("CO"),
                          name="H2/CO=4/3")
        if mode == "aC":
            problem.constrain_carbon_activity(ReactOut, CARBON_ACTIVITY_DESIGN, T=T, P=P)
        else:
            problem.constrain(ReactOut.flow_expr("CH4"),
                              (1.0 - mode) * Feed.flow_expr("CH4"),
                              name=f"CH4 転化率={mode}")
    return problem, (Feed, Mixed, ReactOut, DryGas, Condensate)


def _solve(problem, fixed_feed):
    """A は供給既知なので root で素直に解ける。B〜D は非負制約つきで解く。"""
    if fixed_feed:
        sol = problem.solve(bounds=None, tol=1e-7)
    else:
        sol = problem.solve(bounds=(0, 1e3), tol=1e-7,
                            ftol=1e-14, xtol=1e-14, gtol=1e-14, x_scale="jac")
    if not sol.success:
        raise RuntimeError(f"{problem.name} が収束せず: {sol}")
    return sol


def solve_case(tag, verbose=False):
    """1 ケース解いて (problem, streams) を返す。

    低圧ケースは初期推定が悪いと縮退解（CH4→0）に落ちるので、目標圧力までの
    CONTINUATION を高圧側から順に解いてウォームスタートを作る。各段で効く制約は
    その段の圧力から binding_mode() が決める（目標圧力のものを流用しない）。
    """
    cfg = CASES[tag]
    if "feed_nl" in cfg:                      # A: 供給を全固定
        pb, st = build(cfg["T"], cfg["P"], mode="fixed", feed_nl=cfg["feed_nl"])
        _solve(pb, True)
        return pb, st

    target = cfg["P"]
    seed = None
    for P in CONTINUATION:
        pb, st = build(cfg["T"], P, mode=binding_mode(P), seed=seed)
        try:
            _solve(pb, False)
        except RuntimeError:
            if P == target:
                raise
            continue                      # 途中の圧力は失敗しても構わない
        if st[0].flow_of("CH4") < 0.1:    # 縮退解は seed にしない
            continue
        seed = {s.name: s.molar_flows.copy() for s in pb.streams}
        if verbose:
            print(f"    継続: {P}", flush=True)
        if P == target:
            return pb, st
    # CONTINUATION に目標圧が無い場合は直接解く
    pb, st = build(cfg["T"], target, mode=binding_mode(target), seed=seed)
    _solve(pb, False)
    return pb, st


def summary(tag, problem, streams):
    """1 ケースの全指標を dict で返す。"""
    Feed, _Mixed, ReactOut, DryGas, Condensate = streams
    cfg = CASES[tag]
    f = {c: Feed.flow_of(c) for c in FEED_C}
    C = f["CH4"] + f["CO2"]
    O = 2 * f["CO2"] + f["H2O"]
    H = 4 * f["CH4"] + 2 * f["H2O"] + 2 * f["H2"]
    dg = float(DryGas.total_flow.eval())
    return {
        "温度 [℃]": cfg["T"], "圧力": cfg["P"],
        **{f"供給 {c} [NL/h]": f[c] * NV for c in FEED_C},
        "供給 合計 [NL/h]": sum(f.values()) * NV,
        "CH4/(CH4+H2) [%]": f["CH4"] / (f["CH4"] + f["H2"]) * 100,
        "H2O/CH4 [-]": f["H2O"] / f["CH4"], "CO2/CH4 [-]": f["CO2"] / f["CH4"],
        "O/C [-]": O / C, "H/C [-]": H / C,
        "CH4 転化率 [%]": (f["CH4"] - ReactOut.flow_of("CH4")) / f["CH4"] * 100,
        "CO2 転化率 [%]": (f["CO2"] - DryGas.flow_of("CO2")) / f["CO2"] * 100,
        "炭素活量 a_C [-]": carbon_activity_of(ReactOut, T=cfg["T"], P=cfg["P"]),
        "改質 C 効率 [%]": DryGas.flow_of("CO") / C * 100,
        "DryGas H2/CO [-]": DryGas.flow_of("H2") / DryGas.flow_of("CO"),
        "DryGas 合計 [NL/h]": dg * NV,
        **{f"DryGas {c} [NL/h]": DryGas.flow_of(c) * NV for c in ("CO", "H2", "CO2", "CH4", "H2O")},
        **{f"DryGas {c} [mol%]": DryGas.flow_of(c) / dg * 100 for c in ("CO", "H2", "CO2", "CH4")},
        "凝縮水 [NL/h]": float(Condensate.total_flow.eval()) * NV,
    }


def plant_feed(DryGas):
    """DryGas を下流プラントの新鮮供給の形（mol/h）に直す。"""
    return {c: DryGas.flow_of(c) for c in DryGas.formulas}


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"=== 改質器 A〜D（供給合計 {TOTAL_NL} NL/h 固定）===\n")
    res = {}
    for tag in CASES:
        print(f"--- {tag}: {CASES[tag]['note']}", flush=True)
        pb, st = solve_case(tag, verbose=True)
        res[tag] = summary(tag, pb, st)
        print(stream_table(pb.streams, basis=["normal_volume", "mole_frac"]))
        to_excel(pb.streams, os.path.join(OUT, f"reformer_{tag}.xlsx"),
                 sheet=f"reformer_{tag}", basis=BASIS, components=SPECIES)
        export_mermaid(pb, os.path.join(OUT, f"reformer_{tag}.html"),
                       title=f"Reformer {tag} ({CASES[tag]['T']:.0f}C {CASES[tag]['P']})",
                       style="diamond")
        print()

    print("=" * 78)
    print("=== A〜D 比較")
    print("=" * 78)
    print(f"{'項目':>22s}" + "".join(f"{t:>14s}" for t in res))
    for k in next(iter(res.values())):
        vals = []
        for t in res:
            v = res[t][k]
            vals.append(f"{v:14.4f}" if isinstance(v, float) else f"{str(v):>14s}")
        print(f"{k:>22s}" + "".join(vals))
    print(f"\n出力: {OUT}")


if __name__ == "__main__":
    main()
