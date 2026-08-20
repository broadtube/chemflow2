"""サンプル: example_plant3 の反応器をタンデムから**2段（中間で水/MeOH を除去）**に分けた版。

plant3 との違い（反応網・供給組成・パージ率は同じ）:
    plant3:  M1 → R1(合成+脱水) → R2(カルボニル化) → Condenser → …
    plant4:  M1 → R1(合成+脱水) → Cond1 ┬ ガス → R2(カルボニル化) → Cond2 ┬ ガス → M1
                                        └ 液(H2O,MeOH) → Column1 → Column2  └ 液 → Column3 → MA

狙い: **H-MOR のカルボニル化は水で強く阻害される**。R1 の直後に凝縮器を置いて H2O と
MeOH を抜いてから H-MOR 床に入れることで、水の分圧を下げた状態でカルボニル化させる。
未反応 DME は Cond2 のガス側に残して M1 へ戻し、R1→Cond1→R2 と再循環させる。

触媒量は **各床とも自分の入口基準で SV 5000/h**（plant3 は「タンデム合計」を反応器入口
基準で SV 5000/h だった）。R1 の入口は ReactorIn、R2 の入口は Cond1 のガスなので、
2つの触媒体積を同時に満たす必要がある → 成分ごとの割線法で外側反復する。

要 reaction_rate: pip install -e ../reaction_rate
実行: PYTHONPATH=.:../reaction_rate/src python3 examples/example_plant4.py
"""

import os
import time

import numpy as np

from chemflow2 import (
    KineticReactor,
    Mixer,
    Problem,
    Separator,
    Splitter,
    Stream,
    StreamCondition,
    export_mermaid,
    generate_mermaid,
    stream_table,
    to_csv,
    to_excel,
)

# --- 成分（反応網7種 + 不活性・微量成分）: plant3 と同一 ---
C = ["H2", "CO", "CO2", "CH4", "H2O", "CH3OH", "CH3OCH3", "CH3COOCH3",
     "CH3CHO", "CH3COOH", "N2"]

# --- 新鮮供給 Steam1 [mol/h]: plant3 と同一 ---
STEAM1 = {"H2": 4.911989, "CO": 3.549828, "CO2": 1.675074, "CH4": 0.243186,
          "H2O": 0.0288, "CH3CHO": 0.0, "CH3COOH": 0.0, "N2": 0.0}

# --- 触媒の設計基準（reaction_rate の fig_tandem_sty_stack_zsm5 と同一）---
SV = 5000.0                                      # [1/h] ← 各床とも自分の入口基準
RHO_SYN, RHO_DEH, RHO_MA = 1300.0, 800.0, 837.0  # 充填密度 [kg/m³]
F_SYN = 1.0 / 1.9          # ハイブリッド床の体積分率（synthesis 側）
VM_STP = 22.414e-3         # 標準状態モル体積 [m³/mol]
PURGE = 0.05               # Cond2 ガスのパージ率

gas = StreamCondition(T=250, P="5MPaG", phase="gas")   # 5 MPaG + atm = 51.0 bar
liq = StreamCondition(T=40, P="5MPaG", phase="liquid")

NAMES = ["1. Steam1", "2. ReactorIn", "3. R1Out", "4. Cond1Gas", "5. Cond1Liq",
         "6. Col1Top", "7. Sep1Liq", "8. RecycleMeOH", "9. Water", "10. R2Out",
         "11. Cond2Gas", "12. Purge", "13. RecycleGas", "14. Cond2Liq",
         "15. MethylAcetate", "16. Col3Bottoms",
         # 以下は co2_removal を指定したときだけ現れる（既定の 16 本には影響しない）
         "17. ScrubbedGas", "18. CO2Vent"]


def hybrid_masses(v1: float) -> tuple[float, float]:
    """ハイブリッド床の体積 [m³] → (m_synthesis, m_dehydration) [kg]。"""
    return RHO_SYN * v1 * F_SYN, RHO_DEH * v1 * (1 - F_SYN)


def volume_for(n_molh: float) -> float:
    """入口モル流量 [mol/h] → SV を満たす触媒体積 [m³]。"""
    return n_molh * VM_STP / SV


def build(v1: float, v2: float, seed: dict[str, np.ndarray] | None = None,
          feed: dict[str, float] | None = None,
          co2_removal: float | None = None, purge: float = PURGE,
          co2_position: str = "recycle"):
    """ハイブリッド床 v1 [m³]・カルボニル化床 v2 [m³] のフローシートを組む。

    feed を渡すと新鮮供給の組成を差し替えられる（前工程 pattern1 の DryGas を
    そのまま流し込むため。既定は STEAM1）。

    co2_removal に除去率 η を渡すと CO2 除去塔を 1 段挿す。plant4 は凝縮器が
    2 つあるので、**どちらに付けるかを co2_position で選ぶ**:

        "recycle"      Cond2Gas と SP-gas の間（= plant3 の (b) と同じ位置）
                       Cond2Gas ─[T-101]→ ScrubbedGas ─[SP-gas]→ Purge / RecycleGas

        "interstage"   Cond1Gas と R2 の間（**plant4 でのみ可能な位置**）
                       Cond1 は既に水を抜く場所なので、そこで CO2 も抜けば
                       カルボニル化床に入る CO 分圧を直接上げられる。H-MOR は
                       水で強く阻害されるうえ CO2 で希釈されると不利なので、
                       plant3 には無いこの位置に効く見込みがある。
                       Cond1Gas ─[T-102]→ ScrubbedGas → R2

    塔はガス側に対しては「CO2 だけを η の割合で抜く分離器」でしかないので、
    Separator 1 個 + 回収率制約 1 本で表す（変数・式とも +12 で閉じる）。
    苛性ソーダの量論は循環系に影響しないため、ここには入れない。
    """
    if co2_position not in ("recycle", "interstage"):
        raise ValueError('co2_position は "recycle" か "interstage"')

    def S(i: int, cond: StreamCondition = gas) -> Stream:
        nm = NAMES[i - 1]
        return Stream(C, name=nm, order=i, condition=cond,
                      guess=seed.get(nm) if seed else None)

    Steam1      = Stream(C, name=NAMES[0], order=1, condition=gas,
                         flows=feed if feed is not None else STEAM1)
    ReactorIn   = S(2)
    R1Out       = S(3)          # ハイブリッド床（合成＋脱水）の出口
    Cond1Gas    = S(4)          # → カルボニル化床へ（水/MeOH を抜いた後）
    Cond1Liq    = S(5, liq)     # H2O + MeOH
    Col1Top     = S(6, liq)     # plant3 の Column1 留出に相当（plant4 では 0 になる）
    Sep1Liq     = S(7, liq)
    RecycleMeOH = S(8)          # → M1（パージ無し）
    Water       = S(9, liq)
    R2Out       = S(10)         # カルボニル化床の出口
    Cond2Gas    = S(11)
    Purge       = S(12)
    RecycleGas  = S(13)         # → M1（未反応 DME もここに乗る）
    Cond2Liq    = S(14, liq)
    MA_Product  = S(15, liq)    # Column3 留出 = 酢酸メチル
    Col3Bottoms = S(16, liq)    # Column3 缶出 = 系外へ

    # --- CO2 除去塔（オプション）---
    # 塔を入れると、その下流ユニットの入口が ScrubbedGas に差し替わる。
    scrub_streams: list[Stream] = []
    scrub_units: list = []
    r2_feed, split_feed = Cond1Gas, Cond2Gas
    if co2_removal is not None:
        ScrubbedGas = S(17)
        CO2Vent     = Stream(["CO2"], name=NAMES[17], order=18, condition=gas,
                             guess=seed.get(NAMES[17]) if seed else None)
        scrub_streams = [ScrubbedGas, CO2Vent]
        src = Cond1Gas if co2_position == "interstage" else Cond2Gas
        tag = "T-102 CO2除去塔(段間)" if co2_position == "interstage" else "T-101 CO2除去塔(循環)"
        scrub_units = [Separator(src, [ScrubbedGas, CO2Vent], name=tag)]
        if co2_position == "interstage":
            r2_feed = ScrubbedGas
        else:
            split_feed = ScrubbedGas

    m_syn, m_deh = hybrid_masses(v1)
    m_ma = RHO_MA * v2

    M1   = Mixer([Steam1, RecycleGas, RecycleMeOH], ReactorIn, name="M1")
    R1   = KineticReactor(inlet=ReactorIn, outlet=R1Out,
                          masses={"synthesis": m_syn, "dehydration": m_deh},
                          models={"synthesis": "KOGAS", "dehydration": "ZSM5"},
                          k_eq3="thermo", n_points=2, name="R1 hybrid")
    Cond1 = Separator(R1Out, [Cond1Gas, Cond1Liq], name="Cond1")
    Col1  = Separator(Cond1Liq, [Col1Top, Sep1Liq], name="Column1")
    Col2  = Separator(Sep1Liq, [RecycleMeOH, Water], name="Column2")
    R2    = KineticReactor(inlet=r2_feed, outlet=R2Out,
                           masses={"carbonylation": m_ma},
                           models={"carbonylation": "DTU-Cheung2007-2"},
                           T=250, P="5MPaG", n_points=2, name="R2 MA bed")
    Cond2 = Separator(R2Out, [Cond2Gas, Cond2Liq], name="Cond2")
    SPg   = Splitter(split_feed, [Purge, RecycleGas], ratios=[purge, 1 - purge], name="SP-gas")
    Col3  = Separator(Cond2Liq, [MA_Product, Col3Bottoms], name="Column3")

    streams = [Steam1, ReactorIn, R1Out, Cond1Gas, Cond1Liq, Col1Top, Sep1Liq,
               RecycleMeOH, Water, R2Out, Cond2Gas, Purge, RecycleGas, Cond2Liq,
               MA_Product, Col3Bottoms] + scrub_streams
    problem = Problem(streams=streams,
                      units=[M1, R1, Cond1, Col1, Col2, R2, Cond2, SPg, Col3] + scrub_units,
                      name="Syngas plant (2-stage: hybrid -> knockout -> carbonylation)")

    # 除去塔の指定: CO2 のみ η の割合で CO2Vent へ。CO2Vent は 1 成分ストリームなので
    # 「他成分は回収率 0」を書く必要がない（書くと変数の無い式が増えて自由度が崩れる）。
    if co2_removal is not None:
        src = Cond1Gas if co2_position == "interstage" else Cond2Gas
        problem.constrain_recovery(src, scrub_streams[1], {"CO2": float(co2_removal)},
                                   name=f"CO2除去率={co2_removal}")

    # --- 分離指定（各分離器の「第2出口」への回収率）---
    # Cond1: H2/CO/CO2/DME/CH4/N2 はガス、H2O/MeOH は液（指定どおり）
    problem.constrain_recovery(R1Out, Cond1Liq, {
        "H2": 0, "CO": 0, "CO2": 0, "CH3OCH3": 0, "CH4": 0, "N2": 0,
        "H2O": 1, "CH3OH": 1, "CH3COOCH3": 1, "CH3CHO": 1, "CH3COOH": 1,
    })
    # Column1: plant3 と同じ仕様（MA と CH3CHO を留出）。plant4 では留出が 0 になる
    problem.constrain_recovery(Cond1Liq, Sep1Liq, {
        "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "N2": 0, "CH3OCH3": 0,
        "CH3COOCH3": 0, "CH3CHO": 0, "H2O": 1, "CH3OH": 1, "CH3COOH": 1,
    })
    # Column2: メタノールを留出（全量リサイクル）、水と酢酸は缶出
    problem.constrain_recovery(Sep1Liq, Water, {
        "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "N2": 0, "CH3OCH3": 0,
        "CH3COOCH3": 0, "CH3CHO": 0, "CH3OH": 0, "H2O": 1, "CH3COOH": 1,
    })
    # Cond2: 酢酸メチル等は液、**未反応 DME はガス側**に残して M1 へ戻す
    problem.constrain_recovery(R2Out, Cond2Liq, {
        "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "N2": 0, "CH3OCH3": 0,
        "CH3COOCH3": 1, "H2O": 1, "CH3OH": 1, "CH3CHO": 1, "CH3COOH": 1,
    })
    # Column3: 酢酸メチルだけを留出、それ以外は缶出（系外へ）
    problem.constrain_recovery(Cond2Liq, Col3Bottoms, {
        "CH3COOCH3": 0,
        "H2": 1, "CO": 1, "CO2": 1, "CH4": 1, "N2": 1, "CH3OCH3": 1,
        "H2O": 1, "CH3OH": 1, "CH3CHO": 1, "CH3COOH": 1,
    })
    return problem, streams


# 反応器入口 / 新鮮供給 の初期推定（循環ぶん）。外側反復は割線法なので外れても収束する。
RECYCLE_GUESS = 6.5
GAS_FRACTION_GUESS = 0.95    # Cond1 のガス側に残る割合の初期推定


def solve_with_sv(max_outer: int = 10, tol: float = 1e-4, verbose: bool = True,
                  feed: dict[str, float] | None = None,
                  co2_removal: float | None = None, purge: float = PURGE,
                  co2_position: str = "recycle", solve_tol: float = 1e-4,
                  progress_every: int = 0):
    """2つの床がそれぞれ入口基準 SV を満たすまで触媒体積 (V1, V2) を外側反復する。

    g(V) = [volume_for(ReactorIn) − V1, volume_for(R2の入口) − V2] の零点を、
    成分ごとの割線法で探す（各評価が1回のフローシート求解なので回数を切り詰める）。

    feed / co2_removal / purge / co2_position は build() にそのまま渡す。
    solve_tol は Problem.solve の残差ノルム判定（既定 1e-4 の理由は g() のコメント）。
    """
    state: dict[str, object] = {}

    def g(v: np.ndarray, it: int) -> np.ndarray:
        t0 = time.time()
        problem, streams = build(v[0], v[1], state.get("seed"), feed=feed,
                                 co2_removal=co2_removal, purge=purge,
                                 co2_position=co2_position)
        if it == 1 and verbose:
            print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
        # tol は残差ノルムの合格判定。PFR は rtol=1e-8 で積分しているので残差には
        # 流量 × 1e-8 程度のノイズ床があり、方程式 165 本ではノルムが ~3e-6 で頭打ち
        # になる。1e-4 mol/h は流量 60 mol/h に対して相対 2e-6 で、十分すぎる精度。
        sol = problem.solve(bounds=(0, np.inf), tol=solve_tol,
                            progress_every=progress_every, progress_label=f"外側{it}",
                            stop_at_tol=stop_at_tol,
                            ftol=solver_tols, xtol=solver_tols, gtol=solver_tols)
        if not sol.success:
            raise RuntimeError(f"外側反復 {it} で収束せず: {sol}")
        n1 = float(streams[1].total_flow.eval())    # ReactorIn
        # 第2床の SV 基準は「R2 の実際の入口」。段間に除去塔を挿すと Cond1Gas では
        # なく ScrubbedGas が R2 の入口になるので、そちらを見ないと SV がずれる。
        r2_in = next(u.inlets[0] for u in problem.units if u.name == "R2 MA bed")
        n2 = float(r2_in.total_flow.eval())
        v_new = np.array([volume_for(n1), volume_for(n2)])
        state["seed"] = {s.name: s.molar_flows.copy() for s in streams}
        state["problem"], state["streams"] = problem, streams
        n_pfr = sum(u.n_pfr_calls for u in problem.units if hasattr(u, "n_pfr_calls"))
        if verbose:
            print(f"  [{it}] V=({v[0]*1e6:7.3f}, {v[1]*1e6:7.3f}) mL → "
                  f"入口=({n1:7.3f}, {n2:7.3f}) mol/h → "
                  f"V_SV=({v_new[0]*1e6:7.3f}, {v_new[1]*1e6:7.3f}) mL   "
                  f"({time.time()-t0:.1f}s, PFR積分 {n_pfr}回, nfev={sol.nfev})")
        return v_new - v

    n_feed = sum((feed if feed is not None else STEAM1).values())
    v_prev = np.array([volume_for(RECYCLE_GUESS * n_feed),
                       volume_for(GAS_FRACTION_GUESS * RECYCLE_GUESS * n_feed)])
    g_prev = g(v_prev, 1)
    if np.all(np.abs(g_prev) <= tol * v_prev):
        return state["problem"], state["streams"], v_prev

    v = v_prev + g_prev                     # 1 歩目は固定点ステップ
    for it in range(2, max_outer + 1):
        g_cur = g(v, it)
        if np.all(np.abs(g_cur) <= tol * v):
            return state["problem"], state["streams"], v
        # 成分ごとの割線更新（分母が潰れたら固定点ステップにフォールバック）
        v_next = np.empty_like(v)
        for i in range(len(v)):
            d = g_cur[i] - g_prev[i]
            step = -g_cur[i] * (v[i] - v_prev[i]) / d if d != 0.0 else g_cur[i]
            v_next[i] = np.clip(v[i] + step, 0.2 * v[i], 5.0 * v[i])
        v_prev, g_prev = v, g_cur
        v = v_next
    raise RuntimeError("触媒体積の外側反復が収束しませんでした")


def main():
    print(f"=== plant4: 2段反応器（250℃ / 51.0 bar / 各床とも入口基準 SV {SV:.0f}/h）===")
    problem, streams, v = solve_with_sv()
    (Steam1, ReactorIn, R1Out, Cond1Gas, Cond1Liq, Col1Top, Sep1Liq, RecycleMeOH,
     Water, R2Out, Cond2Gas, Purge, RecycleGas, Cond2Liq, MA_Product,
     Col3Bottoms) = streams

    m_syn, m_deh = hybrid_masses(v[0])
    m_ma = RHO_MA * v[1]
    print("\n--- 触媒（収束値）---")
    print(f"R1 ハイブリッド床 : {v[0]*1e6:.4f} mL  "
          f"(synthesis {m_syn*1e3:.4f} g / dehydration {m_deh*1e3:.4f} g)")
    print(f"R2 カルボニル化床 : {v[1]*1e6:.4f} mL  (carbonylation {m_ma*1e3:.4f} g)")
    print(f"合計              : {(v[0]+v[1])*1e6:.4f} mL")

    print()
    print(stream_table(problem.streams, basis=["mol", "mole_frac"]))

    # --- 性能 ---
    co_feed = Steam1.flow_of("CO")
    ma = MA_Product.flow_of("CH3COOCH3")
    print("\n--- 性能 ---")
    print(f"R1 CO 1パス転化率 : "
          f"{(ReactorIn.flow_of('CO') - R1Out.flow_of('CO'))/ReactorIn.flow_of('CO')*100:.2f} %")
    print(f"R2 DME 1パス転化率: "
          f"{(Cond1Gas.flow_of('CH3OCH3') - R2Out.flow_of('CH3OCH3'))/Cond1Gas.flow_of('CH3OCH3')*100:.2f} %")
    print(f"酢酸メチル生成     : {ma:.4f} mol/h（新鮮 CO 基準の炭素収率 {ma*3/co_feed*100:.2f} %）")
    print(f"CO パージ損失      : {Purge.flow_of('CO'):.4f} mol/h")

    # --- 水の除去効果（plant4 の狙い）---
    print("\n--- カルボニル化床入口の水分（plant4 の狙い）---")
    print(f"R1 出口の H2O      : {R1Out.flow_of('H2O'):.4f} mol/h "
          f"({R1Out.flow_of('H2O')/float(R1Out.total_flow.eval())*100:.3f} mol%)")
    print(f"Cond1 ガスの H2O   : {Cond1Gas.flow_of('H2O'):.4f} mol/h "
          f"({Cond1Gas.flow_of('H2O')/float(Cond1Gas.total_flow.eval())*100:.3f} mol%)")

    # --- 不活性・CO2 の蓄積 ---
    print("\n--- 循環による蓄積（パージ率 {:.0%}）---".format(PURGE))
    for f in ("CH4", "CO2"):
        print(f"{f:4s}: 供給 {Steam1.flow_of(f):7.4f} → 反応器入口 {ReactorIn.flow_of(f):8.4f} mol/h "
              f"({ReactorIn.flow_of(f)/float(ReactorIn.total_flow.eval())*100:5.2f} mol%)"
              f"  パージ {Purge.flow_of(f):7.4f} mol/h")

    # --- メタノール収支 ---
    print("\n--- メタノール収支 ---")
    print(f"ハイブリッド床 正味生成: {R1Out.flow_of('CH3OH') - ReactorIn.flow_of('CH3OH'):.4f} mol/h")
    print(f"循環 MeOH             : {RecycleMeOH.flow_of('CH3OH'):.4f} mol/h（パージ無し）")
    print(f"Column1 留出（Col1Top）: {float(Col1Top.total_flow.eval()):.4e} mol/h"
          "  ← plant4 では 0 になる（MA が Cond1 より下流で生成するため）")

    # --- ファイル出力 ---
    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)
    csv_path = os.path.join(out, "plant4_streams.csv")
    to_csv(problem.streams, csv_path, basis=["mol", "mole_frac", "mass", "normal_volume"])
    print(f"\nCSV 出力  : {csv_path}")
    try:
        xlsx = os.path.join(out, "plant4_streams.xlsx")
        to_excel(problem.streams, xlsx, basis=["mol", "mole_frac", "mass", "normal_volume"])
        print(f"Excel 出力: {xlsx}")
    except ImportError:
        print("Excel 出力スキップ（openpyxl 未導入）")

    print("\n=== Mermaid ===")
    print(generate_mermaid(problem))
    path = os.path.join(out, "plant4.html")
    export_mermaid(problem, path, title="Syngas Plant (2-stage kinetic reactors)",
                   style="diamond")
    print(f"\nフロー図(HTML): {path}")


if __name__ == "__main__":
    main()
