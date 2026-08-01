"""サンプル: example_plant2 の反応器を KineticReactor（速度論 PFR）に置き換えた版。

plant2 との違い:
  - 反応器: 固定転化率・固定選択率の Reactor（R1a/R1b）→ **速度論 PFR**（KineticReactor）
      R1 = ハイブリッド床 {synthesis: Cu/ZnO + dehydration: ZSM-5}
      R2 = カルボニル化床 {carbonylation: H-MOR}
    転化率は「与える」ものではなく、触媒質量・温度・圧力・組成から**決まる**。
  - 供給: 指定組成の Steam1（H2/CO=1.384, CO2 16.1%, **不活性 CH4 2.3%** を含む）。
  - 成分: 反応網7種（CO, CO2, H2, H2O, CH3OH, CH3OCH3, CH3COOCH3）＋ CH4, CH3CHO,
    CH3COOH, N2。後ろ3種は供給0・どの反応でも生成しないので全ストリームで0になるが、
    将来の拡張用に成分リストへ残してある。
    ※ 酢酸メチルの式は plant2 の "CH3OCOCH3" ではなく **"CH3COOCH3"**。
      KineticReactor の species_alias（CH3COOCH3 → reaction_rate の "MA"）に合わせる必要がある。

触媒量の決め方（**反応器入口基準の SV 5000/h**）:
    reaction_rate の例（fig_tandem_sty_stack_zsm5.py）と同じ定義
        N[mol/s] = SV/3600 · V_cat[m³] / 22.414e-3      →  V_cat = N[mol/h]·22.414e-3/SV
    循環があるため反応器入口流量は解いてみないと分からない。そこで
        「V_cat を仮定 → フローシートを解く → 入口流量から V_cat を再計算」
    を収束するまで繰り返す（外側の固定点反復・減衰0.5）。

要 reaction_rate: pip install -e ../reaction_rate
実行: PYTHONPATH=.:../reaction_rate/src python3 examples/example_plant3.py
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

# --- 成分（反応網7種 + 不活性・微量成分）---
C = ["H2", "CO", "CO2", "CH4", "H2O", "CH3OH", "CH3OCH3", "CH3COOCH3",
     "CH3CHO", "CH3COOH", "N2"]

# --- 新鮮供給 Steam1 [mol/h]（指定値）---
STEAM1 = {"H2": 4.911989, "CO": 3.549828, "CO2": 1.675074, "CH4": 0.243186,
          "H2O": 0.0288, "CH3CHO": 0.0, "CH3COOH": 0.0, "N2": 0.0}

# --- 触媒の設計基準（reaction_rate の fig_tandem_sty_stack_zsm5 と同一）---
SV = 5000.0                                    # [1/h] ← 反応器入口（循環込み）基準
RHO_SYN, RHO_DEH, RHO_MA = 1300.0, 800.0, 837.0  # 充填密度 [kg/m³]
F_SYN = 1.0 / 1.9          # ハイブリッド床の体積分率（synthesis 側）
R_MA = 2.0                 # V_MA / V_hybrid（1 mL : 2 mL）
VM_STP = 22.414e-3         # 標準状態モル体積 [m³/mol]
PURGE = 0.05               # 凝縮器ガスのパージ率

gas = StreamCondition(T=250, P="5MPaG", phase="gas")   # 5 MPaG + atm = 51.0 bar
liq = StreamCondition(T=40, P="5MPaG", phase="liquid")

NAMES = ["1. Steam1", "2. ReactorIn", "3. HybridOut", "4. ReactorOut", "5. CondGas",
         "6. CondLiq", "7. Purge", "8. RecycleGas", "9. MethylAcetate", "10. Sep1Liq",
         "11. RecycleMeOH", "12. Water"]


def catalyst_masses(v_tot: float) -> tuple[float, float, float]:
    """全触媒体積 [m³] → (m_synthesis, m_dehydration, m_carbonylation) [kg]。"""
    v_hyb = v_tot / (1.0 + R_MA)
    v_ma = v_tot - v_hyb
    return RHO_SYN * v_hyb * F_SYN, RHO_DEH * v_hyb * (1 - F_SYN), RHO_MA * v_ma


def volume_for(n_molh: float) -> float:
    """入口モル流量 [mol/h] → SV を満たす全触媒体積 [m³]。"""
    return n_molh * VM_STP / SV


def build(v_tot: float, seed: dict[str, np.ndarray] | None = None,
          feed: dict[str, float] | None = None):
    """全触媒体積 v_tot [m³] のフローシートを組む。

    seed は前回解（初期推定）。feed を渡すと新鮮供給の組成を差し替えられる
    （前工程 pattern1 の DryGas をそのまま流し込むため。既定は STEAM1）。
    """

    def S(i: int, cond: StreamCondition = gas) -> Stream:
        nm = NAMES[i - 1]
        return Stream(C, name=nm, order=i, condition=cond,
                      guess=seed.get(nm) if seed else None)

    Steam1      = Stream(C, name=NAMES[0], order=1, condition=gas,
                         flows=feed if feed is not None else STEAM1)
    ReactorIn   = S(2)
    HybridOut   = S(3)          # ハイブリッド床（合成＋脱水）の出口
    ReactorOut  = S(4)          # カルボニル化床の出口
    CondGas     = S(5)
    CondLiq     = S(6, liq)
    Purge       = S(7)
    RecycleGas  = S(8)
    MA_Product  = S(9, liq)
    Sep1Liq     = S(10, liq)
    RecycleMeOH = S(11)         # 分離塔2 ガス = 全量リサイクル（パージ無し）
    Water_Out   = S(12, liq)

    m_syn, m_deh, m_ma = catalyst_masses(v_tot)

    M1   = Mixer([Steam1, RecycleGas, RecycleMeOH], ReactorIn, name="M1")
    R1   = KineticReactor(inlet=ReactorIn, outlet=HybridOut,
                          masses={"synthesis": m_syn, "dehydration": m_deh},
                          models={"synthesis": "KOGAS", "dehydration": "ZSM5"},
                          k_eq3="thermo", n_points=2, name="R1 hybrid")
    R2   = KineticReactor(inlet=HybridOut, outlet=ReactorOut,
                          masses={"carbonylation": m_ma},
                          models={"carbonylation": "DTU-Cheung2007-2"},
                          n_points=2, name="R2 MA bed")
    Cond = Separator(ReactorOut, [CondGas, CondLiq], name="Condenser")
    SPg  = Splitter(CondGas, [Purge, RecycleGas], ratios=[PURGE, 1 - PURGE], name="SP-gas")
    Col1 = Separator(CondLiq, [MA_Product, Sep1Liq], name="Column1")
    Col2 = Separator(Sep1Liq, [RecycleMeOH, Water_Out], name="Column2")

    streams = [Steam1, ReactorIn, HybridOut, ReactorOut, CondGas, CondLiq, Purge,
               RecycleGas, MA_Product, Sep1Liq, RecycleMeOH, Water_Out]
    problem = Problem(streams=streams, units=[M1, R1, R2, Cond, SPg, Col1, Col2],
                      name="Syngas plant (kinetic tandem reactor)")

    # --- 分離指定（各分離器の液側への回収率）---
    # 凝縮器: 軽質ガス（H2/CO/CO2/CH4/N2/DME）は気相、含酸素の重質分は液相へ。
    problem.constrain_recovery(ReactorOut, CondLiq, {
        "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "N2": 0, "CH3OCH3": 0,
        "H2O": 1, "CH3OH": 1, "CH3COOCH3": 1, "CH3CHO": 1, "CH3COOH": 1,
    })
    # 塔1: 酢酸メチル（＋より軽いアセトアルデヒド）を留出、H2O/MeOH/酢酸は缶出。
    problem.constrain_recovery(CondLiq, Sep1Liq, {
        "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "N2": 0, "CH3OCH3": 0,
        "CH3COOCH3": 0, "CH3CHO": 0, "H2O": 1, "CH3OH": 1, "CH3COOH": 1,
    })
    # 塔2: メタノールを留出（全量リサイクル）、水と酢酸は缶出。
    problem.constrain_recovery(Sep1Liq, Water_Out, {
        "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "N2": 0, "CH3OCH3": 0,
        "CH3COOCH3": 0, "CH3CHO": 0, "CH3OH": 0, "H2O": 1, "CH3COOH": 1,
    })
    return problem, streams


# 反応器入口 / 新鮮供給 の初期推定（循環ぶん）。この供給・パージ率での収束値は約 6.5。
# 外側反復は割線法なので外れても収束するが、近いほど反復回数が減る。
RECYCLE_GUESS = 6.5


def solve_with_sv(max_outer: int = 8, tol: float = 1e-4, verbose: bool = True,
                  feed: dict[str, float] | None = None):
    """SV 一定（反応器入口基準）になるまで触媒体積 V を外側反復して解く。

    g(V) = volume_for(反応器入口流量(V)) − V の零点を**割線法**で探す。
    単純な減衰付き固定点反復だと 1 回の求解が数分かかるのに 6〜8 回必要になるため。
    各求解は前回解をそのまま初期推定に使う（ウォームスタート）。
    feed で新鮮供給の組成を差し替えられる（既定は STEAM1）。
    """
    feed = feed if feed is not None else STEAM1
    state: dict[str, object] = {}

    def g(v: float, it: int) -> float:
        t0 = time.time()
        problem, streams = build(v, state.get("seed"), feed=feed)
        if it == 1 and verbose:
            print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
        sol = problem.solve(bounds=(0, np.inf), tol=1e-6,
                            ftol=1e-12, xtol=1e-12, gtol=1e-12)
        if not sol.success:
            raise RuntimeError(f"外側反復 {it} で収束せず: {sol}")
        n_in = float(streams[1].total_flow.eval())
        v_new = volume_for(n_in)
        state["seed"] = {s.name: s.molar_flows.copy() for s in streams}
        state["problem"], state["streams"], state["v"] = problem, streams, v
        n_pfr = sum(u.n_pfr_calls for u in problem.units if hasattr(u, "n_pfr_calls"))
        if verbose:
            print(f"  [{it}] V_cat={v*1e6:8.4f} mL → 反応器入口 {n_in:8.4f} mol/h "
                  f"→ V_cat(SV={SV:.0f}/h)={v_new*1e6:8.4f} mL   "
                  f"({time.time()-t0:.1f}s, PFR積分 {n_pfr}回, nfev={sol.nfev})")
        return v_new - v

    v_prev = volume_for(RECYCLE_GUESS * sum(feed.values()))
    g_prev = g(v_prev, 1)
    if abs(g_prev) <= tol * v_prev:
        return state["problem"], state["streams"], v_prev

    v = v_prev + g_prev                       # 1 歩目は固定点ステップ
    for it in range(2, max_outer + 1):
        g_cur = g(v, it)
        if abs(g_cur) <= tol * v:
            return state["problem"], state["streams"], v
        denom = g_cur - g_prev
        if denom == 0.0:
            raise RuntimeError("割線法の分母が 0 になりました")
        v_next = v - g_cur * (v - v_prev) / denom
        v_prev, g_prev = v, g_cur
        v = float(np.clip(v_next, 0.2 * v, 5.0 * v))   # 暴れ止め
    raise RuntimeError("触媒体積の外側反復が収束しませんでした")


def main():
    print(f"=== plant3: 速度論タンデム反応器（250℃ / 51.0 bar / SV {SV:.0f}/h・反応器入口基準）===")
    problem, streams, v_tot = solve_with_sv()
    (Steam1, ReactorIn, HybridOut, ReactorOut, CondGas, CondLiq, Purge,
     RecycleGas, MA_Product, Sep1Liq, RecycleMeOH, Water_Out) = streams

    m_syn, m_deh, m_ma = catalyst_masses(v_tot)
    print(f"\n--- 触媒（収束値）---")
    print(f"全触媒体積      : {v_tot*1e6:.4f} mL")
    print(f"synthesis   (Cu/ZnO) : {m_syn*1e3:.4f} g")
    print(f"dehydration (ZSM-5)  : {m_deh*1e3:.4f} g")
    print(f"carbonylation(H-MOR) : {m_ma*1e3:.4f} g")

    print()
    print(stream_table(problem.streams, basis=["mol", "mole_frac"]))

    # --- 性能指標 ---
    co_in, co_r_in, co_r_out = (Steam1.flow_of("CO"), ReactorIn.flow_of("CO"),
                                ReactorOut.flow_of("CO"))
    ma = MA_Product.flow_of("CH3COOCH3")
    print("\n--- 性能 ---")
    print(f"CO 1パス転化率      : {(co_r_in - co_r_out)/co_r_in*100:.2f} %")
    print(f"酢酸メチル生成       : {ma:.4f} mol/h（新鮮 CO 基準の炭素収率 "
          f"{ma*3/co_in*100:.2f} %）")
    print(f"CO パージ損失        : {Purge.flow_of('CO'):.4f} mol/h")

    # --- 不活性 CH4 の蓄積（5% パージでどこまで溜まるか）---
    print("\n--- 不活性 CH4 の蓄積（パージ率 {:.0%}）---".format(PURGE))
    print(f"供給 CH4            : {Steam1.flow_of('CH4'):.4f} mol/h")
    print(f"循環ガス中 CH4      : {RecycleGas.flow_of('CH4'):.4f} mol/h "
          f"（{RecycleGas.flow_of('CH4')/float(RecycleGas.total_flow.eval())*100:.2f} mol%）")
    print(f"パージ CH4          : {Purge.flow_of('CH4'):.4f} mol/h（供給と一致すれば収支OK）")
    print(f"反応器入口 CH4 濃度 : "
          f"{ReactorIn.flow_of('CH4')/float(ReactorIn.total_flow.eval())*100:.2f} mol%")

    # --- メタノール収支（plant2 と同じ確認）---
    meoh_made = HybridOut.flow_of("CH3OH") - ReactorIn.flow_of("CH3OH")
    print("\n--- メタノール収支 ---")
    print(f"ハイブリッド床 正味生成: {meoh_made:.4f} mol/h")
    print(f"循環 MeOH             : {RecycleMeOH.flow_of('CH3OH'):.4f} mol/h（パージ無し）")

    # --- ファイル出力 ---
    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)

    csv_path = os.path.join(out, "plant3_streams.csv")
    to_csv(problem.streams, csv_path,
           basis=["mol", "mole_frac", "mass", "normal_volume"])
    print(f"\nCSV 出力  : {csv_path}")
    try:
        xlsx = os.path.join(out, "plant3_streams.xlsx")
        to_excel(problem.streams, xlsx,
                 basis=["mol", "mole_frac", "mass", "normal_volume"])
        print(f"Excel 出力: {xlsx}")
    except ImportError:
        print("Excel 出力スキップ（openpyxl 未導入）")

    # --- Mermaid 出力 ---
    print("\n=== Mermaid ===")
    print(generate_mermaid(problem))
    path = os.path.join(out, "plant3.html")
    export_mermaid(problem, path, title="Syngas Plant (kinetic tandem reactor)",
                   style="diamond")
    print(f"\nフロー図(HTML): {path}")


if __name__ == "__main__":
    main()
