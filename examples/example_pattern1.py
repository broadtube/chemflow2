"""サンプル: ChemFlow(旧) の pattern1（改質器 + 水凝縮）を chemflow2 で再現する。

plant3 / plant4 の新鮮供給 Steam1 は、この pattern1 の DryGas（Gibbs 平衡後に水を
凝縮分離したガス）そのもの。前工程を chemflow2 側に持ってくることで、**改質条件を
振ったときに下流のプラントがどう動くかを一気通貫で見られる**ようにする。

    RG_feed(H2,CH4) ┐
    CO2_feed        ├→ M1 → Mixed(850℃) → G1(Gibbs) → ReactOut → Cond ┬ DryGas    (25℃)
    H2O_feed        ┘                                                  └ Condensate(25℃)

構成:
  - Mixed は組成（vol%）と全流量を固定し、**各 Feed の流量は逆算**する（Mixer が逆算式）。
  - G1 は Cantera のギブズ自由エネルギー最小化。平衡種は CO2/CH4/H2O/CO/H2 の5種。
  - 凝縮は Separator（収支のみ）+ 物性制約で閉じる:
        水     → constrain_saturation（Antoine 式の飽和）      1 本
        溶存ガス → constrain_henry（Henry 則, Sander 2023）    4 本
    シャープスプリット（回収率 0/1）では**溶存ガスぶんがずれる**ため、DryGas を厳密に
    再現するには Henry 則が要る（CO2 で 0.16% の差になる）。

要 Cantera: pip install chemflow2[gibbs]
実行: PYTHONPATH=. python3 examples/example_pattern1.py

⚠️ N2_feed は Mixer に含めない。含めると N2 の収支式（0 = 0）が 1 本増えて過剰決定に
   なる（旧 ChemFlow 側も同じ理由で除外している）。
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
    to_csv,
    to_excel,
)

# --- 条件 ---
T_REACT, T_COND = 850.0, 25.0
PRESSURE = "1.04MPaG"
SPECIES = ["CO2", "CH4", "H2O", "CO", "H2"]      # Gibbs 平衡種
DISSOLVED = ["H2", "CO", "CO2", "CH4"]           # 液水に溶けるガス（H2O 以外）

# Mixed の指定: 全流量 204.72 NL/h、組成 vol%
MIXED_TOTAL_NL = 204.72
MIXED_VOL_FRAC = {"H2": 0.097, "CO2": 0.365, "CH4": 0.234, "H2O": 0.304}

hot = StreamCondition(T=T_REACT, P=PRESSURE, phase="gas")
cold_gas = StreamCondition(T=T_COND, P=PRESSURE, phase="gas")
cold_liq = StreamCondition(T=T_COND, P=PRESSURE, phase="liquid")
feed_cond = StreamCondition(T=25, P=PRESSURE, phase="gas")

# --- ストリーム ---
# 各 Feed は含む成分だけを持つ未知ストリーム（流量は Mixer が逆算する）
RG_feed  = Stream(["H2", "CH4"], name="1. RG_feed", order=1, condition=feed_cond)
CO2_feed = Stream(["CO2"],       name="2. CO2_feed", order=2, condition=feed_cond)
H2O_feed = Stream(["H2O"],       name="3. H2O_feed", order=3, condition=feed_cond)

Mixed = Stream(list(MIXED_VOL_FRAC), name="4. Mixed", order=4, condition=hot,
               flows=MIXED_VOL_FRAC, basis="volume_frac", total=MIXED_TOTAL_NL)

ReactOut   = Stream(SPECIES, name="5. ReactOut", order=5, condition=hot)
DryGas     = Stream(SPECIES, name="6. DryGas", order=6, condition=cold_gas)
Condensate = Stream(SPECIES, name="7. Condensate", order=7, condition=cold_liq)

# --- 装置 ---
M1   = Mixer([RG_feed, CO2_feed, H2O_feed], Mixed, name="M1")
G1   = GibbsReactor(inlet=Mixed, outlet=ReactOut, species=SPECIES,
                    T=T_REACT, P=PRESSURE, name="G1")
Cond = Separator(ReactOut, [DryGas, Condensate], name="Condenser")

problem = Problem(streams=[RG_feed, CO2_feed, H2O_feed, Mixed,
                           ReactOut, DryGas, Condensate],
                  units=[M1, G1, Cond],
                  name="Pattern1 (reformer + water knockout)")

# --- 凝縮の分配を物性で閉じる ---
problem.constrain_saturation(DryGas, "H2O", T=T_COND, P=PRESSURE)          # 1 本
problem.constrain_henry(DryGas, Condensate, DISSOLVED,
                        T=T_COND, P=PRESSURE)                              # 4 本

# --- 旧 ChemFlow の pattern1_result.csv の値（照合用）---
REFERENCE = {
    "5. ReactOut":   {"H2": 4.912179, "CO": 3.549996, "CO2": 1.677814,
                      "CH4": 0.243203, "H2O": 2.538494},
    "7. Condensate": {"H2": 0.000190, "CO": 0.000167, "CO2": 0.002740,
                      "CH4": 0.000017, "H2O": 2.509694},
    "6. DryGas":     {"H2": 4.911989, "CO": 3.549828, "CO2": 1.675074,
                      "CH4": 0.243186, "H2O": 0.028800},
}


def main():
    print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
    # 非負制約は付けない。溶存量が 1e-5 mol/h オーダーと小さく、bounds 付きの
    # least_squares だと境界へ漸近するだけで gtol 停止して残差が残る。
    sol = problem.solve()
    print(sol)
    print()
    print(stream_table(problem.streams, basis=["mol", "mole_frac"]))

    print("\n--- 逆算された Feed [mol/h] ---")
    for s in (RG_feed, CO2_feed, H2O_feed):
        inner = "  ".join(f"{f} {s.flow_of(f):.6f}" for f in s.formulas)
        print(f"{s.name:14s} {inner}   (計 {float(s.total_flow.eval()):.6f})")

    print(f"\n--- 旧 ChemFlow pattern1 との照合（許容差 1e-5 mol/h）---")
    worst = 0.0
    for stream in (ReactOut, Condensate, DryGas):
        ref = REFERENCE[stream.name]
        diffs = {f: stream.flow_of(f) - v for f, v in ref.items()}
        worst = max(worst, max(abs(d) for d in diffs.values()))
        cells = "  ".join(f"{f}:{d:+.2e}" for f, d in diffs.items())
        print(f"{stream.name:14s} {cells}")
    print(f"\n最大差 {worst:.2e} mol/h → {'一致' if worst < 1e-5 else '不一致'}"
          "（旧側 CSV は小数6桁なので 5e-7 程度の丸め差は残る）")

    # --- 下流プラントの供給になるか確認 ---
    print("\n--- plant3 / plant4 の Steam1 との対応 ---")
    steam1 = {"H2": 4.911989, "CO": 3.549828, "CO2": 1.675074,
              "CH4": 0.243186, "H2O": 0.028800}
    for f, v in steam1.items():
        print(f"  {f:5s} DryGas {DryGas.flow_of(f):10.6f}  Steam1 {v:10.6f}")
    print(f"  {'計':5s} DryGas {float(DryGas.total_flow.eval()):10.6f}  "
          f"Steam1 {sum(steam1.values()):10.6f}")

    # --- 出力 ---
    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)
    to_csv(problem.streams, os.path.join(out, "pattern1_streams.csv"),
           basis=["mol", "mole_frac", "mass", "normal_volume"])
    print(f"\nCSV 出力  : {os.path.join(out, 'pattern1_streams.csv')}")
    try:
        to_excel(problem.streams, os.path.join(out, "pattern1_streams.xlsx"),
                 basis=["mol", "mole_frac", "mass", "normal_volume"])
        print(f"Excel 出力: {os.path.join(out, 'pattern1_streams.xlsx')}")
    except ImportError:
        print("Excel 出力スキップ（openpyxl 未導入）")

    path = os.path.join(out, "pattern1.html")
    export_mermaid(problem, path, title="Pattern1 (reformer + water knockout)",
                   style="diamond")
    print(f"フロー図(HTML): {path}")


if __name__ == "__main__":
    main()
