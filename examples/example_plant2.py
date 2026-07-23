"""サンプル: example_plant に「反応器内のメタノール脱水」を追加した版。

追加反応:  2 CH3OH -> CH3OCH3 + H2O   （メタノール脱水で DME）

この反応は CO を含まないので CO-key の反応器には入れられない。そこで
**CH3OH を key にした反応器 R1b を直列に 1 段追加**する（R1a: CO 系、R1b: MeOH 脱水）。

効果:
    example_plant.py では、メタノールが循環ループ内で生成されるだけで消費されず、
    100% 循環だと無限に溜まる（定常解が特異）。そのため 5% パージが必須だった。
    本例では R1b がメタノールを DME に変換するため、メタノールに反応器内シンクができ、
    **MeOH ガスをパージ無しで 100% 循環しても系が閉じる**（DME は凝縮器ガスの
    5% パージから抜ける）。元の指定どおり「分離塔2ガスは全量リサイクル」を実現できる。
"""

import os

import numpy as np

from chemflow2 import (
    Mixer,
    Problem,
    Reaction,
    Reactor,
    Separator,
    Splitter,
    Stream,
    StreamCondition,
    export_mermaid,
    generate_mermaid,
    stream_table,
)

C = ["H2", "CO", "CO2", "CH4", "H2O", "CH3OH", "CH3OCH3", "CH3OCOCH3"]
gas = StreamCondition(T=250, P="5MPaG", phase="gas")
liq = StreamCondition(T=40, P="5MPaG", phase="liquid")


def S(name, order, cond=gas, **kw):
    return Stream(C, name=name, order=order, condition=cond, **kw)


# --- ストリーム ---
FreshFeed   = S("1. FreshFeed", 1, flows={"H2": 150, "CO": 100})   # H2/CO = 1.5
ReactorIn   = S("2. ReactorIn", 2)
R1a_out     = S("3. R1a_out", 3)          # CO 系反応の後、MeOH 脱水の前
ReactorOut  = S("4. ReactorOut", 4)       # R1b（MeOH 脱水）の後
CondGas     = S("5. CondGas", 5)
CondLiq     = S("6. CondLiq", 6, cond=liq)
Purge       = S("7. Purge", 7)
RecycleGas  = S("8. RecycleGas", 8)
MA_Product  = S("9. MethylAcetate", 9, cond=liq)
Sep1Liq     = S("10. Sep1Liq", 10, cond=liq)
RecycleMeOH = S("11. RecycleMeOH", 11)    # 分離塔2 ガス = 全量リサイクル（パージ無し）
Water_Out   = S("12. Water", 12, cond=liq)

# --- 反応 ---
co_rxns = [
    Reaction({"CO": -1, "H2": -2, "CH3OH": 1}, name="MeOH"),
    Reaction({"CO": -2, "H2": -4, "CH3OCH3": 1, "H2O": 1}, name="DME"),
    Reaction({"CO": -3, "H2": -4, "CH3OCOCH3": 1, "H2O": 1}, name="MethylAcetate"),
    Reaction({"CO": -1, "H2O": -1, "CO2": 1, "H2": 1}, name="WGS"),
]
meoh_dehydration = Reaction({"CH3OH": -2, "CH3OCH3": 1, "H2O": 1}, name="MeOH->DME")

# --- 装置 ---
M1   = Mixer([FreshFeed, RecycleGas, RecycleMeOH], ReactorIn, name="M1")
R1a  = Reactor(ReactorIn, R1a_out, co_rxns, key_component="CO",
               conversion=0.30, selectivities=[0.4, 0.3, 0.2, 0.1], name="R1a")
R1b  = Reactor(R1a_out, ReactorOut, [meoh_dehydration], key_component="CH3OH",
               conversion=0.50, selectivities=[1], name="R1b")
Cond = Separator(ReactorOut, [CondGas, CondLiq], name="Condenser")
SPg  = Splitter(CondGas, [Purge, RecycleGas], ratios=[0.05, 0.95], name="SP-gas")
Col1 = Separator(CondLiq, [MA_Product, Sep1Liq], name="Column1")
Col2 = Separator(Sep1Liq, [RecycleMeOH, Water_Out], name="Column2")   # ガス=MeOH を全量循環

problem = Problem(
    streams=[FreshFeed, ReactorIn, R1a_out, ReactorOut, CondGas, CondLiq, Purge,
             RecycleGas, MA_Product, Sep1Liq, RecycleMeOH, Water_Out],
    units=[M1, R1a, R1b, Cond, SPg, Col1, Col2],
    name="Syngas plant (with MeOH dehydration)",
)

# --- 分離指定（各分離器の液側への回収率）---
problem.constrain_recovery(ReactorOut, CondLiq, {
    "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "CH3OCH3": 0,
    "H2O": 1, "CH3OH": 1, "CH3OCOCH3": 1,
})
problem.constrain_recovery(CondLiq, Sep1Liq, {
    "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "CH3OCH3": 0, "CH3OCOCH3": 0,
    "H2O": 1, "CH3OH": 1,
})
problem.constrain_recovery(Sep1Liq, Water_Out, {
    "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "CH3OCH3": 0, "CH3OCOCH3": 0, "CH3OH": 0,
    "H2O": 1,
})

print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
sol = problem.solve(bounds=(0, np.inf), ftol=1e-15, xtol=1e-15, gtol=1e-15, max_nfev=50000)
print(sol)
print()
print(stream_table(problem.streams, basis=["mol", "mole_frac"]))

# メタノール収支の確認
meoh_made = R1a_out.flow_of("CH3OH") - ReactorIn.flow_of("CH3OH")   # R1a での正味生成
meoh_used = R1a_out.flow_of("CH3OH") - ReactorOut.flow_of("CH3OH")  # R1b での消費
print("\n--- メタノール収支 ---")
print(f"R1a 生成: {meoh_made:.3f} mol/h")
print(f"R1b 消費: {meoh_used:.3f} mol/h")
print(f"循環 MeOH: {RecycleMeOH.flow_of('CH3OH'):.3f} mol/h（パージ無しでも有限に収束）")

# --- Mermaid 出力 ---
print("\n=== Mermaid ===")
print(generate_mermaid(problem))
out = os.path.join(os.path.dirname(__file__), "output")
os.makedirs(out, exist_ok=True)
path = os.path.join(out, "plant2.html")
export_mermaid(problem, path, title="Syngas Plant + MeOH dehydration", style="diamond")
print(f"\nフロー図(HTML): {path}")
