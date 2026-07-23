"""サンプル: 合成ガス（H2/CO）からの多段プロセス（循環 2 系統）。

フロー:
    FreshFeed(H2/CO=1.5) ─┐
    RecycleGas(95%) ──────┤→ Mixer → Reactor → Condenser ─┬─ Gas → Splitter ─┬─ Purge(5%) Out
    RecycleMeOH(95%) ─────┘                               │                  └─ Recycle(95%) ↩
                                                          └─ Liq → 分離塔1 ─┬─ Gas: 酢酸メチル Out
                                                                            └─ Liq → 分離塔2 ─┬─ Gas: MeOH → Splitter ─┬─ Purge(5%) Out
                                                                                              └─ Liq: H2O Out         └─ Recycle(95%) ↩

反応（すべて CO を key、原子収支は Reaction が自動検証）:
    1. CO + 2H2 → CH3OH                       (メタノール)
    2. 2CO + 4H2 → CH3OCH3 + H2O              (DME)
    3. 3CO + 4H2 → CH3OCOCH3 + H2O            (酢酸メチル)
    4. CO + H2O → CO2 + H2                    (WGS で CO2)

分離:
    凝縮器  … H2,CO,CO2,CH3OCH3 → ガス / H2O,CH3OH,CH3OCOCH3 → 液
    分離塔1 … CH3OCOCH3 → ガス(製品) / H2O,CH3OH → 液
    分離塔2 … CH3OH → ガス(循環) / H2O → 液(製品)

注: 循環ループ内で「生成・消費される成分」は必ずパージが要る（無いと成分が無限に
溜まり定常解が存在しない/特異になる）。そのため凝縮器ガスと MeOH ガスの双方に 5% パージを置く。
"""

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
ReactorOut  = S("3. ReactorOut", 3)
CondGas     = S("4. CondGas", 4)
CondLiq     = S("5. CondLiq", 5, cond=liq)
Purge       = S("6. Purge", 6)
RecycleGas  = S("7. RecycleGas", 7)
MA_Product  = S("8. MethylAcetate", 8, cond=liq)      # 分離塔1 ガス（製品）
Sep1Liq     = S("9. Sep1Liq", 9, cond=liq)
MeOH_Gas    = S("10. MeOH_Gas", 10)
Water_Out   = S("11. Water", 11, cond=liq)            # 分離塔2 液（製品）
MeOH_Purge  = S("12. MeOH_Purge", 12)
RecycleMeOH = S("13. RecycleMeOH", 13)

# --- 反応 ---
rxns = [
    Reaction({"CO": -1, "H2": -2, "CH3OH": 1}, name="MeOH"),
    Reaction({"CO": -2, "H2": -4, "CH3OCH3": 1, "H2O": 1}, name="DME"),
    Reaction({"CO": -3, "H2": -4, "CH3OCOCH3": 1, "H2O": 1}, name="MethylAcetate"),
    Reaction({"CO": -1, "H2O": -1, "CO2": 1, "H2": 1}, name="WGS"),
]

# --- 装置 ---
M1  = Mixer([FreshFeed, RecycleGas, RecycleMeOH], ReactorIn, name="M1")
R1  = Reactor(ReactorIn, ReactorOut, rxns, key_component="CO",
              conversion=0.30, selectivities=[0.4, 0.3, 0.2, 0.1], name="R1")
Cond = Separator(ReactorOut, [CondGas, CondLiq], name="Condenser")
SP1 = Splitter(CondGas, [Purge, RecycleGas], ratios=[0.05, 0.95], name="SP-gas")
Sep1 = Separator(CondLiq, [MA_Product, Sep1Liq], name="Column1")
Sep2 = Separator(Sep1Liq, [MeOH_Gas, Water_Out], name="Column2")
SP2 = Splitter(MeOH_Gas, [MeOH_Purge, RecycleMeOH], ratios=[0.05, 0.95], name="SP-MeOH")

problem = Problem(
    streams=[FreshFeed, ReactorIn, ReactorOut, CondGas, CondLiq, Purge, RecycleGas,
             MA_Product, Sep1Liq, MeOH_Gas, Water_Out, MeOH_Purge, RecycleMeOH],
    units=[M1, R1, Cond, SP1, Sep1, Sep2, SP2],
    name="Syngas plant",
)

# --- 分離指定（各分離器の液側への回収率）---
problem.constrain_recovery(ReactorOut, CondLiq, {
    "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "CH3OCH3": 0,     # ガスへ
    "H2O": 1, "CH3OH": 1, "CH3OCOCH3": 1,                    # 液へ
})
problem.constrain_recovery(CondLiq, Sep1Liq, {
    "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "CH3OCH3": 0, "CH3OCOCH3": 0,  # ガスへ
    "H2O": 1, "CH3OH": 1,                                                 # 液へ
})
problem.constrain_recovery(Sep1Liq, Water_Out, {
    "H2": 0, "CO": 0, "CO2": 0, "CH4": 0, "CH3OCH3": 0, "CH3OCOCH3": 0, "CH3OH": 0,  # ガスへ
    "H2O": 1,                                                                          # 液へ
})

print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
# メタノールが 95% 循環するためループゲインが大きく硬い系。許容値を締める。
sol = problem.solve(bounds=(0, np.inf), ftol=1e-15, xtol=1e-15, gtol=1e-15, max_nfev=50000)
print(sol)
print()
print(stream_table(problem.streams, basis=["mol", "mole_frac"]))

# --- Mermaid 出力 ---
print("\n=== Mermaid ===")
print(generate_mermaid(problem))
import os
out = os.path.join(os.path.dirname(__file__), "output")
os.makedirs(out, exist_ok=True)
path = os.path.join(out, "plant.html")
export_mermaid(problem, path, title="Syngas Plant", style="diamond")
print(f"\nフロー図(HTML): {path}")
