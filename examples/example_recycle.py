"""サンプル: リサイクルあり（宣言的スタイル）

フロー:
    Feed(N2O4=10) → Mixer → Reactor(N2O4->2NO2, 50%) → Separator ─┬─→ Product (70%)
                     ↑                                             └─→ Recycle (30%) ↩

循環（テア）ストリーム Recycle は「未知ストリームを 1 つ宣言し、
Mixer の入口と Separator の出口の両方で同じオブジェクトを共有する」だけで表現する。
（成分は molmass が解決できる実在の示性式を使う。ここでは N2O4 <-> 2 NO2 の解離。）
"""

import numpy as np

from chemflow2 import (
    Expr,
    Mixer,
    Problem,
    Reaction,
    Reactor,
    Separator,
    Stream,
    StreamCondition,
)

components = ["N2O4", "NO2"]
cond = StreamCondition(T=25, P=1, phase="gas")

# --- ストリーム宣言 ---
Feed    = Stream(components, name="1. Feed",     order=1, flows={"N2O4": 10, "NO2": 0}, condition=cond)
Recycle = Stream(components, name="5. Recycle",  order=5, condition=cond)   # テア（未知）
Mixed   = Stream(components, name="2. Mixed",    order=2, condition=cond)
Rout    = Stream(components, name="3. ReactOut", order=3, condition=cond)
Product = Stream(components, name="4. Product",  order=4, condition=cond)

# --- 反応 ---
rxn = Reaction(stoich={"N2O4": -1, "NO2": 2}, name="N2O4->2NO2")

# --- 装置 ---
M1   = Mixer(inlet=[Feed, Recycle], outlet=Mixed, name="M1")
R1   = Reactor(inlet=Mixed, outlet=Rout, reactions=[rxn],
               key_component="N2O4", conversion=0.5, selectivities=[1], name="R1")
Sep1 = Separator(inlet=Rout, outlet=[Product, Recycle], name="Sep1")

# --- 問題 ---
problem = Problem(
    streams=[Feed, Recycle, Mixed, Rout, Product],
    units=[M1, R1, Sep1],
    name="Recycle Example",
)

# 分離器の分配: Recycle は各成分の 30%（成分ごとのベクトル制約）
problem.constrain(
    Expr(lambda: Recycle.molar_flows - 0.3 * Rout.molar_flows),
    name="Recycle = 30% of ReactOut",
)

# --- 求解 ---
print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
sol = problem.solve(bounds=(0, np.inf))
print(sol)
print()
sol.print_report()

print("\n--- 全体物質収支 ---")
print(f"Feed N2O4:    {Feed.flow_of('N2O4'):.4f}")
print(f"Product N2O4: {Product.flow_of('N2O4'):.4f}")
print(f"Product NO2:  {Product.flow_of('NO2'):.4f}")
print(f"Product total flow: {Product.total_flow:.4f}")   # Expr を数値として表示
