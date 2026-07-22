"""サンプル: Splitter で循環系を直截に書く。

フロー:
    Feed(N2O4=10) → Mixer → Reactor(N2O4->2NO2, 50%) → Splitter ─┬─→ Product (70%)
                     ↑                                             └─→ Recycle (30%) ↩

example_recycle.py では Separator + 生の Expr 制約
    problem.constrain(Expr(lambda: Recycle.molar_flows - 0.3 * Rout.molar_flows))
が必要だったが、Splitter を使うと分割比を ratios=[0.7, 0.3] と書くだけで済む。
"""

import numpy as np

from chemflow2 import (
    Mixer,
    Problem,
    Reaction,
    Reactor,
    Splitter,
    Stream,
    StreamCondition,
)

components = ["N2O4", "NO2"]
cond = StreamCondition(T=25, P=1, phase="gas")

Feed    = Stream(components, name="1. Feed",     order=1, flows={"N2O4": 10, "NO2": 0}, condition=cond)
Recycle = Stream(components, name="5. Recycle",  order=5, condition=cond)   # テア（未知）
Mixed   = Stream(components, name="2. Mixed",    order=2, condition=cond)
Rout    = Stream(components, name="3. ReactOut", order=3, condition=cond)
Product = Stream(components, name="4. Product",  order=4, condition=cond)

rxn = Reaction(stoich={"N2O4": -1, "NO2": 2}, name="N2O4->2NO2")

M1  = Mixer(inlet=[Feed, Recycle], outlet=Mixed, name="M1")
R1  = Reactor(inlet=Mixed, outlet=Rout, reactions=[rxn],
              key_component="N2O4", conversion=0.5, selectivities=[1], name="R1")
SP1 = Splitter(inlet=Rout, outlet=[Product, Recycle], ratios=[0.7, 0.3], name="SP1")

problem = Problem(
    streams=[Feed, Recycle, Mixed, Rout, Product],
    units=[M1, R1, SP1],
    name="Recycle via Splitter",
)

print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
sol = problem.solve(bounds=(0, np.inf))
print(sol)
print()
sol.print_report()

print("\n--- 検証 ---")
print(f"Recycle/ReactOut 比: {Recycle.total_flow / Rout.total_flow:.4f} (target 0.30)")
print(f"Product/ReactOut 比: {Product.total_flow / Rout.total_flow:.4f} (target 0.70)")
print(f"Feed N2O4 = {Feed.flow_of('N2O4'):.2f}, Product 総流量 = {Product.total_flow:.4f}")
