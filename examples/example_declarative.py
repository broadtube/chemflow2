"""サンプル: あなたのイメージに沿った宣言的な書き方（循環なし）。

Feed(組成 + 総流量を制約で指定) → Reactor(タンデム反応) → Separator。
constrain_fracs と constrain(total_flow, ...) の使い方を示す。
"""

from chemflow2 import (
    Problem,
    Reaction,
    Reactor,
    Separator,
    Stream,
    StreamCondition,
)

components = ["H2", "CO", "CO2", "CH3OCH3", "H2O"]

S1 = Stream(components, name="1. Feed",     order=1, condition=StreamCondition(T=300, P=1, phase="gas"))
S3 = Stream(components, name="3. ReactOut", order=3, condition=StreamCondition(T=300, P=1, phase="gas"))
S4 = Stream(components, name="4. Gas",      order=4)
S5 = Stream(components, name="5. Liquid",   order=5)

rxn1 = Reaction(stoich={"CO": -2, "H2": -4, "CH3OCH3": 1, "H2O": 1}, name="Tandem")

R1 = Reactor(inlet=S1, outlet=S3, reactions=[rxn1],
             key_component="CO", conversion=0.8, selectivities=[1], name="R1")
Sep1 = Separator(inlet=S3, outlet=[S4, S5], components=components, name="Sep1")

problem = Problem(streams=[S1, S3, S4, S5], units=[R1, Sep1], name="DME Tandem")

# Feed 組成（H2, CO を指定 / CO2 は残り。CH3OCH3, H2O は 0）
problem.constrain_fracs(S1, {"H2": 0.48, "CO": 0.24})
problem.constrain(S1.frac_of("CH3OCH3"), 0.0)
problem.constrain(S1.frac_of("H2O"), 0.0)
# Feed 総流量
problem.constrain(S1.total_flow, 165, name="S1 total flow specified")

# 分離: H2O は全量 Liquid へ、それ以外は全量 Gas へ（成分ごとの完全分離）
# flow_expr で「その成分の流量」を直截にゼロ指定できる。
for f in ["H2", "CO", "CO2", "CH3OCH3"]:
    problem.constrain(S5.flow_expr(f), 0.0, name=f"{f} -> gas")
problem.constrain(S4.flow_expr("H2O"), 0.0, name="H2O -> liquid")

print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
# 循環なし・非負解になる素直な系なので root で厳密に解く（bounds 不要）
sol = problem.solve()
print(sol)
print()
sol.print_report()
print(f"\nS1 total flow = {S1.total_flow:.2f}  (target 165)")
print(f"CO 転化率     = {1 - S3.flow_of('CO') / S1.flow_of('CO'):.2%}")
