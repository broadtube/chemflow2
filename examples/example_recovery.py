"""サンプル: 成分回収率で分離を指定する（constrain_recovery）。

Separator は物質収支だけを課すノード。どの成分がどの出口へ行くかは
constrain_recovery(inlet, outlet, {成分: 回収率}) で指定する。

    Feed → Reactor(タンデム反応) → Separator ─┬─→ Gas
                                              └─→ Liquid（H2O を全量回収）
"""

from chemflow2 import Problem, Reaction, Reactor, Separator, Stream, StreamCondition

components = ["H2", "CO", "CO2", "CH3OCH3", "H2O"]

S1 = Stream(components, name="1. Feed",     order=1, condition=StreamCondition(T=300, P=1, phase="gas"))
S3 = Stream(components, name="3. ReactOut", order=3)
S4 = Stream(components, name="4. Gas",      order=4)
S5 = Stream(components, name="5. Liquid",   order=5)

rxn1 = Reaction(stoich={"CO": -2, "H2": -4, "CH3OCH3": 1, "H2O": 1}, name="Tandem")
R1   = Reactor(inlet=S1, outlet=S3, reactions=[rxn1],
               key_component="CO", conversion=0.8, selectivities=[1], name="R1")
Sep1 = Separator(inlet=S3, outlet=[S4, S5], name="Sep1")

problem = Problem(streams=[S1, S3, S4, S5], units=[R1, Sep1], name="DME + Recovery")

# Feed 指定（組成 + 総流量）
problem.constrain_fracs(S1, {"H2": 0.48, "CO": 0.24})
problem.constrain(S1.frac_of("CH3OCH3"), 0.0)
problem.constrain(S1.frac_of("H2O"), 0.0)
problem.constrain(S1.total_flow, 165, name="S1 total flow")

# 分離: Liquid への回収率で指定（H2O は全量 Liquid、他は全量 Gas = Liquid 回収率 0）
problem.constrain_recovery(S3, S5, {
    "H2O": 1.0,       # H2O は 100% Liquid へ
    "H2": 0.0,
    "CO": 0.0,
    "CO2": 0.0,
    "CH3OCH3": 0.0,
})

print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
sol = problem.solve()
print(sol)
print()
sol.print_report()

print(f"\nH2O 回収率 (Liquid) = {S5.flow_of('H2O') / S3.flow_of('H2O'):.1%}")
print(f"CO 転化率           = {1 - S3.flow_of('CO') / S1.flow_of('CO'):.1%}")
