"""サンプル: Gibbs 平衡反応器（水蒸気メタン改質）。

CH4 + H2O を高温で平衡させると、改質反応
    CH4 + H2O <-> CO + 3H2
    CO  + H2O <-> CO2 + H2
が同時に平衡に達する。転化率を与えるのではなく、指定 T・P での
平衡組成そのものを Cantera が解く。

要 Cantera: pip install chemflow2[gibbs]
"""

from chemflow2 import GibbsReactor, Problem, Stream, StreamCondition

species = ["CH4", "H2O", "CO", "CO2", "H2"]

# 入口: CH4 1 mol + H2O 2 mol（S/C = 2）、850°C, 常圧
Feed = Stream(
    species, name="1. Feed", order=1,
    flows={"CH4": 1.0, "H2O": 2.0, "CO": 0, "CO2": 0, "H2": 0},
    condition=StreamCondition(T=850, P="0.1MPa", phase="gas"),
)
Out = Stream(species, name="2. Equilibrium", order=2,
             condition=StreamCondition(T=850, P="0.1MPa", phase="gas"))

# T・P は StreamCondition から取る（引数で上書きも可能）
G1 = GibbsReactor(inlet=Feed, outlet=Out, species=species, name="G1")

problem = Problem(streams=[Feed, Out], units=[G1], name="Steam Methane Reforming")

print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
sol = problem.solve()
print(sol)
print()
sol.print_report()

ch4_conv = 1 - Out.flow_of("CH4") / Feed.flow_of("CH4")
print(f"\nCH4 転化率 = {ch4_conv:.1%}")
print("平衡モル分率:")
for sp in species:
    print(f"  {sp:5s}: {Out.flow_of(sp) / float(Out.total_flow):.3f}")
