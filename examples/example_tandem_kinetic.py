"""サンプル: KineticReactor による単管タンデム DME→MA 合成。

reaction_rate の速度論 PFR を chemflow2 の Unit として使い、
    床① ハイブリッド（Cu/ZnO 合成 + ZSM-5 脱水）→ 床② H-MOR カルボニル化
を 2 つの KineticReactor 直列で解く。条件は reaction_rate の
examples/fig_tandem_sty_stack_zsm5.py と同一（250℃, 51 bar, SV 5000/h）。

要 reaction_rate: pip install -e ../reaction_rate
（本ファイルは PYTHONPATH=../reaction_rate/src でも動く）。
"""

from chemflow2 import KineticReactor, Problem, Stream, StreamCondition

# 反応網の全 7 種（DME=CH3OCH3, MA=CH3COOCH3 は実式で。KineticReactor が短縮キーに変換）
SPECIES = ["CO", "CO2", "H2", "H2O", "CH3OH", "CH3OCH3", "CH3COOCH3"]

# 触媒質量 [kg]（fig_tandem_sty_stack_zsm5 と同一: ρ×V, f_syn=1/1.9）
M_SYN, M_DEH, M_MA = 6.842e-4, 3.789e-4, 1.674e-3

# 供給モル流 [mol/h]（N=1.859e-4 mol/s × 3600 = 0.6693 mol/h; H2/CO=1.5, CO2 3%）
N_MOLH = 1.859e-4 * 3600.0
Y_CO = (1 - 0.03) / (1 + 1.5)

cond = StreamCondition(T=250, P="5MPaG", phase="gas")   # 5 MPa ゲージ + atm = 51.0 bar
Feed = Stream(SPECIES, name="1. Feed", order=1, condition=cond,
              flows={"CO": Y_CO * N_MOLH, "CO2": 0.03 * N_MOLH, "H2": 1.5 * Y_CO * N_MOLH,
                     "H2O": 0, "CH3OH": 0, "CH3OCH3": 0, "CH3COOCH3": 0})
Mid = Stream(SPECIES, name="2. hybrid out", order=2, condition=cond)   # 未知
Out = Stream(SPECIES, name="3. MA out", order=3, condition=cond)       # 未知

R1 = KineticReactor(inlet=Feed, outlet=Mid,
                    masses={"synthesis": M_SYN, "dehydration": M_DEH},
                    models={"synthesis": "KOGAS", "dehydration": "ZSM5"},
                    k_eq3="thermo", name="hybrid bed")
R2 = KineticReactor(inlet=Mid, outlet=Out,
                    masses={"carbonylation": M_MA},
                    models={"carbonylation": "DTU-Cheung2007-2"},
                    name="MA bed")

problem = Problem(streams=[Feed, Mid, Out], units=[R1, R2], name="Tandem DME->MA (kinetic)")


def main():
    print("自由度 (変数, 方程式):", problem.degrees_of_freedom())
    sol = problem.solve()
    print(sol)
    sol.print_report()

    print("\n出口モル流 [mol/h]:")
    for sp in SPECIES:
        print(f"  {sp:10s} {Out.flow_of(sp):.4e}")

    # 直接 reaction_rate と突き合わせ（同一結果になるはず）
    from reaction_rate.reactors import CatalystBed, pfr
    f = {"CO": Y_CO * N_MOLH / 3600, "CO2": 0.03 * N_MOLH / 3600, "H2": 1.5 * Y_CO * N_MOLH / 3600,
         "H2O": 0, "CH3OH": 0, "DME": 0, "MA": 0}
    r1 = pfr(f, 523.15, 51.01325, CatalystBed({"synthesis": M_SYN, "dehydration": M_DEH}),
             models={"synthesis": "KOGAS", "dehydration": "ZSM5"}, k_eq3="thermo")
    r2 = pfr(r1.outlet(), 523.15, 51.01325, CatalystBed({"carbonylation": M_MA}),
             models={"carbonylation": "DTU-Cheung2007-2"})
    o = r2.outlet()
    alias = {"CH3OCH3": "DME", "CH3COOCH3": "MA"}
    print("\n直接 reaction_rate との差 [mol/h]:")
    for sp in SPECIES:
        direct = o.get(alias.get(sp, sp), 0.0) * 3600
        print(f"  {sp:10s} chemflow2={Out.flow_of(sp):.4e}  direct={direct:.4e}")


if __name__ == "__main__":
    main()
