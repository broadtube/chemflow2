"""KineticReactor（reaction_rate バックエンド）のテスト。

reaction_rate が未導入なら skip（`pip install -e ../reaction_rate` で有効化）。
"""

import pytest

pytest.importorskip("reaction_rate")

from chemflow2 import KineticReactor, Problem, Stream, StreamCondition

SPECIES = ["CO", "CO2", "H2", "H2O", "CH3OH", "CH3OCH3", "CH3COOCH3"]
M_SYN, M_DEH, M_MA = 6.842e-4, 3.789e-4, 1.674e-3
N_MOLH = 1.859e-4 * 3600.0
Y_CO = (1 - 0.03) / (1 + 1.5)


def _build():
    cond = StreamCondition(T=250, P="5MPaG", phase="gas")
    Feed = Stream(SPECIES, name="F", condition=cond,
                  flows={"CO": Y_CO * N_MOLH, "CO2": 0.03 * N_MOLH, "H2": 1.5 * Y_CO * N_MOLH,
                         "H2O": 0, "CH3OH": 0, "CH3OCH3": 0, "CH3COOCH3": 0})
    Mid = Stream(SPECIES, name="M", condition=cond)
    Out = Stream(SPECIES, name="O", condition=cond)
    R1 = KineticReactor(inlet=Feed, outlet=Mid,
                        masses={"synthesis": M_SYN, "dehydration": M_DEH},
                        models={"synthesis": "KOGAS", "dehydration": "ZSM5"}, k_eq3="thermo")
    R2 = KineticReactor(inlet=Mid, outlet=Out,
                        masses={"carbonylation": M_MA},
                        models={"carbonylation": "DTU-Cheung2007-2"})
    problem = Problem(streams=[Feed, Mid, Out], units=[R1, R2])
    return problem, Feed, Mid, Out


def _direct():
    from reaction_rate.reactors import CatalystBed, pfr
    f = {"CO": Y_CO * N_MOLH / 3600, "CO2": 0.03 * N_MOLH / 3600, "H2": 1.5 * Y_CO * N_MOLH / 3600,
         "H2O": 0, "CH3OH": 0, "DME": 0, "MA": 0}
    r1 = pfr(f, 523.15, 51.01325, CatalystBed({"synthesis": M_SYN, "dehydration": M_DEH}),
             models={"synthesis": "KOGAS", "dehydration": "ZSM5"}, k_eq3="thermo")
    r2 = pfr(r1.outlet(), 523.15, 51.01325, CatalystBed({"carbonylation": M_MA}),
             models={"carbonylation": "DTU-Cheung2007-2"})
    return r2.outlet()


def test_dof_balanced_and_converges():
    problem, *_ = _build()
    assert problem.degrees_of_freedom() == (14, 14)
    sol = problem.solve()
    assert sol.success


def test_matches_direct_reaction_rate():
    problem, Feed, Mid, Out = _build()
    assert problem.solve().success
    o = _direct()
    alias = {"CH3OCH3": "DME", "CH3COOCH3": "MA"}
    for sp in SPECIES:
        direct_molh = o.get(alias.get(sp, sp), 0.0) * 3600.0
        assert Out.flow_of(sp) == pytest.approx(direct_molh, abs=1e-6)


def test_mass_conserved_and_MA_formed():
    problem, Feed, Mid, Out = _build()
    assert problem.solve().success
    # 質量保存（chemflow2 の molmass MW で総質量流を比較）
    m_in = float(Feed.total_mass_flow)
    m_out = float(Out.total_mass_flow)
    assert m_out == pytest.approx(m_in, rel=1e-5)
    # MA（酢酸メチル）が生成している
    assert Out.flow_of("CH3COOCH3") > 0.0
