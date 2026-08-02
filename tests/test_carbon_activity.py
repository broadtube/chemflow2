"""炭素活量 a_C と、それを使う制約のテスト。"""

import numpy as np
import pytest

pytest.importorskip("cantera")

from chemflow2 import GibbsReactor, Problem, Stream, StreamCondition
from chemflow2.core.carbon_activity import (
    carbon_activity,
    carbon_activity_of,
    equilibrium_constants,
)
from chemflow2.core.errors import CanteraError, ConstraintError

SPECIES = ["CO2", "CH4", "H2O", "CO", "H2"]
P = "1.04MPaG"


# --- 平衡定数 --------------------------------------------------------- #
def test_boudouard_K_is_one_near_700c():
    """2CO ⇌ C + CO2 の K = 1 は約 700 °C（Boudouard 平衡の古典的な基準温度）。"""
    K_B, _, _ = equilibrium_constants(700.0)
    assert K_B == pytest.approx(1.0, rel=1e-3)


def test_boudouard_is_exothermic_and_cracking_endothermic():
    """発熱の Boudouard は高温で K 低下、吸熱のメタン分解は高温で K 上昇。"""
    KB_lo, KM_lo, _ = equilibrium_constants(500.0)
    KB_hi, KM_hi, _ = equilibrium_constants(1000.0)
    assert KB_hi < KB_lo
    assert KM_hi > KM_lo


def test_constants_are_cached_but_temperature_dependent():
    assert equilibrium_constants(850.0) == equilibrium_constants(850.0)
    assert equilibrium_constants(850.0) != equilibrium_constants(900.0)


# --- 3 経路の一致 ----------------------------------------------------- #
def _reformer(T, flows=None):
    hot = StreamCondition(T=T, P=P, phase="gas")
    inlet = Stream(SPECIES, name="Mixed", condition=hot,
                   flows=flows or {"H2": 0.885957, "CH4": 2.137257,
                                   "CO2": 3.333756, "H2O": 2.776607})
    out = Stream(SPECIES, name="ReactOut", condition=hot,
                 guess={"CO": 3.5, "H2": 4.9, "CO2": 1.7, "CH4": 0.25, "H2O": 1.0})
    p = Problem(streams=[inlet, out],
                units=[GibbsReactor(inlet=inlet, outlet=out, species=SPECIES, T=T, P=P)])
    assert p.solve().success
    return p, inlet, out


def test_three_routes_agree_at_gibbs_equilibrium():
    """気相が内部平衡していれば 3 経路の a_C は一致する（差は WGS / SMR / ドライ改質）。"""
    for T in (850.0, 900.0):
        _, _, out = _reformer(T)
        routes = carbon_activity_of(out, per_route=True)
        values = list(routes.values())
        assert len(routes) == 3
        assert max(values) - min(values) < 1e-9 * max(values)


def test_routes_differ_for_non_equilibrated_gas():
    """非平衡の気相では経路ごとにばらつく（原料には CO が無いので 2 経路が 0）。"""
    _, inlet, _ = _reformer(850.0)
    routes = carbon_activity(
        {f: inlet.flow_of(f) / float(inlet.molar_flows.sum()) * 11.264 for f in SPECIES},
        850.0, per_route=True)
    assert routes["Boudouard"] == 0.0
    assert routes["CO_reduction"] == 0.0
    assert routes["CH4_cracking"] > 1.0        # 原料はメタン分解に対して過飽和


def test_known_values_for_pattern1_outlet():
    """pattern1 の改質条件での既知値（炭素析出線図の検証で使った値）。"""
    _, _, out850 = _reformer(850.0)
    _, _, out900 = _reformer(900.0)
    assert carbon_activity_of(out850) == pytest.approx(0.3944, abs=1e-4)
    assert carbon_activity_of(out900) == pytest.approx(0.2345, abs=1e-4)


def test_higher_temperature_is_safer():
    """改質では高温ほど a_C が下がる（炭素析出しにくい）。"""
    acts = [carbon_activity_of(_reformer(T)[2]) for T in (800.0, 850.0, 900.0)]
    assert acts[0] > acts[1] > acts[2]


def test_less_steam_raises_activity():
    """水蒸気を減らすと a_C が上がる。"""
    base = {"H2": 0.885957, "CH4": 2.137257, "CO2": 3.333756, "H2O": 2.776607}
    lean = dict(base, H2O=0.5)
    assert carbon_activity_of(_reformer(850.0, lean)[2]) > \
           carbon_activity_of(_reformer(850.0, base)[2])


# --- エラー ------------------------------------------------------------ #
def test_requires_T_and_P():
    s = Stream(SPECIES, name="s", flows={"CO": 1.0, "CO2": 1.0})
    with pytest.raises(CanteraError):
        carbon_activity_of(s)               # condition に T・P が無い


def test_zero_flow_stream_is_rejected():
    s = Stream(SPECIES, name="s", condition=StreamCondition(T=850, P=P),
               flows={f: 0.0 for f in SPECIES})
    with pytest.raises(CanteraError):
        carbon_activity_of(s)


# --- 制約 -------------------------------------------------------------- #
def test_constrain_rejects_non_positive_target():
    _, _, out = _reformer(850.0)
    p = Problem(streams=[out], units=[])
    with pytest.raises(ConstraintError):
        p.constrain_carbon_activity(out, 0.0)


def test_constraint_solves_the_carbon_limit():
    """a_C = 1 を課し、水蒸気供給を未知にすると「炭素析出しない最小の水蒸気量」が出る。"""
    from chemflow2 import Mixer

    T = 900.0
    hot = StreamCondition(T=T, P=P, phase="gas")
    rg = Stream(["H2", "CH4"], name="RG", condition=hot,
                flows={"H2": 0.885957, "CH4": 2.137257})
    co2f = Stream(["CO2"], name="CO2f", condition=hot, guess={"CO2": 2.0})
    h2of = Stream(["H2O"], name="H2Of", condition=hot, guess={"H2O": 1.0})  # ← 未知
    mixed = Stream(["H2", "CO2", "CH4", "H2O"], name="Mixed", condition=hot)
    out = Stream(SPECIES, name="ReactOut", condition=hot,
                 guess={"CO": 3.5, "H2": 4.9, "CO2": 1.7, "CH4": 0.25, "H2O": 1.0})

    p = Problem(streams=[rg, co2f, h2of, mixed, out],
                units=[Mixer([rg, co2f, h2of], mixed),
                       GibbsReactor(inlet=mixed, outlet=out, species=SPECIES, T=T, P=P)])
    # CO2 供給を決めるための指定（H2/CO 比）と、水蒸気供給を決めるための炭素活量
    p.constrain(out.flow_expr("H2"), 1.3837259 * out.flow_expr("CO"))
    p.constrain_carbon_activity(out, 1.0, T=T, P=P)

    n_var, n_eq = p.degrees_of_freedom()
    assert n_var == n_eq
    assert p.solve().success
    assert carbon_activity_of(out, T=T, P=P) == pytest.approx(1.0, rel=1e-6)
    assert h2of.flow_of("H2O") > 0.0


def test_constraint_residual_is_logarithmic():
    """残差は log(a_C) − log(target)。a_C が桁で振れても overflow しない。"""
    _, _, out = _reformer(850.0)                      # a_C ≈ 0.394
    p = Problem(streams=[out], units=[])
    p.constrain_carbon_activity(out, 1.0)
    r = float(p.constraints[-1].residuals()[0])
    assert r == pytest.approx(np.log(carbon_activity_of(out)), rel=1e-12)
    assert r < 0                                       # a_C < 1 なので負
