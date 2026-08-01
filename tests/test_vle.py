"""気液平衡の物性式（Antoine / Henry）と、それを使う制約のテスト。"""

import numpy as np
import pytest

from chemflow2 import Mixer, Problem, Separator, Stream, StreamCondition
from chemflow2.core.errors import ConstraintError
from chemflow2.core.vle import antoine_water_psat, henry_pa


# --- 物性式 ------------------------------------------------------------- #
def test_antoine_water_psat_at_boiling_point():
    """100 °C の飽和蒸気圧は 1 atm。"""
    assert antoine_water_psat(100.0) == pytest.approx(101325.0, rel=2e-3)


def test_antoine_water_psat_at_25c():
    """25 °C は約 3.17 kPa（NIST）。"""
    assert antoine_water_psat(25.0) == pytest.approx(3170.0, rel=5e-3)


def test_antoine_monotonic():
    ps = [antoine_water_psat(t) for t in (0, 25, 50, 75, 100)]
    assert all(a < b for a, b in zip(ps, ps[1:]))


def test_henry_unknown_species_is_none():
    assert henry_pa("CH3OCH3", 25.0) is None


def test_henry_co2_more_soluble_than_h2():
    """H = p/x なので、よく溶ける CO2 の方が H は小さい。"""
    assert henry_pa("CO2", 25.0) < henry_pa("H2", 25.0)


def test_henry_decreases_with_temperature():
    """昇温すると溶解度が下がる = H が大きくなる。"""
    assert henry_pa("CO2", 60.0) > henry_pa("CO2", 25.0)


# --- 制約 --------------------------------------------------------------- #
def _flash(dissolved, T=25.0, P="1.04MPaG", henry=None, sp=None, flows=None):
    """ReactOut を固定した 1 段フラッシュを組む。"""
    sp = sp or ["CO2", "CH4", "H2O", "CO", "H2"]
    flows = flows or {"CO2": 1.677814, "CH4": 0.243203, "H2O": 2.538494,
                      "CO": 3.549996, "H2": 4.912179}
    feed = Stream(sp, name="in", flows=flows)
    gas = Stream(sp, name="gas", condition=StreamCondition(T=T, P=P, phase="gas"))
    liq = Stream(sp, name="liq", condition=StreamCondition(T=T, P=P, phase="liquid"))
    p = Problem(streams=[feed, gas, liq], units=[Separator(feed, [gas, liq])])
    p.constrain_saturation(gas, "H2O", T=T, P=P)
    p.constrain_henry(gas, liq, dissolved, T=T, P=P, henry=henry)
    return p, feed, gas, liq


def _solve(p):
    """非負制約は付けない。

    「溶けない成分の液相 = 0」のように解が境界に乗る場合、bounds 付きの
    least_squares は境界へ漸近するだけで gtol 停止し、残差が 1e-8 台で残る。
    この手の系は全流量が正なので、素直に root で解いた方が桁違いに正確。
    """
    return p.solve()


def test_saturation_and_henry_close_the_dof():
    p, *_ = _flash(["H2", "CO", "CO2", "CH4"])
    assert p.degrees_of_freedom() == (10, 10)


def test_flash_matches_reference():
    """旧 ChemFlow pattern1 の凝縮器と一致する（CSV は小数6桁）。"""
    p, feed, gas, liq = _flash(["H2", "CO", "CO2", "CH4"])
    assert _solve(p).success
    expected_gas = {"H2": 4.911989, "CO": 3.549828, "CO2": 1.675074,
                    "CH4": 0.243186, "H2O": 0.028800}
    for f, v in expected_gas.items():
        assert gas.flow_of(f) == pytest.approx(v, abs=1e-6)
    # 成分収支
    for f in feed.formulas:
        assert gas.flow_of(f) + liq.flow_of(f) == pytest.approx(feed.flow_of(f), abs=1e-9)


def test_saturation_gives_expected_mole_fraction():
    """ガス側の水は y_sat = Psat/P になる。"""
    p, _, gas, _ = _flash(["H2", "CO", "CO2", "CH4"])
    assert _solve(p).success
    y = gas.flow_of("H2O") / float(gas.total_flow.eval())
    assert y == pytest.approx(antoine_water_psat(25.0) / (1.04e6 + 101325.0), rel=1e-6)


def test_henry_override_makes_insoluble():
    """henry= で H を大きく上書きすると、その成分は溶けなくなる。"""
    dissolved = ["H2", "CO", "CO2", "CH4"]
    base, _, _, base_liq = _flash(dissolved)
    assert _solve(base).success
    assert base_liq.flow_of("CO2") > 1e-4          # 既定では溶ける

    over, _, _, liq = _flash(dissolved, henry={"CO2": 1e15})
    assert _solve(over).success
    assert liq.flow_of("CO2") == pytest.approx(0.0, abs=1e-9)
    assert liq.flow_of("H2") > 0                   # 上書きしていない成分は溶けたまま


def test_henry_without_data_does_not_dissolve():
    """内蔵テーブルに無く上書きも無い成分は、液相 0 になる。"""
    sp = ["CO2", "H2O", "CH3OCH3"]
    p, _, _, liq = _flash(["CO2", "CH3OCH3"], sp=sp,
                          flows={"CO2": 1.5, "H2O": 2.5, "CH3OCH3": 0.8})
    assert p.degrees_of_freedom() == (6, 6)
    assert _solve(p).success
    assert henry_pa("CH3OCH3", 25.0) is None
    assert liq.flow_of("CH3OCH3") == pytest.approx(0.0, abs=1e-12)
    assert liq.flow_of("CO2") > 0


def test_saturation_requires_psat_for_non_water():
    sp = ["H2O", "CH3OH"]
    gas = Stream(sp, name="g")
    p = Problem(streams=[gas], units=[])
    with pytest.raises(ConstraintError):
        p.constrain_saturation(gas, "CH3OH", T=25, P="1atm")


def test_henry_rejects_solvent_itself():
    sp = ["H2O", "CO2"]
    gas, liq = Stream(sp, name="g"), Stream(sp, name="l")
    p = Problem(streams=[gas, liq], units=[])
    with pytest.raises(ConstraintError):
        p.constrain_henry(gas, liq, ["H2O"], T=25, P="1atm")


# --- pattern1 全体 ------------------------------------------------------ #
def test_pattern1_reproduces_chemflow():
    """改質器 + 水凝縮の一気通貫で、旧 ChemFlow の pattern1 と一致する。"""
    pytest.importorskip("cantera")
    import examples.example_pattern1 as p1

    assert p1.problem.degrees_of_freedom() == (19, 19)
    assert p1.problem.solve(bounds=(0, np.inf), ftol=1e-14, xtol=1e-14, gtol=1e-14).success

    for name, ref in p1.REFERENCE.items():
        stream = next(s for s in p1.problem.streams if s.name == name)
        for f, v in ref.items():
            assert stream.flow_of(f) == pytest.approx(v, abs=1e-5), f"{name}/{f}"

    # 逆算された Feed
    assert p1.RG_feed.flow_of("H2") == pytest.approx(0.885957, abs=1e-5)
    assert p1.RG_feed.flow_of("CH4") == pytest.approx(2.137257, abs=1e-5)
    assert p1.CO2_feed.flow_of("CO2") == pytest.approx(3.333756, abs=1e-5)
    assert p1.H2O_feed.flow_of("H2O") == pytest.approx(2.776607, abs=1e-5)

    # DryGas = plant3 / plant4 の Steam1
    steam1 = {"H2": 4.911989, "CO": 3.549828, "CO2": 1.675074,
              "CH4": 0.243186, "H2O": 0.028800}
    for f, v in steam1.items():
        assert p1.DryGas.flow_of(f) == pytest.approx(v, abs=1e-5)
