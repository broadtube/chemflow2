"""基本動作の最小テスト。

    PYTHONPATH=. python3 -m pytest tests/ -q
"""

import numpy as np
import pytest

from chemflow2 import (
    Expr,
    Mixer,
    Problem,
    Reaction,
    ReactionError,
    Reactor,
    Separator,
    Stream,
    SolveError,
)


def test_component_mw():
    s = Stream(["H2O"])
    assert abs(s.components[0].mw - 18.015) < 0.01


def test_mixer_balance():
    a = Stream(["N2", "H2"], flows={"N2": 1, "H2": 3})
    b = Stream(["N2", "H2"], flows={"N2": 2, "H2": 1})
    out = Stream(["N2", "H2"], name="out")
    p = Problem([a, b, out], [Mixer([a, b], out)])
    sol = p.solve()
    assert sol.success
    assert out.flow_of("N2") == pytest.approx(3.0)
    assert out.flow_of("H2") == pytest.approx(4.0)


def test_reactor_conversion():
    feed = Stream(["N2O4", "NO2"], flows={"N2O4": 10, "NO2": 0})
    out = Stream(["N2O4", "NO2"], name="out")
    rxn = Reaction({"N2O4": -1, "NO2": 2})
    p = Problem([feed, out], [Reactor(feed, out, [rxn], key_component="N2O4", conversion=0.5)])
    sol = p.solve()
    assert sol.success
    assert out.flow_of("N2O4") == pytest.approx(5.0)
    assert out.flow_of("NO2") == pytest.approx(10.0)


def test_recycle_converges_and_balances():
    comps = ["N2O4", "NO2"]
    feed = Stream(comps, flows={"N2O4": 10, "NO2": 0})
    rec = Stream(comps, name="rec")
    mixed = Stream(comps, name="mixed")
    rout = Stream(comps, name="rout")
    prod = Stream(comps, name="prod")
    rxn = Reaction({"N2O4": -1, "NO2": 2})
    p = Problem(
        [feed, rec, mixed, rout, prod],
        [
            Mixer([feed, rec], mixed),
            Reactor(mixed, rout, [rxn], key_component="N2O4", conversion=0.5),
            Separator(rout, [prod, rec]),
        ],
    )
    p.constrain(Expr(lambda: rec.molar_flows - 0.3 * rout.molar_flows))
    sol = p.solve(bounds=(0, np.inf))
    assert sol.success
    # 元素収支: 入った N は出た N（Product）と等しい（N2O4=2N, NO2=1N）
    n_in = 2 * feed.flow_of("N2O4") + feed.flow_of("NO2")
    n_out = 2 * prod.flow_of("N2O4") + prod.flow_of("NO2")
    assert n_in == pytest.approx(n_out)


def test_reaction_atom_balance_ok():
    # 2 CO + 4 H2 -> CH3OCH3 + H2O は原子保存する
    Reaction({"CO": -2, "H2": -4, "CH3OCH3": 1, "H2O": 1})


def test_reaction_atom_balance_bad_raises():
    # H2O の係数をタイポ（原子が消える）
    with pytest.raises(ReactionError):
        Reaction({"CO": -2, "H2": -4, "CH3OCH3": 1})


def test_selectivity_sum_must_be_one():
    feed = Stream(["N2O4", "NO2"], flows={"N2O4": 10, "NO2": 0})
    out = Stream(["N2O4", "NO2"], name="out")
    rxn = Reaction({"N2O4": -1, "NO2": 2})
    with pytest.raises(ValueError):
        Reactor(feed, out, [rxn, rxn], key_component="N2O4",
                conversion=0.5, selectivities=[0.7, 0.7])


def test_flow_expr_constraint():
    feed = Stream(["N2O4", "NO2"], flows={"N2O4": 10, "NO2": 5})
    a = Stream(["N2O4", "NO2"], name="a")
    b = Stream(["N2O4", "NO2"], name="b")
    p = Problem([feed, a, b], [Separator(feed, [a, b])])
    p.constrain(b.flow_expr("N2O4"), 0.0)
    p.constrain(b.flow_expr("NO2"), 2.0)
    sol = p.solve()
    assert sol.success
    assert a.flow_of("N2O4") == pytest.approx(10.0)
    assert a.flow_of("NO2") == pytest.approx(3.0)


def test_dof_mismatch_raises():
    feed = Stream(["N2O4", "NO2"], flows={"N2O4": 10, "NO2": 0})
    out = Stream(["N2O4", "NO2"], name="out")
    # 未知 out に対しユニットも制約も無い → 自由度不足
    p = Problem([feed, out], [])
    with pytest.raises(SolveError):
        p.solve()
