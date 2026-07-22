"""chemflow2: シンプルな定常状態プロセスシミュレーター。

宣言的に書ける:

    from chemflow2 import Stream, StreamCondition, Reaction, Reactor, Mixer, Separator, Problem

    components = ["H2", "CO", "CO2"]
    S1 = Stream(order=1, components=components, name="1. Feed",
                condition=StreamCondition(T=300, P=1, phase="gas"))
    ...
    problem = Problem(streams=[...], units=[...])
    problem.constrain(S1.total_flow, 165, name="S1 total flow")
    sol = problem.solve()
    sol.print_report()

設計:
    core/   ... 不変の核（Stream / Reaction / Unit契約 / Problem・求解）
    units/  ... 差し替え・追加可能な装置（Mixer / Reactor / Separator ...）
    io/     ... 出力（表 / CSV ...）
"""

from chemflow2.core import (
    BasisError,
    ChemflowError,
    Component,
    ComponentError,
    Constraint,
    ConstraintError,
    Expr,
    Problem,
    Reaction,
    ReactionError,
    Solution,
    SolveError,
    get_composition,
    Stream,
    StreamCondition,
    Unit,
    component_union,
    flows_on,
    get_component,
)
from chemflow2.core.pressure import parse_pressure
from chemflow2.io import export_mermaid, generate_mermaid, stream_table, to_csv, to_excel
from chemflow2.units import CanteraError, GibbsReactor, Mixer, Reactor, Separator, Splitter

__all__ = [
    # core
    "Stream",
    "StreamCondition",
    "Reaction",
    "Problem",
    "Solution",
    "Constraint",
    "Unit",
    "Expr",
    "Component",
    "get_component",
    "get_composition",
    "component_union",
    "flows_on",
    # units
    "Mixer",
    "Reactor",
    "Separator",
    "Splitter",
    "GibbsReactor",
    "parse_pressure",
    # io
    "stream_table",
    "to_csv",
    "to_excel",
    "generate_mermaid",
    "export_mermaid",
    # errors
    "ChemflowError",
    "ComponentError",
    "BasisError",
    "ConstraintError",
    "ReactionError",
    "SolveError",
    "CanteraError",
]
