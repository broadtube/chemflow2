"""chemflow2 の不変の核。

ここには「ストリーム・反応・装置の契約・求解」という、変わりにくい概念だけを置く。
個々の装置は units/、出力は io/ に分離する。
"""

from chemflow2.core.component import Component, get_component, get_composition
from chemflow2.core.errors import (
    BasisError,
    ChemflowError,
    ComponentError,
    ConstraintError,
    ReactionError,
    SolveError,
)
from chemflow2.core.expr import Expr
from chemflow2.core.problem import Constraint, Problem, Solution
from chemflow2.core.reaction import Reaction
from chemflow2.core.stream import Stream, StreamCondition
from chemflow2.core.unit import Unit, component_union, flows_on

__all__ = [
    "Component",
    "get_component",
    "get_composition",
    "Expr",
    "Stream",
    "StreamCondition",
    "Reaction",
    "Unit",
    "component_union",
    "flows_on",
    "Problem",
    "Solution",
    "Constraint",
    "ChemflowError",
    "ComponentError",
    "BasisError",
    "ConstraintError",
    "ReactionError",
    "SolveError",
]
