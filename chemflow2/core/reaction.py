"""Reaction: 反応の化学量論を表す独立オブジェクト。

Reactor から切り離しているので、1 つの反応を複数の Reactor で使い回したり、
複数反応を 1 つの Reactor に渡したりできる。

生成時に **原子収支（Σ_i ν_i · 元素ベクトル_i == 0）を検査** する。
係数のタイポで原子が生成・消滅する「閉じない物質収支」を黙って通さないための安全網。
"""

from __future__ import annotations

from dataclasses import dataclass, field

from chemflow2.core.component import get_composition
from chemflow2.core.errors import ReactionError


@dataclass
class Reaction:
    """化学反応。

    Parameters
    ----------
    stoich : dict[str, float]
        化学量論係数。反応物は負、生成物は正。
        例: {"CO": -2, "H2": -4, "CH3OCH3": 1, "H2O": 1}
    name : str | None
        表示用の名前。
    check : bool
        True（既定）なら原子収支を検査し、崩れていれば ReactionError。
    atol : float
        原子収支許容誤差（分数係数対応）。
    """

    stoich: dict[str, float]
    name: str | None = None
    check: bool = True
    atol: float = 1e-9
    species: list[str] = field(init=False)

    def __post_init__(self):
        self.species = list(self.stoich.keys())
        if self.check:
            self.validate_atom_balance()

    def coeff(self, formula: str) -> float:
        """指定成分の化学量論係数（無ければ 0）。"""
        return float(self.stoich.get(formula, 0.0))

    def element_balance(self) -> dict[str, float]:
        """元素ごとの収支 Σ_i ν_i · 原子数_i を返す（0 が理想）。"""
        balance: dict[str, float] = {}
        for formula, nu in self.stoich.items():
            for element, count in get_composition(formula).items():
                balance[element] = balance.get(element, 0.0) + nu * count
        return balance

    def validate_atom_balance(self) -> None:
        """原子収支が崩れていれば ReactionError を送出する。"""
        bad = {el: v for el, v in self.element_balance().items() if abs(v) > self.atol}
        if bad:
            detail = ", ".join(f"{el}: {v:+g}" for el, v in bad.items())
            raise ReactionError(
                f"反応 {self.name or self.stoich!r} の原子収支が崩れています ({detail})。"
                f" 係数を確認してください。"
            )

    def __repr__(self) -> str:
        return f"Reaction({self.name!r}, {self.stoich})"
