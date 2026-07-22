"""Splitter: 1 入口を、組成そのままに比率分割する（ティー / パージ分岐）。

Separator との違い:
    Separator ... 物質収支だけを課す一般分離ノード（組成が変わる。分配は制約で指定）
    Splitter  ... 各出口 = 入口 × 比率。全出口が入口と同じ組成になる（分流器・パージ）

循環・パージがこれ 1 つで直截に書ける（生の Expr 制約が不要になる）。
"""

from __future__ import annotations

import numpy as np

from chemflow2.core.stream import Stream
from chemflow2.core.unit import Unit, component_union, flows_on


class Splitter(Unit):
    """分流器。

        SP1 = Splitter(inlet=Rout, outlet=[Product, Recycle], ratios=[0.7, 0.3], name="SP1")

    残差: 出口ごと・成分ごとに ``出口_k - 比率_k × 入口 = 0``。
    比率の和は 1 でなければならない（検査する）。
    """

    def __init__(
        self,
        inlet: Stream,
        outlet: list[Stream],
        ratios: list[float],
        *,
        name: str | None = None,
    ):
        self._inlet = inlet
        self._outlets = list(outlet)
        self.ratios = [float(r) for r in ratios]
        if len(self._outlets) != len(self.ratios):
            raise ValueError("outlet と ratios は同じ長さで指定してください")
        total = sum(self.ratios)
        if abs(total - 1.0) > 1e-9:
            raise ValueError(f"ratios の和は 1 である必要があります（現在 {total:g}）")
        self.name = name

    @property
    def inlets(self) -> list[Stream]:
        return [self._inlet]

    @property
    def outlets(self) -> list[Stream]:
        return self._outlets

    def residuals(self) -> np.ndarray:
        formulas = component_union([self._inlet] + self._outlets)
        inflow = flows_on(self._inlet, formulas)
        parts = [flows_on(out, formulas) - r * inflow for out, r in zip(self._outlets, self.ratios)]
        return np.concatenate(parts)
