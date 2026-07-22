"""Separator: 1 入口を複数出口に分ける一般分離ノード。

このユニット自体は「全体の物質収支（成分ごとに 入口 = Σ出口）」だけを課す。
どの成分がどの出口へどれだけ行くか（分配・分離率）は、あえてここで決めず、
Problem 側の constrain() / constrain_fracs() で指定する。

これにより、単純な分流器（Splitter）から、成分ごとに分離率の異なる分離塔まで、
同じ 1 つのユニットで表現でき、拡張が制約の追加だけで済む。
"""

from __future__ import annotations

import numpy as np

from chemflow2.core.stream import Stream
from chemflow2.core.unit import Unit, component_union, flows_on


class Separator(Unit):
    """分離器 / 分流器。

        Sep1 = Separator(inlet=S1, outlet=[S4, S5], name="Sep1")

    残差: 成分ごとに ``入口 - Σ出口 = 0``。
    出口間の分配は制約で閉じる（未指定だと自由度不足になる）。
    """

    def __init__(
        self,
        inlet: Stream,
        outlet: list[Stream],
        *,
        components: list[str] | None = None,
        name: str | None = None,
    ):
        self._inlet = inlet
        self._outlets = list(outlet)
        self.components = components
        self.name = name

    @property
    def inlets(self) -> list[Stream]:
        return [self._inlet]

    @property
    def outlets(self) -> list[Stream]:
        return self._outlets

    def residuals(self) -> np.ndarray:
        formulas = component_union([self._inlet] + self._outlets, extra=self.components)
        outflow = np.zeros(len(formulas))
        for s in self._outlets:
            outflow += flows_on(s, formulas)
        return flows_on(self._inlet, formulas) - outflow
