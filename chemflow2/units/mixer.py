"""Mixer: 複数入口を 1 つの出口に合流する。"""

from __future__ import annotations

import numpy as np

from chemflow2.core.stream import Stream
from chemflow2.core.unit import Unit, component_union, flows_on


class Mixer(Unit):
    """混合器。

        M1 = Mixer(inlet=[S1, S14], outlet=S2, name="M1")

    残差: 成分ごとに ``出口 - Σ入口 = 0``（物質収支）。
    """

    def __init__(self, inlet: list[Stream], outlet: Stream, *, name: str | None = None):
        self._inlets = list(inlet)
        self._outlet = outlet
        self.name = name

    @property
    def inlets(self) -> list[Stream]:
        return self._inlets

    @property
    def outlets(self) -> list[Stream]:
        return [self._outlet]

    def residuals(self) -> np.ndarray:
        formulas = component_union(self._inlets + [self._outlet])
        inflow = np.zeros(len(formulas))
        for s in self._inlets:
            inflow += flows_on(s, formulas)
        return flows_on(self._outlet, formulas) - inflow
