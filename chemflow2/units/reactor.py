"""Reactor: 転化率 + 選択率で複数反応を扱う反応器。

反応進行度モデル:
    key 成分の総消費量  = conversion * (入口の key 流量)
    反応 j への配分     = 総消費量 * selectivities[j]
    反応 j の進行度 ξ_j = (配分量) / |stoich_j[key]|
    出口成分 i          = 入口_i + Σ_j ξ_j * stoich_j[i]

selectivities は「key 消費量の各反応への配分割合」。単一反応なら [1.0]。
"""

from __future__ import annotations

import numpy as np

from chemflow2.core.reaction import Reaction
from chemflow2.core.stream import Stream
from chemflow2.core.unit import Unit, component_union, flows_on


class Reactor(Unit):
    """反応器。

        rxn1 = Reaction(stoich={"CO": -2, "H2": -4, "CH3OCH3": 1, "H2O": 1}, name="Tandem")
        R1 = Reactor(inlet=S2, outlet=S3, reactions=[rxn1],
                     key_component="CO", conversion=0.8, selectivities=[1], name="R1")

    残差: 成分ごとに ``出口 - (入口 + Σ 進行度*係数) = 0``。
    出口ストリームは全生成物成分を components に含めておくこと。
    """

    def __init__(
        self,
        inlet: Stream,
        outlet: Stream,
        reactions: list[Reaction],
        *,
        key_component: str,
        conversion: float,
        selectivities: list[float] | None = None,
        name: str | None = None,
    ):
        self._inlet = inlet
        self._outlet = outlet
        self.reactions = list(reactions)
        self.key = key_component
        self.conversion = float(conversion)
        if selectivities is None:
            selectivities = [1.0] if len(self.reactions) == 1 else None
        if selectivities is None or len(selectivities) != len(self.reactions):
            raise ValueError("selectivities は reactions と同じ長さで指定してください")
        self.selectivities = [float(s) for s in selectivities]
        # Σsel != 1 だと指定 conversion と実効 conversion がずれる（footgun）
        total_sel = sum(self.selectivities)
        if abs(total_sel - 1.0) > 1e-9:
            raise ValueError(
                f"selectivities の和は 1 である必要があります（現在 {total_sel:g}）。"
                f" key 消費量の各反応への配分割合として指定してください。"
            )
        self.name = name

    @property
    def inlets(self) -> list[Stream]:
        return [self._inlet]

    @property
    def outlets(self) -> list[Stream]:
        return [self._outlet]

    def residuals(self) -> np.ndarray:
        species = [s for rxn in self.reactions for s in rxn.species]
        formulas = component_union([self._inlet, self._outlet], extra=species)

        predicted = flows_on(self._inlet, formulas).astype(float)

        key_in = self._inlet.flow_of(self.key)
        total_key_consumed = self.conversion * key_in

        for rxn, sel in zip(self.reactions, self.selectivities):
            nu_key = rxn.coeff(self.key)  # 反応物なので負
            if nu_key == 0:
                raise ValueError(f"反応 {rxn.name!r} に key 成分 {self.key!r} が含まれていません")
            extent = (total_key_consumed * sel) / (-nu_key)
            for i, f in enumerate(formulas):
                predicted[i] += extent * rxn.coeff(f)

        return flows_on(self._outlet, formulas) - predicted
