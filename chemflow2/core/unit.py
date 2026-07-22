"""Unit: すべての装置の基底クラス（拡張ポイント）。

新しい装置を追加する手順はこれだけ:

1. `chemflow2/units/` に新しいファイルを作る
2. `Unit` を継承し、`residuals()` を実装する
   （= その装置が満たすべき「残差 = 0」のベクトルを返す）
3. `inlets` / `outlets` を返せるようにする（フロー図・検証用）

Problem は Unit の中身を知らず、`unit.residuals()` を集めるだけ。
これが「あとから Unit を足せる」設計の核。
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from chemflow2.core.stream import Stream


class Unit(ABC):
    """装置の基底クラス。

    Attributes
    ----------
    name : str | None
        装置名。
    """

    name: str | None = None

    @property
    @abstractmethod
    def inlets(self) -> list[Stream]:
        """入口ストリーム一覧。"""

    @property
    @abstractmethod
    def outlets(self) -> list[Stream]:
        """出口ストリーム一覧。"""

    @abstractmethod
    def residuals(self) -> np.ndarray:
        """満たすべき残差ベクトル（= 0 になるべき値）を返す。"""

    @property
    def streams(self) -> list[Stream]:
        return list(self.inlets) + list(self.outlets)

    def __repr__(self) -> str:
        ins = ",".join(s.name or "?" for s in self.inlets)
        outs = ",".join(s.name or "?" for s in self.outlets)
        return f"{type(self).__name__}({self.name!r}: [{ins}] -> [{outs}])"


# ---------------------------------------------------------------------- #
# 物質収支のための共通ヘルパー
# ---------------------------------------------------------------------- #
def component_union(streams: list[Stream], extra: list[str] | None = None) -> list[str]:
    """複数ストリーム（+ 追加成分）に現れる成分を、出現順で重複なく並べる。"""
    seen: dict[str, None] = {}
    for s in streams:
        for f in s.formulas:
            seen.setdefault(f, None)
    for f in extra or []:
        seen.setdefault(f, None)
    return list(seen)


def flows_on(stream: Stream, formulas: list[str]) -> np.ndarray:
    """ストリームのモル流量を、指定した成分順のベクトルに整列する（無い成分は 0）。"""
    return np.array([stream.flow_of(f) for f in formulas])
