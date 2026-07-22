"""Expr: 遅延評価スカラー/ベクトル。

制約 DSL の核。`S1.total_flow` のようなプロパティは即座に float を返すのではなく
Expr を返す。これにより:

    problem.constrain(S1.total_flow, 165)

と書いた時点では「S1 の総流量」という *参照* を捕まえておき、solve 中に
その時々の molar_flows から値を計算できる。

同時に __float__ / __format__ を実装しているので、solve 後の表示では

    print(f"{S1.total_flow:.2f}")   # 解いた値が出る
    float(S1.total_flow)

のように通常の数値として扱える。ビルド時は遅延・表示時は数値、同じ属性で両立する。
"""

from __future__ import annotations

from typing import Callable

import numpy as np


def value_of(x) -> float | np.ndarray:
    """Expr なら評価し、そうでなければそのまま返す。"""
    if isinstance(x, Expr):
        return x.eval()
    return x


class Expr:
    """molar_flows に対する遅延評価式。

    Parameters
    ----------
    fn : Callable[[], float | np.ndarray]
        評価時に現在の状態から値を計算する関数。
    label : str
        表示・デバッグ用のラベル。
    """

    def __init__(self, fn: Callable[[], float | np.ndarray], label: str = ""):
        self._fn = fn
        self.label = label

    def eval(self) -> float | np.ndarray:
        return self._fn()

    # --- 数値としての振る舞い（solve 後の表示用）---
    def __float__(self) -> float:
        return float(self.eval())

    def __format__(self, spec: str) -> str:
        return format(float(self.eval()), spec)

    def __repr__(self) -> str:
        try:
            v = self.eval()
        except Exception:
            return f"Expr({self.label!r})"
        return f"{v:.6g}" if np.isscalar(v) else f"Expr({self.label!r}={v})"

    # --- 制約を組み立てるための最小限の算術 ---
    def __sub__(self, other) -> "Expr":
        return Expr(lambda: value_of(self) - value_of(other), f"({self.label}-{_lab(other)})")

    def __rsub__(self, other) -> "Expr":
        return Expr(lambda: value_of(other) - value_of(self), f"({_lab(other)}-{self.label})")

    def __add__(self, other) -> "Expr":
        return Expr(lambda: value_of(self) + value_of(other), f"({self.label}+{_lab(other)})")

    __radd__ = __add__

    def __mul__(self, other) -> "Expr":
        return Expr(lambda: value_of(self) * value_of(other), f"({self.label}*{_lab(other)})")

    __rmul__ = __mul__

    def __truediv__(self, other) -> "Expr":
        return Expr(lambda: value_of(self) / value_of(other), f"({self.label}/{_lab(other)})")


def _lab(x) -> str:
    return x.label if isinstance(x, Expr) else repr(x)
