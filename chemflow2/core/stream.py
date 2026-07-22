"""Stream: プロセスストリーム。

内部状態は各成分の **モル流量ベクトル** `molar_flows` のみ。
質量・体積・分率はすべてここから導出する（既存 chemflow の方針を踏襲）。

- `flows` を与える → 固定ストリーム（solve 対象外）。
- `flows` を与えない → 未知ストリーム（solve 対象）。制約やユニットの残差で決まる。

制約 DSL 用のプロパティ（total_flow, mole_fractions ...）は Expr を返す。
ユニットが物質収支を組むときは生の `molar_flows` / `flow_of()` を使う。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemflow2.core.component import Component, get_component
from chemflow2.core.errors import BasisError
from chemflow2.core.expr import Expr

_ABS_BASES = {"mol", "mole", "mass", "normal_volume"}
_FRAC_BASES = {"mole_frac", "mass_frac", "volume_frac"}


@dataclass
class StreamCondition:
    """ストリームの状態量（v1 ではメタ情報。物質収支には未使用）。

    温度・圧力・相を保持し、表示や将来のエネルギー収支・相平衡拡張で使う。
    """

    T: float | None = None       # 温度 [°C]
    P: float | None = None       # 圧力（数値 or "3MPaG" 等の文字列）
    phase: str | None = None     # "gas" / "liquid" / "solid" / "mixed"

    def __repr__(self) -> str:
        parts = []
        if self.T is not None:
            parts.append(f"T={self.T}")
        if self.P is not None:
            parts.append(f"P={self.P}")
        if self.phase is not None:
            parts.append(f"phase={self.phase!r}")
        return "StreamCondition(" + ", ".join(parts) + ")"


class Stream:
    """プロセスストリーム。

    Examples
    --------
    未知（solve 対象）::

        S1 = Stream(components=["H2", "CO", "CO2"], name="1. Feed", order=1,
                    condition=StreamCondition(T=300, P=1, phase="gas"))

    固定（basis 変換対応）::

        S0 = Stream(components=["H2", "N2"], flows={"H2": 20, "N2": 60})
        S0 = Stream(components=["H2", "N2"], flows={"H2": 0.75, "N2": 0.25},
                    basis="mole_frac", total=100)
    """

    def __init__(
        self,
        components: list[str],
        *,
        name: str | None = None,
        order: int | None = None,
        condition: StreamCondition | None = None,
        flows: dict[str, float] | None = None,
        basis: str = "mol",
        total: float | None = None,
        guess: dict[str, float] | list[float] | np.ndarray | None = None,
    ):
        self.name = name
        self.order = order
        self.condition = condition or StreamCondition()

        self.components: list[Component] = [get_component(f) for f in components]
        self.formulas: list[str] = [c.formula for c in self.components]
        self.index: dict[str, int] = {f: i for i, f in enumerate(self.formulas)}
        self.n = len(self.components)

        if flows is None:
            # 未知ストリーム：初期推定（既定は 1.0、guess で上書き可能）
            self.molar_flows = self._make_guess(guess)
            self.fixed = False
        else:
            self.molar_flows = self._to_molar(flows, basis, total)
            self.fixed = True

    def _make_guess(self, guess) -> np.ndarray:
        """未知ストリームの初期推定値ベクトルを作る。"""
        if guess is None:
            return np.ones(self.n)
        if isinstance(guess, dict):
            g = np.ones(self.n)
            for f, v in guess.items():
                if f not in self.index:
                    raise BasisError(f"guess のキー {f!r} が components に含まれていません")
                g[self.index[f]] = float(v)
            return g
        g = np.asarray(guess, dtype=float)
        if g.shape != (self.n,):
            raise BasisError(f"guess の長さ {g.shape} が成分数 {self.n} と一致しません")
        return g

    # ------------------------------------------------------------------ #
    # 初期化補助
    # ------------------------------------------------------------------ #
    def _to_molar(self, flows: dict[str, float], basis: str, total: float | None) -> np.ndarray:
        vals = np.zeros(self.n)
        for f, v in flows.items():
            if f not in self.index:
                raise BasisError(f"flows のキー {f!r} が components に含まれていません")
            vals[self.index[f]] = float(v)

        mws = np.array([c.mw for c in self.components])
        nvols = np.array([c.normal_volume for c in self.components])

        if basis in ("mol", "mole"):
            return vals
        if basis == "mass":
            return vals / mws
        if basis == "normal_volume":
            return vals / nvols
        if basis in _FRAC_BASES:
            if total is None:
                raise BasisError(f"basis={basis!r} には total が必要です")
            frac = vals / vals.sum()
            if basis == "mole_frac":
                return frac * total
            if basis == "mass_frac":
                return (frac * total) / mws
            return (frac * total) / nvols  # volume_frac
        raise BasisError(f"未知の basis: {basis!r}")

    # ------------------------------------------------------------------ #
    # ユニットが使う生アクセサ
    # ------------------------------------------------------------------ #
    def flow_of(self, formula: str) -> float:
        """指定成分のモル流量（無ければ 0）。数値で欲しいとき用。"""
        i = self.index.get(formula)
        return float(self.molar_flows[i]) if i is not None else 0.0

    # ------------------------------------------------------------------ #
    # 制約 DSL 用（Expr を返す）
    # ------------------------------------------------------------------ #
    def flow_expr(self, formula: str) -> Expr:
        """指定成分のモル流量を Expr で返す（制約に直截に使える）。

            problem.constrain(S5.flow_expr("CO"), 0)   # CO を全量ゼロに
        """
        i = self.index.get(formula)
        if i is None:
            return Expr(lambda: 0.0, f"{self.name}.flow[{formula}]")
        return Expr(lambda: float(self.molar_flows[i]), f"{self.name}.flow[{formula}]")

    # ------------------------------------------------------------------ #
    # 制約 DSL 用プロパティ（Expr を返す）
    # ------------------------------------------------------------------ #
    @property
    def total_flow(self) -> Expr:
        """総モル流量 [mol/時]。"""
        return Expr(lambda: float(np.sum(self.molar_flows)), f"{self.name}.total_flow")

    @property
    def total_mass_flow(self) -> Expr:
        mws = np.array([c.mw for c in self.components])
        return Expr(lambda: float(np.sum(self.molar_flows * mws)), f"{self.name}.total_mass_flow")

    @property
    def total_normal_volume_flow(self) -> Expr:
        nv = np.array([c.normal_volume for c in self.components])
        return Expr(lambda: float(np.sum(self.molar_flows * nv)), f"{self.name}.total_nvol_flow")

    @property
    def mole_fractions(self) -> Expr:
        def fn():
            t = np.sum(self.molar_flows)
            return self.molar_flows / t if t else self.molar_flows * 0.0
        return Expr(fn, f"{self.name}.mole_fractions")

    def frac_of(self, formula: str) -> Expr:
        """指定成分のモル分率（制約に使う）。"""
        i = self.index[formula]

        def fn():
            t = np.sum(self.molar_flows)
            return float(self.molar_flows[i] / t) if t else 0.0

        return Expr(fn, f"{self.name}.frac[{formula}]")

    # ------------------------------------------------------------------ #
    def __repr__(self) -> str:
        tag = "fixed" if self.fixed else "var"
        return f"Stream({self.name!r}, {self.formulas}, {tag})"
