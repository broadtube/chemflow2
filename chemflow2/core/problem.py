"""Problem: ストリーム・ユニット・制約を集約し、連立方程式として解く。

グローバル状態を持たず、ユーザーが明示的にストリームとユニットを渡す。
循環系は「同じ Stream オブジェクトを複数ユニットで共有する」だけで表現できる
（テアストリームを未知ストリームとして宣言し、初期推定 1.0 から solve が収束させる）。

方程式の数え方（自由度）:
    変数   = 未知ストリームの成分数の総和
    方程式 = 各ユニットの residuals() + 各 constrain() の残差
両者が一致していないと SolveError を送出する。
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.optimize import least_squares, root

from chemflow2.core.errors import SolveError
from chemflow2.core.expr import Expr, value_of
from chemflow2.core.stream import Stream
from chemflow2.core.unit import Unit


class Constraint:
    """名前付きの残差（= 0 になるべき値を返す）。"""

    def __init__(self, fn: Callable[[], np.ndarray], name: str | None = None):
        self._fn = fn
        self.name = name

    def residuals(self) -> np.ndarray:
        return np.atleast_1d(self._fn())


class Solution:
    """solve() の結果。"""

    def __init__(self, success: bool, message: str, problem: "Problem", nfev: int = 0):
        self.success = success
        self.message = message
        self.problem = problem
        self.nfev = nfev

    def __repr__(self) -> str:
        state = "OK" if self.success else "FAILED"
        return f"Solution({state}, nfev={self.nfev}, {self.message!r})"

    def report(self) -> str:
        """ストリーム表（order 順）を文字列で返す。"""
        from chemflow2.io.table import stream_table

        return stream_table(self.problem.streams)

    def print_report(self) -> None:
        print(self.report())


class Problem:
    """フローシート全体を表す求解問題。"""

    def __init__(
        self,
        streams: list[Stream],
        units: list[Unit],
        name: str | None = None,
    ):
        self.name = name
        self.streams = list(streams)
        self.units = list(units)
        self.constraints: list[Constraint] = []

    # ------------------------------------------------------------------ #
    # 制約 API
    # ------------------------------------------------------------------ #
    def constrain(self, lhs, rhs=0.0, *, name: str | None = None) -> None:
        """`lhs == rhs` を課す（残差 = lhs - rhs）。

        lhs / rhs は Expr（例: ``S1.total_flow``）でも数値・配列でもよい。

            problem.constrain(S1.total_flow, 165, name="S1 total flow")
            problem.constrain(S1.total_mass_flow, S2.total_mass_flow)
        """
        self.constraints.append(
            Constraint(lambda: np.atleast_1d(value_of(lhs) - value_of(rhs)), name)
        )

    def constrain_recovery(
        self,
        inlet: Stream,
        outlet: Stream,
        fracs: dict[str, float],
        *,
        name: str | None = None,
    ) -> None:
        """成分回収率を指定する: ``outlet.flow(c) = frac · inlet.flow(c)``。

        Separator（収支のみ課すノード）と組み合わせて分離を簡潔に閉じる。

            # H2O は全量 Liquid へ、それ以外は全量 Gas へ（= Liquid への回収率 0）
            problem.constrain_recovery(S3, S5, {"H2O": 1.0, "H2": 0.0, "CO": 0.0,
                                                "CO2": 0.0, "CH3OCH3": 0.0})
        """
        for formula, frac in fracs.items():
            label = name or f"recovery[{formula}] {inlet.name}->{outlet.name}={frac}"
            self.constraints.append(
                Constraint(
                    lambda f=formula, r=frac: np.atleast_1d(outlet.flow_of(f) - r * inlet.flow_of(f)),
                    label,
                )
            )

    def constrain_fracs(self, stream: Stream, fracs: dict[str, float], *, name: str | None = None) -> None:
        """ストリームのモル分率を指定する。

            problem.constrain_fracs(S1, {"H2": 0.48, "CO": 0.24, "CO2": 0.28})

        注意: 分率は和が 1 になるため、全成分を指定すると過剰拘束になりやすい。
        通常は 1 成分を残し、その分を total などで閉じる。
        """
        for formula, target in fracs.items():
            expr = stream.frac_of(formula)
            label = name or f"{stream.name}.frac[{formula}]={target}"
            self.constraints.append(
                Constraint(lambda e=expr, t=target: np.atleast_1d(e.eval() - t), label)
            )

    # ------------------------------------------------------------------ #
    # 求解
    # ------------------------------------------------------------------ #
    def _variables(self) -> list[Stream]:
        return [s for s in self.streams if not s.fixed]

    def _pack(self) -> np.ndarray:
        arrs = [s.molar_flows for s in self._variables()]
        return np.concatenate(arrs) if arrs else np.array([])

    def _unpack(self, x: np.ndarray) -> None:
        i = 0
        for s in self._variables():
            s.molar_flows = x[i : i + s.n].copy()
            i += s.n

    def _residuals(self, x: np.ndarray) -> np.ndarray:
        self._unpack(x)
        parts: list[np.ndarray] = []
        for u in self.units:
            r = u.residuals()
            if r is not None and len(r):
                parts.append(np.atleast_1d(r))
        for c in self.constraints:
            r = c.residuals()
            if r is not None and len(r):
                parts.append(np.atleast_1d(r))
        return np.concatenate(parts) if parts else np.array([])

    def degrees_of_freedom(self) -> tuple[int, int]:
        """(変数の数, 方程式の数) を返す。"""
        x0 = self._pack()
        return len(x0), len(self._residuals(x0))

    def solve(self, *, bounds: tuple | None = None, tol: float = 1e-8, **kwargs) -> Solution:
        """連立方程式を解く。

        Parameters
        ----------
        bounds : tuple | None
            (下限, 上限)。例 ``(0, np.inf)`` で非負制約。指定時は least_squares を使う。
        tol : float
            収束判定に使う最終残差ノルムの閾値。solver の早期停止フラグではなく、
            実際に ‖residual‖ が十分小さいかを直接確認する。
        **kwargs
            scipy.optimize.root / least_squares に渡す追加引数。
        """
        x0 = self._pack()
        if len(x0) == 0:
            return Solution(True, "変数がありません", self)

        n_var, n_eq = self.degrees_of_freedom()
        if n_var != n_eq:
            kind = "過剰決定" if n_eq > n_var else "自由度不足"
            raise SolveError(f"{kind}: 変数 {n_var} 個 / 方程式 {n_eq} 個")

        if bounds is not None:
            res = least_squares(self._residuals, x0, bounds=bounds, **kwargs)
            self._unpack(res.x)
            resid_norm = float(np.linalg.norm(res.fun))
            ok = resid_norm < tol  # status ではなく最終残差そのもので判定
            msg = f"‖residual‖={resid_norm:.2e} ({res.message})"
            return Solution(ok, msg, self, res.nfev)

        # root を複数メソッドで試行（method 明示時はフォールバックしない）
        methods = [kwargs.pop("method")] if "method" in kwargs else ["hybr", "lm", "df-sane"]
        last = None
        for m in methods:
            try:
                res = root(self._residuals, x0, method=m, **kwargs)
            except Exception:
                continue
            last = res
            if res.success:
                self._unpack(res.x)
                return Solution(True, res.message, self, res.get("nfev", 0))

        self._unpack(last.x if last is not None else x0)
        raise SolveError(f"収束しませんでした: {last.message if last else 'すべての手法が失敗'}")

    def __repr__(self) -> str:
        return f"Problem({self.name!r}, streams={len(self.streams)}, units={len(self.units)})"
