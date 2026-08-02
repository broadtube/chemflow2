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

from chemflow2.core.carbon_activity import carbon_activity_of
from chemflow2.core.errors import ConstraintError, SolveError
from chemflow2.core.expr import Expr, value_of
from chemflow2.core.pressure import parse_pressure
from chemflow2.core.stream import Stream
from chemflow2.core.unit import Unit
from chemflow2.core.vle import antoine_water_psat, henry_pa


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

    def constrain_saturation(
        self,
        gas: Stream,
        formula: str = "H2O",
        *,
        T: float,
        P: float | str,
        psat: float | None = None,
        name: str | None = None,
    ) -> None:
        """ガス側の成分を飽和量に固定する（**方程式 1 本**）。

            gas(formula) = y_sat/(1 − y_sat) · Σ gas(formula 以外)
            y_sat = P_sat(T) / P

        凝縮器を Separator（収支のみ）で組み、その分配を物理で閉じるための制約。
        ``constrain_recovery`` の物性版で、混ぜて使える（例: 水は飽和で決め、
        有機物は回収率で指定）。

            Cond = Separator(ReactOut, [DryGas, Condensate])
            problem.constrain_saturation(DryGas, "H2O", T=25, P="1.04MPaG")

        psat を渡せば内蔵テーブルに無い成分にも使える（既定は水の Antoine 式）。
        """
        if formula not in gas.index:
            raise ConstraintError(f"{gas.name!r} に成分 {formula!r} がありません")
        if psat is None:
            if formula != "H2O":
                raise ConstraintError(
                    f"{formula!r} の飽和蒸気圧は内蔵していません。psat=[Pa] を渡してください"
                )
            psat = antoine_water_psat(T)
        y_sat = psat / parse_pressure(P)
        if not 0.0 < y_sat < 1.0:
            raise ConstraintError(f"飽和モル分率が範囲外です: y_sat={y_sat:.4g}（T・P を確認）")
        ratio = y_sat / (1.0 - y_sat)
        i = gas.index[formula]

        def fn() -> np.ndarray:
            others = gas.molar_flows.sum() - gas.molar_flows[i]
            return np.atleast_1d(gas.molar_flows[i] - ratio * others)

        self.constraints.append(
            Constraint(fn, name or f"saturation[{formula}] {gas.name} @{T}°C")
        )

    def constrain_henry(
        self,
        gas: Stream,
        liquid: Stream,
        formulas: list[str],
        *,
        T: float,
        P: float | str,
        solvent: str = "H2O",
        henry: dict[str, float] | None = None,
        name: str | None = None,
    ) -> None:
        """Henry 則で気体の液相への溶解量を決める（**方程式 = len(formulas) 本**）。

            liquid(i) = x_i · liquid(solvent),  x_i = p_i / H_i,  p_i = gas(i)/Σgas · P

        H_i は Sander 2023 の内蔵テーブルから van't Hoff 式で温度補正して求める
        （`chemflow2.core.vle.HENRY_DATA`）。henry={成分: H[Pa]} で上書きできる。
        テーブルに無く上書きも無い成分は「溶けない」= liquid(i) = 0 とする。
        """
        P_pa = parse_pressure(P)
        overrides = henry or {}
        if solvent not in liquid.index:
            raise ConstraintError(f"{liquid.name!r} に溶媒 {solvent!r} がありません")
        j_solvent = liquid.index[solvent]

        for f in formulas:
            if f == solvent:
                raise ConstraintError(f"溶媒 {solvent!r} 自身は constrain_henry の対象にできません")
            if f not in liquid.index:
                raise ConstraintError(f"{liquid.name!r} に成分 {f!r} がありません")
            H = overrides.get(f, henry_pa(f, T))
            j = liquid.index[f]
            i = gas.index.get(f)

            def fn(j=j, i=i, H=H) -> np.ndarray:
                liq_i = liquid.molar_flows[j]
                if H is None:                       # 溶解データ無し → 溶けない
                    return np.atleast_1d(liq_i)
                gas_i = gas.molar_flows[i] if i is not None else 0.0
                gas_total = gas.molar_flows.sum()
                # 反復の途中でガスが空になっても壊れないようにする
                y_i = gas_i / gas_total if abs(gas_total) > 1e-12 else 0.0
                dissolved = y_i * (P_pa / H) * abs(liquid.molar_flows[j_solvent])
                return np.atleast_1d(liq_i - dissolved)

            self.constraints.append(
                Constraint(fn, name or f"henry[{f}] {gas.name}->{liquid.name} @{T}°C")
            )

    def constrain_carbon_activity(
        self,
        stream: Stream,
        value: float = 1.0,
        *,
        T: float | None = None,
        P: float | str | None = None,
        name: str | None = None,
    ) -> None:
        """ストリームの**炭素活量を指定値に固定する**（方程式 1 本）。

        炭素析出の判定は本来 ``a_C ≤ 1`` という不等式だが、chemflow2 は自由度が
        釣り合った等式系なので不等式は載らない。代わりに **``a_C = 1`` を課して
        「限界そのもの」を解く**。自由変数を 1 つ開けておけば、その変数の限界値が
        直接求まる（例: 水蒸気供給を未知にすれば、炭素析出しない最小の水蒸気量）。

            # H2/CO を保ったまま、炭素析出限界まで水蒸気を絞る条件を解く
            problem.constrain_carbon_activity(ReactOut, 1.0, T=900, P="1.04MPaG")

        評価するのは**反応器出口など析出が起きる場所**の組成・温度・圧力。
        T・P を省略すると ``stream.condition`` を使う。
        解いた後の検算には `chemflow2.core.carbon_activity.carbon_activity_of` を
        直接呼べばよい（制約にせず値だけ見る場合）。

        残差は **log(a_C) − log(value)** を使う。a_C は組成しだいで 1e±30 まで振れ、
        差の形だと反復の途中で overflow して解けない（分母の分圧が 0 に近づくため）。
        零点は同じで、条件数だけが桁違いに良くなる。
        """
        if value <= 0:
            raise ConstraintError(f"炭素活量の目標値は正の数にしてください: {value}")
        log_target = np.log(value)
        tiny = 1e-300

        def fn() -> np.ndarray:
            a = carbon_activity_of(stream, T=T, P=P)
            return np.atleast_1d(np.log(max(a, tiny)) - log_target)

        self.constraints.append(
            Constraint(fn, name or f"carbon_activity[{stream.name}]={value}")
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
