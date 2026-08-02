"""chemflow2 の例外階層。"""

from __future__ import annotations


class ChemflowError(Exception):
    """chemflow2 の基底例外。"""


class ComponentError(ChemflowError):
    """示性式の解決・分子量計算に失敗した。"""


class BasisError(ChemflowError):
    """basis 指定と入力の組み合わせが不正。"""


class ReactionError(ChemflowError):
    """反応の化学量論が原子保存を満たさない等、反応定義が不正。"""


class ConstraintError(ChemflowError):
    """制約式の次元・内容が不正。"""


class SolveError(ChemflowError):
    """連立方程式が収束しない・自由度が合わない。"""


class CanteraError(ChemflowError):
    """Cantera の読み込み・平衡計算に失敗した。

    GibbsReactor（units/gibbs.py）と炭素活量（core/carbon_activity.py）の両方が
    使うため、単体の置き場ではなくここに置く。
    """
