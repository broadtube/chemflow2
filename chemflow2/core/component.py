"""成分（示性式 → 分子量）。

molmass で示性式から分子量を自動計算し、結果をキャッシュする。
ユーザーが Component を直接触ることは基本的にない（Stream が内部で使う）。
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from chemflow2.core.errors import ComponentError

#: 理想気体のノルマルモル体積 [L/mol]（0°C, 1atm）
NORMAL_MOLAR_VOLUME = 22.414


@dataclass(frozen=True)
class Component:
    """化学成分。

    Attributes
    ----------
    formula : str
        示性式（"H2", "CO2", "CH3OCH3" 等）。
    mw : float
        分子量 [g/mol]。
    normal_volume : float
        ノルマルモル体積 [L/mol]。
    """

    formula: str
    mw: float
    normal_volume: float = NORMAL_MOLAR_VOLUME

    def __repr__(self) -> str:
        return f"Component({self.formula!r}, mw={self.mw:.4g})"


@lru_cache(maxsize=None)
def get_component(formula: str) -> Component:
    """示性式から Component を取得する（キャッシュあり）。"""
    try:
        import molmass
    except ImportError as e:  # pragma: no cover
        raise ComponentError(
            "molmass が必要です: pip install molmass"
        ) from e
    try:
        mw = float(molmass.Formula(formula).mass)
    except Exception as e:
        raise ComponentError(f"示性式を解決できません: {formula!r}") from e
    return Component(formula=formula, mw=mw)


@lru_cache(maxsize=None)
def get_composition(formula: str) -> dict[str, int]:
    """示性式の元素組成 {元素記号: 原子数} を返す（キャッシュあり）。"""
    try:
        import molmass
    except ImportError as e:  # pragma: no cover
        raise ComponentError("molmass が必要です: pip install molmass") from e
    try:
        comp = molmass.Formula(formula).composition()
    except Exception as e:
        raise ComponentError(f"示性式を解決できません: {formula!r}") from e
    return {symbol: item.count for symbol, item in comp.items()}
