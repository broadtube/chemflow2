"""GibbsReactor: Cantera によるギブズ自由エネルギー最小化（化学平衡）。

転化率指定の Reactor と違い、指定 T・P での平衡組成そのものを解く。
残差方式に統合するため、残差 = 出口流量 - 平衡流量(入口) とする
（循環で入口が未知でも解ける。入口が既知なら 1 反復で収束）。

Cantera はオプション依存: `pip install chemflow2[gibbs]`
"""

from __future__ import annotations

import numpy as np

from chemflow2.core.errors import ChemflowError
from chemflow2.core.pressure import parse_pressure
from chemflow2.core.stream import Stream
from chemflow2.core.unit import Unit


class CanteraError(ChemflowError):
    """Cantera の読み込み・平衡計算に失敗した。"""


class GibbsReactor(Unit):
    """ギブズ平衡反応器。

        G1 = GibbsReactor(inlet=Feed, outlet=Out,
                          T=850, P="0.1MPa",
                          species=["CH4", "H2O", "CO", "CO2", "H2"], name="G1")

    Parameters
    ----------
    T : float | None
        温度 [°C]。None なら inlet.condition.T を使う。
    P : float | str | None
        圧力（Pa もしくは "2MPaG" 等）。None なら inlet.condition.P を使う。
    species : list[str]
        平衡計算に含める化学種（mechanism に存在する名前）。
        出口ストリームの成分はこの species と一致させること。
    mechanism : str
        Cantera の熱力学データファイル（既定 gri30.yaml）。

    残差: species ごとに ``出口 - 平衡流量 = 0``（長さ = len(species)）。
    """

    def __init__(
        self,
        inlet: Stream,
        outlet: Stream,
        *,
        species: list[str],
        T: float | None = None,
        P: float | str | None = None,
        mechanism: str = "gri30.yaml",
        name: str | None = None,
    ):
        self._inlet = inlet
        self._outlet = outlet
        self.species = list(species)
        self.name = name

        T = T if T is not None else inlet.condition.T
        P = P if P is not None else inlet.condition.P
        if T is None or P is None:
            raise CanteraError("GibbsReactor には T と P が必要です（引数か StreamCondition で指定）")
        self.T_kelvin = float(T) + 273.15
        self.P_pascal = parse_pressure(P)

        try:
            import cantera as ct
        except ImportError as e:  # pragma: no cover
            raise CanteraError("Cantera が必要です: pip install chemflow2[gibbs]") from e
        try:
            all_species = ct.Species.list_from_file(mechanism)
            selected = [s for s in all_species if s.name in self.species]
            found = {s.name for s in selected}
            missing = set(self.species) - found
            if missing:
                raise CanteraError(f"{mechanism} に存在しない化学種: {sorted(missing)}")
            self._gas = ct.Solution(thermo="ideal-gas", species=selected)
        except CanteraError:
            raise
        except Exception as e:  # pragma: no cover
            raise CanteraError(f"Cantera Solution の構築に失敗: {e}") from e

    @property
    def inlets(self) -> list[Stream]:
        return [self._inlet]

    @property
    def outlets(self) -> list[Stream]:
        return [self._outlet]

    def residuals(self) -> np.ndarray:
        import cantera as ct

        inlet_molar = np.array([self._inlet.flow_of(sp) for sp in self.species])
        outlet_molar = np.array([self._outlet.flow_of(sp) for sp in self.species])

        total_in = inlet_molar.sum()
        if total_in <= 0:
            return outlet_molar  # 入口ゼロ → 平衡もゼロ

        self._gas.TPX = self.T_kelvin, self.P_pascal, {
            sp: f for sp, f in zip(self.species, inlet_molar / total_in)
        }
        q = ct.Quantity(self._gas, moles=total_in / 1000.0)  # Cantera は kmol
        q.equilibrate("TP")

        eq_molar = np.array([q.moles * q.X[self._gas.species_index(sp)] * 1000.0 for sp in self.species])
        return outlet_molar - eq_molar
