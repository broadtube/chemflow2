"""KineticReactor: reaction_rate の速度論 PFR を chemflow2 の Unit として組み込む。

GibbsReactor と同型で、残差 = 出口流量 − pfr(入口) とする（循環でも解ける。
入口が既知なら 1 反復で収束）。**1 触媒床 = 1 Unit**。タンデム多床は
KineticReactor を複数直列に接続して表す（各床の出口 Stream を次床の入口にする）。

reaction_rate はオプション依存（sibling プロジェクト）:
    pip install -e ../reaction_rate       # ローカル editable
（PyPI 未公開のため chemflow2[...] には含めない）。

⚠️ 単位の橋渡し（本 Unit が境界で吸収する）:
  - 流量: chemflow2 は **mol/h**、reaction_rate.pfr は **mol/s** → flow_to_mol_s(既定 1/3600)。
  - 温度: StreamCondition.T は **°C** → pfr は **K**（+273.15）。
  - 圧力: parse_pressure は **Pa(絶対)** → pfr は **bar**（÷1e5）。
  - 種名: chemflow2 は実式（DME=CH3OCH3, MA=CH3COOCH3）、reaction_rate は短縮キー
    （"DME"/"MA"）。species_alias で相互変換（既定で両者を対応づけ）。
"""

from __future__ import annotations

import numpy as np

from chemflow2.core.errors import ChemflowError
from chemflow2.core.pressure import parse_pressure
from chemflow2.core.stream import Stream
from chemflow2.core.unit import Unit, component_union, flows_on

#: chemflow2 の実式 → reaction_rate の種キー（既定）。対応しない種はそのまま使う。
DEFAULT_SPECIES_ALIAS = {"CH3OCH3": "DME", "CH3COOCH3": "MA"}


class KineticError(ChemflowError):
    """reaction_rate の読み込み・PFR 計算に失敗した。"""


class KineticReactor(Unit):
    """速度論 PFR 反応器（reaction_rate バックエンド・1 触媒床 = 1 Unit）。

        # ハイブリッド床（合成 Cu/ZnO + 脱水 ZSM-5）
        R1 = KineticReactor(
            inlet=S1, outlet=S2,
            masses={"synthesis": 6.84e-4, "dehydration": 3.79e-4},  # [kg]
            models={"synthesis": "KOGAS", "dehydration": "ZSM5"},
            k_eq3="thermo", T=250, P="51bar", name="hybrid bed")
        # カルボニル化床（H-MOR）
        R2 = KineticReactor(
            inlet=S2, outlet=S3,
            masses={"carbonylation": 1.674e-3},
            models={"carbonylation": "DTU-Cheung2007-2"},
            T=250, P="51bar", name="MA bed")

    Parameters
    ----------
    masses : dict[str, float]
        {役割: 触媒質量[kg]}。役割 = 'synthesis'/'dehydration'/'carbonylation'。
    models : dict[str, str]
        {役割: モデル名}（例 dehydration:'ZSM5', carbonylation:'DTU-Cheung2007-2'）。
    k_eq3 : str
        脱水平衡定数 'thermo'(推奨)/'BL'/'KOGAS'(過大・非推奨)。
    acid_site_density : float
        mol Al/kg（DTU 系カルボニル化の per-mol-Al → per-kg 変換）。既定 1.43。
    T : float | None
        温度 [°C]。None なら inlet.condition.T。
    P : float | str | None
        圧力（Pa 数値 or "51bar"/"5MPaG" 等）。None なら inlet.condition.P。
    flow_to_mol_s : float
        ストリーム流量 → mol/s の係数（既定 1/3600 = ストリームが mol/h）。
    species_alias : dict[str, str] | None
        chemflow2 実式 → reaction_rate 種キーの対応（既定 DEFAULT_SPECIES_ALIAS）。
    n_points : int
        pfr の積分出力点数（既定 200）。

    残差: **出口成分ごとに `出口流量 − pfr(入口)流量 = 0`**（長さ = 出口成分数）。
    出口 Stream には反応網の全 7 種（CO, CO2, H2, H2O, CH3OH, CH3OCH3, CH3COOCH3）
    と不活性種を漏れなく含めること（GibbsReactor と同じ約束）。
    """

    def __init__(
        self,
        inlet: Stream,
        outlet: Stream,
        *,
        masses: dict[str, float],
        models: dict[str, str],
        k_eq3: str = "thermo",
        acid_site_density: float = 1.43,
        T: float | None = None,
        P: float | str | None = None,
        flow_to_mol_s: float = 1.0 / 3600.0,
        species_alias: dict[str, str] | None = None,
        n_points: int = 200,
        name: str | None = None,
    ):
        self._inlet = inlet
        self._outlet = outlet
        self.masses = dict(masses)
        self.models = dict(models)
        self.k_eq3 = k_eq3
        self.acid_site_density = float(acid_site_density)
        self.flow_to_mol_s = float(flow_to_mol_s)
        self.species_alias = dict(species_alias) if species_alias is not None else dict(DEFAULT_SPECIES_ALIAS)
        self.n_points = int(n_points)
        self.name = name

        T = T if T is not None else inlet.condition.T
        P = P if P is not None else inlet.condition.P
        if T is None or P is None:
            raise KineticError("KineticReactor には T と P が必要です（引数か StreamCondition で指定）")
        self.T_kelvin = float(T) + 273.15
        self.P_bar = parse_pressure(P) / 1.0e5  # Pa(絶対) → bar

        try:
            from reaction_rate.reactors import CatalystBed, pfr
        except ImportError as e:  # pragma: no cover
            raise KineticError(
                "reaction_rate が必要です: pip install -e ../reaction_rate"
            ) from e
        self._pfr = pfr
        self._bed = CatalystBed(self.masses, acid_site_density=self.acid_site_density)

    # ------------------------------------------------------------------ #
    @property
    def inlets(self) -> list[Stream]:
        return [self._inlet]

    @property
    def outlets(self) -> list[Stream]:
        return [self._outlet]

    def _to_rr(self, formula: str) -> str:
        """chemflow2 実式 → reaction_rate 種キー。"""
        return self.species_alias.get(formula, formula)

    def residuals(self) -> np.ndarray:
        out_formulas = self._outlet.formulas
        all_formulas = component_union([self._inlet, self._outlet])

        # pfr 入口: reaction_rate キーで mol/s（活性種は 0 でも含める）
        F_in = {self._to_rr(f): self._inlet.flow_of(f) * self.flow_to_mol_s for f in all_formulas}

        outlet_molar = flows_on(self._outlet, out_formulas)
        if sum(F_in.values()) <= 0:
            return outlet_molar  # 入口ゼロ → 出口ゼロ

        result = self._pfr(
            F_in, self.T_kelvin, self.P_bar, self._bed,
            models=self.models, k_eq3=self.k_eq3, n_points=self.n_points,
        )
        out = result.outlet()  # reaction_rate キー, mol/s
        predicted = np.array([out.get(self._to_rr(f), 0.0) / self.flow_to_mol_s for f in out_formulas])
        return outlet_molar - predicted
