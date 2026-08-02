"""炭素活量 a_C: 気相と平衡にある固体炭素（黒鉛）の活量。

基準状態は純黒鉛（T, 1 atm）で a = 1。**a_C > 1 なら気相が黒鉛に対して過飽和 =
炭素析出が熱力学的に有利**。改質器の運転条件を切るのに使う。

固体炭素を生む 3 反応それぞれについて、固体炭素の活量を未知数として平衡式を解く:

    B  2CO   ⇌ C(gr) + CO2    a_C = K_B · p_CO² / p_CO2
    M  CH4   ⇌ C(gr) + 2H2    a_C = K_M · p_CH4 / p_H2²
    R  CO+H2 ⇌ C(gr) + H2O    a_C = K_R · p_CO·p_H2 / p_H2O

3 式の差は改質系の反応そのもの（B−R = 水性ガスシフト, M−R = 水蒸気改質,
B−M = ドライ改質）なので、**気相が内部平衡していれば 3 つの a_C は厳密に一致する**。
GibbsReactor の出口はこれに当たり、値は一意。非平衡の気相（改質器入口、急冷後、
循環ガス）ではばらつくため、既定では**最大値**を取る（どれか 1 つでも 1 を超えれば
析出しうる、という保守側の扱い）。

⚠️ 分圧は **atm 基準**。K は無次元なので式中の p は本当は p/p°（p° = 1 atm）。
   Cantera の reference_pressure が 101325 Pa なのでそれに合わせている。

平衡定数は標準生成ギブズエネルギーから K = exp(−ΔG°/RT) で計算する。
気相は `gri30.yaml`、固体炭素は `graphite.yaml`（相 `C(gr)`）。
Cantera はオプション依存: `pip install chemflow2[gibbs]`

⚠️ この判定は**気相が完全平衡に達している前提**の工学的判定（"equilibrated gas"）。
   実際の改質管の入口付近など気相が平衡していない領域の局所的な析出リスクまでは
   保証しない。千代田化工の耐炭素析出触媒が効くのはまさにその速度論的な領域で、
   a_C ≤ 1 は「そうした触媒を前提にしない保守的な条件」にあたる。
   出典と検証は `references/PAPERS_INDEX.md` を参照。
"""

from __future__ import annotations

import math

from chemflow2.core.errors import CanteraError
from chemflow2.core.pressure import parse_pressure

#: 炭素析出を起こす 3 反応。(名前, 生成側の気体, 消費側の気体) は下の式で使う。
CARBON_ROUTES = ("Boudouard", "CH4_cracking", "CO_reduction")

_MECHANISM = "gri30.yaml"
_GRAPHITE = "graphite.yaml"
_SPECIES = ("CH4", "CO", "CO2", "H2", "H2O")

_cache: dict[float, tuple[float, float, float]] = {}


def equilibrium_constants(T_celsius: float) -> tuple[float, float, float]:
    """(K_B, K_M, K_R) を返す（1 atm 基準・無次元）。温度ごとにキャッシュする。"""
    key = round(float(T_celsius), 6)
    if key in _cache:
        return _cache[key]
    try:
        import cantera as ct
    except ImportError as e:  # pragma: no cover
        raise CanteraError(
            "炭素活量の計算には Cantera が必要です: pip install chemflow2[gibbs]"
        ) from e

    T_kelvin = float(T_celsius) + 273.15
    gas = ct.Solution(_MECHANISM)
    graphite = ct.Solution(_GRAPHITE)
    gas.TP = T_kelvin, gas.reference_pressure
    graphite.TP = T_kelvin, graphite.reference_pressure
    g = {s: gas.standard_gibbs_RT[gas.species_index(s)] for s in _SPECIES}
    g["C"] = graphite.standard_gibbs_RT[0]

    K = (math.exp(-((g["C"] + g["CO2"]) - 2 * g["CO"])),      # 2CO -> C + CO2
         math.exp(-((g["C"] + 2 * g["H2"]) - g["CH4"])),      # CH4 -> C + 2H2
         math.exp(-((g["C"] + g["H2O"]) - g["CO"] - g["H2"])))  # CO + H2 -> C + H2O
    _cache[key] = K
    return K


def carbon_activity(p_atm: dict[str, float], T_celsius: float, *,
                    per_route: bool = False):
    """分圧 [atm] の気相の炭素活量。

    per_route=True なら {経路名: a_C} を返す（平衡していない気相の診断用）。
    既定は 3 経路の最大値（保守側）。
    """
    K_B, K_M, K_R = equilibrium_constants(T_celsius)
    p = {s: max(float(p_atm.get(s, 0.0)), 0.0) for s in _SPECIES}
    tiny = 1e-300
    routes = {
        "Boudouard": K_B * p["CO"] ** 2 / max(p["CO2"], tiny),
        "CH4_cracking": K_M * p["CH4"] / max(p["H2"] ** 2, tiny),
        "CO_reduction": K_R * p["CO"] * p["H2"] / max(p["H2O"], tiny),
    }
    return routes if per_route else max(routes.values())


def carbon_activity_of(stream, T: float | None = None,
                       P: float | str | None = None, *, per_route: bool = False):
    """Stream の組成から炭素活量を計算する。

    T [°C] と P を省略すると stream.condition の値を使う。GibbsReactor の出口の
    ように気相が平衡していれば 3 経路は一致するので、値は一意。

        a = carbon_activity_of(ReactOut)          # 出口条件で評価
        a = carbon_activity_of(ReactOut, T=900)   # 温度だけ差し替え
    """
    T = T if T is not None else stream.condition.T
    P = P if P is not None else stream.condition.P
    if T is None or P is None:
        raise CanteraError("炭素活量の計算には T と P が必要です（引数か StreamCondition）")

    total = float(stream.molar_flows.sum())
    if total <= 0:
        raise CanteraError(f"{stream.name!r} の総流量が 0 なので分圧を定義できません")
    P_atm = parse_pressure(P) / 101325.0
    p_atm = {f: stream.flow_of(f) / total * P_atm for f in _SPECIES if f in stream.index}
    return carbon_activity(p_atm, T, per_route=per_route)
