"""気液平衡の物性式: 水の飽和蒸気圧（Antoine）と気体の水への溶解度（Henry）。

`Problem.constrain_saturation` / `Problem.constrain_henry` が使う。装置ではなく物性式
なので core に置く（pressure.py・reaction.py と同じ扱い）。

**外部通信はしない。** Henry 定数は Sander 2023 の推奨値を内蔵テーブルに持ち、
van't Hoff 式で温度補正する。テーブルに無い成分は呼び出し側で値を渡すか、
溶解しないものとして扱う。
"""

from __future__ import annotations

import math

#: 水のモル体積 [m³/mol]
VM_WATER = 18.015e-6
#: Henry 定数の基準温度 [K]
T0_HENRY = 298.15


def antoine_water_psat(T_celsius: float) -> float:
    """Antoine 式で水の飽和蒸気圧を返す [Pa]。

    ``log10(P[mmHg]) = A − B / (C + T[°C])``、A/B/C は NIST の 1–100 °C 用の値。
    """
    A, B, C = 8.07131, 1730.63, 233.426
    return 10.0 ** (A - B / (C + T_celsius)) * 133.322   # mmHg → Pa


#: Henry 溶解度定数（298.15 K）と温度係数。
#: Hcp [mol/(m³·Pa)], Tderiv = d ln(Hcp)/d(1/T) [K]。
#: 出典: Sander, R., Atmos. Chem. Phys., 23, 10901–12440 (2023) の推奨値。
HENRY_DATA: dict[str, dict[str, float]] = {
    "H2":      {"Hcp": 7.8e-6,  "Tderiv": 500.0},
    "N2":      {"Hcp": 6.4e-6,  "Tderiv": 1300.0},
    "O2":      {"Hcp": 1.3e-5,  "Tderiv": 1500.0},
    "CO":      {"Hcp": 9.5e-6,  "Tderiv": 1300.0},
    "CO2":     {"Hcp": 3.3e-4,  "Tderiv": 2400.0},
    "CH4":     {"Hcp": 1.4e-5,  "Tderiv": 1600.0},
    "NH3":     {"Hcp": 5.9e-1,  "Tderiv": 4200.0},
    "H2S":     {"Hcp": 1.0e-3,  "Tderiv": 2100.0},
    "SO2":     {"Hcp": 1.2e-2,  "Tderiv": 2900.0},
    "CH3CHO":  {"Hcp": 1.3e-1,  "Tderiv": 5900.0},
    "CH3COOH": {"Hcp": 4.1e+3,  "Tderiv": 6300.0},
    "CH3OH":   {"Hcp": 2.2e+0,  "Tderiv": 5200.0},
    "C2H5OH":  {"Hcp": 1.9e+0,  "Tderiv": 6600.0},
    "HCHO":    {"Hcp": 3.2e+3,  "Tderiv": 6800.0},
    "HCOOH":   {"Hcp": 8.9e+3,  "Tderiv": 5700.0},
}

#: 極めて溶けにくい成分に使う Henry 定数 [Pa]（実質「溶けない」）。
HENRY_INSOLUBLE = 1.0e15


def henry_pa(formula: str, T_celsius: float) -> float | None:
    """化学式と温度から Henry 定数 H を返す [Pa]。テーブルに無ければ None。

    ``x_i = p_i / H_i`` の形で使う（x_i = 液相モル分率, p_i = 分圧）。
    Hcp は van't Hoff 式 ``Hcp(T) = Hcp(T0)·exp(Tderiv·(1/T − 1/T0))`` で温度補正し、
    ``H = 1 / (Hcp · Vm_water)`` で Pa 基準に直す。
    """
    data = HENRY_DATA.get(formula)
    if data is None:
        return None
    T_kelvin = T_celsius + 273.15
    hcp = data["Hcp"] * math.exp(data["Tderiv"] * (1.0 / T_kelvin - 1.0 / T0_HENRY))
    if hcp <= 0:
        return HENRY_INSOLUBLE
    return 1.0 / (hcp * VM_WATER)
