"""圧力の文字列 → Pa(絶対圧) 変換。

数値はそのまま Pa(絶対圧)。文字列はゲージ圧("G"接尾)/絶対圧を解釈する。
    "2MPaG" → 2e6 + 101325
    "2MPa"  → 2e6
    "50kPaG"/"50kPa"/"3atm" などにも対応。
Gibbs/吸収など、状態量を要する装置でのみ使う。
"""

from __future__ import annotations

import re

ATM = 101325.0

_UNITS = {"mpa": 1e6, "kpa": 1e3, "pa": 1.0, "atm": ATM, "bar": 1e5}


def parse_pressure(P) -> float:
    """圧力を Pa(絶対圧) に変換する。"""
    if isinstance(P, (int, float)):
        return float(P)
    s = str(P).strip()
    m = re.match(r"^([0-9.]+)\s*([a-zA-Z]+?)(G?)$", s)
    if not m:
        raise ValueError(f"圧力を解釈できません: {P!r}")
    value, unit, gauge = float(m.group(1)), m.group(2).lower(), m.group(3)
    if unit not in _UNITS:
        raise ValueError(f"未知の圧力単位: {unit!r}（{P!r}）")
    pa = value * _UNITS[unit]
    if gauge:  # ゲージ圧 → 絶対圧
        pa += ATM
    return pa
