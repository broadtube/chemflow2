"""ストリーム結果の表出力（テキスト / CSV）。

basis を指定して mol/h・g/h・NL/h と各分率（mol%/wt%/vol%）を出せる。
複数指定すると basis ごとにセクションを重ねて出力する（既定は "mol"）。
"""

from __future__ import annotations

import numpy as np

from chemflow2.core.component import get_component
from chemflow2.core.stream import Stream

#: 絶対量 basis → 表示単位
_ABS_UNITS = {"mol": "mol/h", "mass": "g/h", "normal_volume": "NL/h"}
#: 分率 basis → 表示単位
_FRAC_UNITS = {"mole_frac": "mol%", "mass_frac": "wt%", "volume_frac": "vol%"}
_UNITS = {**_ABS_UNITS, **_FRAC_UNITS}


def _ordered(streams: list[Stream]) -> list[Stream]:
    return sorted(streams, key=lambda s: (s.order is None, s.order, s.name or ""))


def _all_formulas(streams: list[Stream]) -> list[str]:
    seen: dict[str, None] = {}
    for s in streams:
        for f in s.formulas:
            seen.setdefault(f, None)
    return list(seen)


def _total_mass(stream: Stream) -> float:
    mws = np.array([c.mw for c in stream.components])
    return float(np.sum(stream.molar_flows * mws))


def _as_bases(basis) -> list[str]:
    bases = [basis] if isinstance(basis, str) else list(basis)
    for b in bases:
        if b not in _UNITS:
            raise ValueError(f"未知の basis: {b!r}（{list(_UNITS)} のいずれか）")
    return bases


def _stream_values(stream: Stream, formulas: list[str], mw, nv, basis: str) -> np.ndarray:
    """指定 basis での成分値ベクトル（formulas 順、無い成分は 0）。"""
    mol = np.array([stream.flow_of(f) for f in formulas])
    if basis == "mol":
        return mol
    if basis == "mass":
        return mol * mw
    if basis == "normal_volume":
        return mol * nv
    # 分率（%）
    if basis == "mole_frac":
        base = mol
    elif basis == "mass_frac":
        base = mol * mw
    else:  # volume_frac
        base = mol * nv
    t = base.sum()
    return base / t * 100.0 if t else base * 0.0


def _formula_maps(formulas: list[str]):
    comps = [get_component(f) for f in formulas]
    return (np.array([c.mw for c in comps]), np.array([c.normal_volume for c in comps]))


def stream_table(streams: list[Stream], basis="mol") -> str:
    """成分 × ストリームの表を文字列で返す。

    basis: "mol"（既定）/"mass"/"normal_volume"/"mole_frac"/"mass_frac"/"volume_frac"、
    もしくはそれらのリスト（複数セクションを重ねる）。
    """
    bases = _as_bases(basis)
    streams = _ordered(streams)
    formulas = _all_formulas(streams)
    mw, nv = _formula_maps(formulas)
    headers = [s.name or f"S{i}" for i, s in enumerate(streams)]

    w = max([10] + [len(h) for h in headers]) + 2
    head = "component".ljust(10) + "".join(h.rjust(w) for h in headers)
    sep = "-" * len(head)
    lines = [head, sep]

    for b in bases:
        lines.append(f"[{_UNITS[b]}]")
        cols = [_stream_values(s, formulas, mw, nv, b) for s in streams]
        for i, f in enumerate(formulas):
            lines.append(f.ljust(10) + "".join(f"{col[i]:{w}.4g}" for col in cols))
        totals = [float(col.sum()) for col in cols]
        lines.append("total".ljust(10) + "".join(f"{t:{w}.4g}" for t in totals))
        lines.append(sep)

    # 質量閉包の確認（mass セクションが無いときのみ補助表示）
    if "mass" not in bases:
        masses = "total[g/h]".ljust(10) + "".join(f"{_total_mass(s):{w}.4g}" for s in streams)
        lines.append(masses)

    return "\n".join(lines)


def to_csv(streams: list[Stream], path: str, basis="mol") -> None:
    """成分 × ストリームの値を CSV 出力する（先頭列に basis の単位を付す）。"""
    import csv

    bases = _as_bases(basis)
    streams = _ordered(streams)
    formulas = _all_formulas(streams)
    mw, nv = _formula_maps(formulas)
    names = [s.name or f"S{i}" for i, s in enumerate(streams)]

    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["basis", "component"] + names)
        for b in bases:
            unit = _UNITS[b]
            cols = [_stream_values(s, formulas, mw, nv, b) for s in streams]
            for i, f in enumerate(formulas):
                writer.writerow([unit, f] + [f"{col[i]:.6g}" for col in cols])
            writer.writerow([unit, "total"] + [f"{float(col.sum()):.6g}" for col in cols])
