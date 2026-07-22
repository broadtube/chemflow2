"""ストリーム結果の表出力（テキスト / CSV）。

v1 は最小限。mermaid / reactflow / excel 等はここに足していく拡張ポイント。
"""

from __future__ import annotations

import numpy as np

from chemflow2.core.stream import Stream


def _ordered(streams: list[Stream]) -> list[Stream]:
    return sorted(streams, key=lambda s: (s.order is None, s.order, s.name or ""))


def _all_formulas(streams: list[Stream]) -> list[str]:
    seen: dict[str, None] = {}
    for s in streams:
        for f in s.formulas:
            seen.setdefault(f, None)
    return list(seen)


def stream_table(streams: list[Stream]) -> str:
    """成分 × ストリームのモル流量表を文字列で返す。"""
    streams = _ordered(streams)
    formulas = _all_formulas(streams)
    headers = [s.name or f"S{i}" for i, s in enumerate(streams)]

    w = max([8] + [len(h) for h in headers]) + 2
    lines = []
    head = "component".ljust(10) + "".join(h.rjust(w) for h in headers)
    lines.append(head)
    lines.append("-" * len(head))
    for f in formulas:
        row = f.ljust(10) + "".join(f"{s.flow_of(f):{w}.4g}" for s in streams)
        lines.append(row)
    lines.append("-" * len(head))
    totals = "total".ljust(10) + "".join(f"{float(np.sum(s.molar_flows)):{w}.4g}" for s in streams)
    lines.append(totals)
    # 質量閉包の確認用（総質量 = Σ molar * MW）
    masses = "total_mass".ljust(10) + "".join(f"{_total_mass(s):{w}.4g}" for s in streams)
    lines.append(masses)
    return "\n".join(lines)


def _total_mass(stream: Stream) -> float:
    mws = np.array([c.mw for c in stream.components])
    return float(np.sum(stream.molar_flows * mws))


def to_csv(streams: list[Stream], path: str) -> None:
    """成分 × ストリームのモル流量を CSV 出力する。"""
    import csv

    streams = _ordered(streams)
    formulas = _all_formulas(streams)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["component"] + [s.name or f"S{i}" for i, s in enumerate(streams)])
        for f in formulas:
            writer.writerow([f] + [f"{s.flow_of(f):.6g}" for s in streams])
        writer.writerow(["total"] + [f"{float(np.sum(s.molar_flows)):.6g}" for s in streams])
