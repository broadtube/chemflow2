"""ストリーム結果の Excel(.xlsx) 出力。

元 chemflow は win32com で「起動中の Excel ブック」に書き込む Windows 専用方式だったが、
chemflow2 ではスタンドアロンの .xlsx ファイルを生成する（クロスプラットフォーム）。

openpyxl はオプション依存: `pip install chemflow2[excel]`
"""

from __future__ import annotations

import numpy as np

from chemflow2.core.stream import Stream
from chemflow2.io.table import _all_formulas, _ordered, _total_mass


def to_excel(streams: list[Stream], path: str, *, sheet: str = "Streams") -> None:
    """成分 × ストリームのモル流量表を .xlsx として書き出す。

    レイアウト:
        ヘッダ（ストリーム名 / T / P / Phase）
        成分行（示性式・MW・各ストリームのモル流量 [mol/h]）
        合計行（総モル流量 / 総質量流量 [g/h]）
    """
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font
    except ImportError as e:  # pragma: no cover
        raise ImportError("openpyxl が必要です: pip install chemflow2[excel]") from e

    streams = _ordered(streams)
    formulas = _all_formulas(streams)
    names = [s.name or f"S{i}" for i, s in enumerate(streams)]

    wb = Workbook()
    ws = wb.active
    ws.title = sheet
    bold = Font(bold=True)

    def _cond(s: Stream, attr: str):
        v = getattr(s.condition, attr, None)
        return "" if v is None else v

    # --- ヘッダ ---
    header = ["Component", "MW"] + names
    ws.append(header)
    for c in ws[ws.max_row]:
        c.font = bold
    ws.append(["T [°C]", ""] + [_cond(s, "T") for s in streams])
    ws.append(["P", ""] + [_cond(s, "P") for s in streams])
    ws.append(["Phase", ""] + [_cond(s, "phase") for s in streams])
    ws.append([])

    # --- 成分行（モル流量 [mol/h]）---
    unit_row = ["[mol/h]", ""] + ["" for _ in streams]
    ws.append(unit_row)
    for f in formulas:
        mw = next((c.mw for s in streams for c in s.components if c.formula == f), "")
        ws.append([f, round(mw, 4) if mw != "" else ""] + [s.flow_of(f) for s in streams])

    # --- 合計行 ---
    ws.append([])
    total_row = ws.max_row + 1
    ws.append(["total [mol/h]", ""] + [float(np.sum(s.molar_flows)) for s in streams])
    ws.append(["total mass [g/h]", ""] + [_total_mass(s) for s in streams])
    for r in (total_row, total_row + 1):
        ws.cell(row=r, column=1).font = bold

    # 列幅を軽く整える
    ws.column_dimensions["A"].width = 18
    for i in range(len(names)):
        ws.column_dimensions[chr(ord("C") + i)].width = 12

    wb.save(path)
