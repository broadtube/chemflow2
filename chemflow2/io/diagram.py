"""フローシートの Mermaid 図出力。

ストリームは **連続線の中央に乗る番号バッジ**（丸数字 + 名前のエッジラベル）として描く。
装置は矩形ノード。Mermaid の構造上、線の途中に図形ノード（ひし形）を「乗せる」ことは
できないため、番号はエッジラベルとして中央に表示する（真のひし形が必要なら別途 SVG 描画）。

- どの装置の出口でもないストリーム = フィード（線の始点に小さな青丸）
- どの装置の入口でもないストリーム = プロダクト（線の終点に小さな緑丸）
- 循環ストリームは装置間の 1 本の連続線（中央に番号）としてループになる

番号 = Stream.order。番号と名前の対応は stream_table / Excel、
および export_mermaid が HTML に出す凡例表を参照。依存追加なし。
"""

from __future__ import annotations

import re

from chemflow2.core.stream import Stream


def _units(source) -> list:
    return source.units if hasattr(source, "units") else list(source)


def _uid(name: str | None) -> str:
    return "U_" + re.sub(r"\W", "_", name or "unit")


def _label(name: str | None) -> str:
    return (name or "").replace('"', "'").replace("|", "/")


def _badge(n) -> str:
    """番号を丸数字グリフに（1..20 は ①..⑳、範囲外は (n)）。"""
    if isinstance(n, int) and 1 <= n <= 20:
        return chr(0x2460 + n - 1)
    return f"({n})"


def _collect_streams(units: list) -> list[Stream]:
    streams: list[Stream] = []
    for u in units:
        for s in u.streams:
            if s not in streams:
                streams.append(s)
    return streams


def _numbered_streams(units: list) -> list[tuple[str, object, Stream]]:
    """(一意キー, 表示番号, Stream) のリスト。order 昇順。キーは端点ノード用に衝突しない連番。"""
    streams = _collect_streams(units)
    ordered = sorted(
        streams,
        key=lambda s: (s.order is None, s.order if s.order is not None else 0, s.name or ""),
    )
    result = []
    for i, s in enumerate(ordered):
        key = f"ST{i + 1}"
        number = s.order if s.order is not None else i + 1
        result.append((key, number, s))
    return result


def generate_mermaid(source, *, direction: str = "LR") -> str:
    """Problem もしくは units リストから Mermaid flowchart を生成する。"""
    units = _units(source)
    lines = [f"flowchart {direction}"]

    # 装置（矩形）
    for u in units:
        lines.append(f'    {_uid(u.name)}["{_label(u.name)}<br/><small>{type(u).__name__}</small>"]')

    # ストリーム（連続線 + 中央の番号バッジ）
    for key, number, s in _numbered_streams(units):
        producers = [u for u in units if s in u.outlets]
        consumers = [u for u in units if s in u.inlets]
        lbl = f"{_badge(number)} {_label(s.name)}".strip()
        if producers and consumers:
            for p in producers:
                for c in consumers:
                    lines.append(f"    {_uid(p.name)} -->|{lbl}| {_uid(c.name)}")
        elif consumers:  # フィード（生産者なし）: 始点に小さな丸
            src = f"IN_{key}"
            lines.append(f"    {src}(( )):::feed")
            for c in consumers:
                lines.append(f"    {src} -->|{lbl}| {_uid(c.name)}")
        elif producers:  # プロダクト（消費者なし）: 終点に小さな丸
            snk = f"OUT_{key}"
            lines.append(f"    {snk}(( )):::product")
            for p in producers:
                lines.append(f"    {_uid(p.name)} -->|{lbl}| {snk}")

    lines.append("    classDef feed fill:#1976d2,stroke:#0d47a1,color:#fff;")
    lines.append("    classDef product fill:#388e3c,stroke:#1b5e20,color:#fff;")
    return "\n".join(lines)


_HTML = """<!doctype html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<script src="https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.min.js"></script>
<script>mermaid.initialize({{ startOnLoad: true, theme: "default" }});</script>
<style>
  body{{font-family:sans-serif;margin:2rem}}
  h1{{font-size:1.2rem}}
  table{{border-collapse:collapse;margin-top:1rem;font-size:.9rem}}
  th,td{{border:1px solid #ccc;padding:2px 10px;text-align:left}}
  th{{background:#f5f5f5}}
  .legend{{margin-top:1.5rem}}
</style>
</head>
<body>
<h1>{title}</h1>
<pre class="mermaid">
{src}
</pre>
<div class="legend">
<h2 style="font-size:1rem">Stream legend</h2>
<table><tr><th>No.</th><th>Stream</th></tr>
{legend}
</table>
</div>
</body>
</html>
"""


def export_mermaid(source, path: str, *, title: str = "chemflow2 flowsheet", direction: str = "LR") -> str:
    """Mermaid 図を自己完結 HTML（番号↔名前の凡例表つき）として書き出す。

    生成した Mermaid ソースを返す。
    """
    units = _units(source)
    src = generate_mermaid(units, direction=direction)
    rows = "\n".join(
        f"<tr><td>{_badge(number)} {number}</td><td>{_label(s.name)}</td></tr>"
        for _key, number, s in _numbered_streams(units)
    )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(_HTML.format(title=title, src=src, legend=rows))
    return src
