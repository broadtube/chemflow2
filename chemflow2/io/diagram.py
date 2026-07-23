"""フローシートの Mermaid 図出力（PFD 式: ストリーム番号をひし形で表示）。

実プロセスフロー図の流儀に合わせ、ストリームを **番号入りのひし形ノード** として描く。
装置は矩形ノード、ストリームはその間のひし形（番号 = Stream.order）。

- どの装置の出口でもないストリーム = フィード（青）
- どの装置の入口でもないストリーム = プロダクト（緑）
- 循環ストリームはひし形を介した「生産装置 → 消費装置」のループになる

番号と名前の対応は stream_table / Excel、および export_mermaid が HTML に出す凡例表を参照。
依存追加なし。generate_mermaid() は Mermaid ソース文字列、
export_mermaid() は CDN 版 mermaid.js を読み込む自己完結 HTML（凡例付き）を書き出す。
"""

from __future__ import annotations

import re

from chemflow2.core.stream import Stream


def _units(source) -> list:
    return source.units if hasattr(source, "units") else list(source)


def _uid(name: str | None) -> str:
    return "U_" + re.sub(r"\W", "_", name or "unit")


def _label(name: str | None) -> str:
    return (name or "").replace('"', "'")


def _collect_streams(units: list) -> list[Stream]:
    streams: list[Stream] = []
    for u in units:
        for s in u.streams:
            if s not in streams:
                streams.append(s)
    return streams


def _numbered_streams(units: list) -> list[tuple[str, object, Stream]]:
    """(ノードID, 表示番号, Stream) のリスト。order 昇順。ノードIDは衝突しない連番。"""
    streams = _collect_streams(units)
    ordered = sorted(
        streams,
        key=lambda s: (s.order is None, s.order if s.order is not None else 0, s.name or ""),
    )
    result = []
    for i, s in enumerate(ordered):
        node_id = f"ST{i + 1}"
        number = s.order if s.order is not None else i + 1
        result.append((node_id, number, s))
    return result


def _stream_class(s: Stream, units: list) -> str:
    producers = [u for u in units if s in u.outlets]
    consumers = [u for u in units if s in u.inlets]
    if not producers:
        return "feed"
    if not consumers:
        return "product"
    return "stream"


def generate_mermaid(source, *, direction: str = "LR") -> str:
    """Problem もしくは units リストから PFD 式 Mermaid flowchart を生成する。"""
    units = _units(source)
    lines = [f"flowchart {direction}"]

    # 装置（矩形）
    for u in units:
        lines.append(f'    {_uid(u.name)}["{_label(u.name)}<br/><small>{type(u).__name__}</small>"]')

    numbered = _numbered_streams(units)

    # ストリーム（番号入りひし形）
    for node_id, number, s in numbered:
        lines.append(f'    {node_id}{{"{number}"}}:::{_stream_class(s, units)}')

    # 結線: 生産装置 → ひし形 → 消費装置
    for node_id, _number, s in numbered:
        for p in [u for u in units if s in u.outlets]:
            lines.append(f"    {_uid(p.name)} --> {node_id}")
        for c in [u for u in units if s in u.inlets]:
            lines.append(f"    {node_id} --> {_uid(c.name)}")

    lines.append("    classDef stream fill:#ffffff,stroke:#555,stroke-width:1px;")
    lines.append("    classDef feed fill:#e3f2fd,stroke:#1976d2,color:#0d47a1;")
    lines.append("    classDef product fill:#e8f5e9,stroke:#388e3c,color:#1b5e20;")
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
    """PFD 式 Mermaid 図を自己完結 HTML（番号↔名前の凡例表つき）として書き出す。

    生成した Mermaid ソースを返す。
    """
    units = _units(source)
    src = generate_mermaid(units, direction=direction)
    rows = "\n".join(
        f"<tr><td>{number}</td><td>{_label(s.name)}</td></tr>"
        for _node_id, number, s in _numbered_streams(units)
    )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(_HTML.format(title=title, src=src, legend=rows))
    return src
