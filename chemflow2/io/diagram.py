"""フローシートの Mermaid 図出力。

ユニットをノード、ストリームをエッジとして flowchart を生成する。
- どのユニットの出口でもないストリーム = フィード（外部入力ノード）
- どのユニットの入口でもないストリーム = プロダクト（外部出力ノード）
- 循環ストリームは「生産ユニット → 消費ユニット」のエッジとして自然にループになる

依存追加なし。generate_mermaid() は Mermaid ソース文字列、
export_mermaid() は CDN 版 mermaid.js を読み込む自己完結 HTML を書き出す。
"""

from __future__ import annotations

import re

from chemflow2.core.stream import Stream


def _units(source) -> list:
    return source.units if hasattr(source, "units") else list(source)


def _uid(name: str | None) -> str:
    return "U_" + re.sub(r"\W", "_", name or "unit")


def _sid(prefix: str, name: str | None) -> str:
    return f"{prefix}_" + re.sub(r"\W", "_", name or "s")


def _label(name: str | None) -> str:
    return (name or "").replace('"', "'")


def _collect_streams(units: list) -> list[Stream]:
    streams: list[Stream] = []
    for u in units:
        for s in u.streams:
            if s not in streams:
                streams.append(s)
    return streams


def generate_mermaid(source, *, direction: str = "LR") -> str:
    """Problem もしくは units リストから Mermaid flowchart ソースを生成する。"""
    units = _units(source)
    lines = [f"flowchart {direction}"]

    for u in units:
        lines.append(f'    {_uid(u.name)}["{_label(u.name)}<br/><small>{type(u).__name__}</small>"]')

    for s in _collect_streams(units):
        producers = [u for u in units if s in u.outlets]
        consumers = [u for u in units if s in u.inlets]
        label = _label(s.name)
        if producers and consumers:
            for p in producers:
                for c in consumers:
                    lines.append(f"    {_uid(p.name)} -->|{label}| {_uid(c.name)}")
        elif consumers:  # フィード（生産者なし）
            fid = _sid("feed", s.name)
            lines.append(f"    {fid}([{label}]):::feed")
            for c in consumers:
                lines.append(f"    {fid} --> {_uid(c.name)}")
        elif producers:  # プロダクト（消費者なし）
            pid = _sid("prod", s.name)
            lines.append(f"    {pid}([{label}]):::product")
            for p in producers:
                lines.append(f"    {_uid(p.name)} --> {pid}")

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
<style>body{{font-family:sans-serif;margin:2rem}} h1{{font-size:1.2rem}}</style>
</head>
<body>
<h1>{title}</h1>
<pre class="mermaid">
{src}
</pre>
</body>
</html>
"""


def export_mermaid(source, path: str, *, title: str = "chemflow2 flowsheet", direction: str = "LR") -> str:
    """Mermaid 図を自己完結 HTML として書き出す。生成した Mermaid ソースを返す。"""
    src = generate_mermaid(source, direction=direction)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(_HTML.format(title=title, src=src))
    return src
