"""サンプル: Mermaid フロー図を出力する。

solve 後の Problem（またはユニットのリスト）から、フローシートの Mermaid 図を
生成・出力する。ユニットがノード、ストリームがエッジになり、
フィード/プロダクト/循環が自動で描き分けられる。
"""

import os

import numpy as np

from chemflow2 import (
    Mixer,
    Problem,
    Reaction,
    Reactor,
    Splitter,
    Stream,
    export_mermaid,
    generate_mermaid,
)

components = ["N2O4", "NO2"]

Feed    = Stream(components, name="Feed",     order=1, flows={"N2O4": 10, "NO2": 0})
Recycle = Stream(components, name="Recycle",  order=5)
Mixed   = Stream(components, name="Mixed",    order=2)
Rout    = Stream(components, name="ReactOut", order=3)
Product = Stream(components, name="Product",  order=4)

rxn = Reaction({"N2O4": -1, "NO2": 2}, name="N2O4->2NO2")
M1  = Mixer([Feed, Recycle], Mixed, name="M1")
R1  = Reactor(Mixed, Rout, [rxn], key_component="N2O4", conversion=0.5, name="R1")
SP1 = Splitter(Rout, [Product, Recycle], ratios=[0.7, 0.3], name="SP1")

problem = Problem([Feed, Recycle, Mixed, Rout, Product], [M1, R1, SP1], name="Recycle")
problem.solve(bounds=(0, np.inf))

# 1) Mermaid ソース文字列を取得（そのまま Markdown 等に貼れる）
print("=== Mermaid source ===")
print(generate_mermaid(problem))

# 2) 自己完結 HTML として書き出す（ブラウザで開くと図が描画される）
out_dir = os.path.join(os.path.dirname(__file__), "output")
os.makedirs(out_dir, exist_ok=True)
path = os.path.join(out_dir, "flow_recycle.html")
export_mermaid(problem, path, title="Recycle Flowsheet")
print(f"\nHTML 出力: {path}")
