"""サンプル: 多 basis 出力（mol/h・g/h・NL/h・各分率）。

table / csv / excel はいずれも basis 引数を取り、単一 or 複数を指定できる。
複数指定すると basis ごとにセクションを重ねて出力する。
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
    StreamCondition,
    stream_table,
    to_csv,
)

components = ["N2O4", "NO2"]
cond = StreamCondition(T=25, P="0.1MPaG", phase="gas")

Feed    = Stream(components, name="1. Feed",     order=1, flows={"N2O4": 10, "NO2": 0}, condition=cond)
Recycle = Stream(components, name="5. Recycle",  order=5, condition=cond)
Mixed   = Stream(components, name="2. Mixed",    order=2, condition=cond)
Rout    = Stream(components, name="3. ReactOut", order=3, condition=cond)
Product = Stream(components, name="4. Product",  order=4, condition=cond)

rxn = Reaction({"N2O4": -1, "NO2": 2}, name="N2O4->2NO2")
M1  = Mixer([Feed, Recycle], Mixed, name="M1")
R1  = Reactor(Mixed, Rout, [rxn], key_component="N2O4", conversion=0.5, name="R1")
SP1 = Splitter(Rout, [Product, Recycle], ratios=[0.7, 0.3], name="SP1")

problem = Problem([Feed, Recycle, Mixed, Rout, Product], [M1, R1, SP1], name="Recycle")
problem.solve(bounds=(0, np.inf))

# モル流量・質量流量・ノルマル体積流量・モル分率を一度に表示
print(stream_table(problem.streams, basis=["mol", "mass", "normal_volume", "mole_frac"]))

# CSV も同様（先頭列に単位が付く）
out_dir = os.path.join(os.path.dirname(__file__), "output")
os.makedirs(out_dir, exist_ok=True)
path = os.path.join(out_dir, "streams_multibasis.csv")
to_csv(problem.streams, path, basis=["mol", "mass", "normal_volume"])
print(f"\nCSV 出力: {path}")

# Excel も basis 引数対応（要 openpyxl）
try:
    from chemflow2 import to_excel

    xlsx = os.path.join(out_dir, "streams_multibasis.xlsx")
    to_excel(problem.streams, xlsx, basis=["mol", "mass", "normal_volume", "mole_frac"])
    print(f"Excel 出力: {xlsx}")
except ImportError:
    print("Excel 出力スキップ（openpyxl 未導入）")
