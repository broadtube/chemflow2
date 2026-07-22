"""サンプル: 結果を Excel(.xlsx) に出力する。

solve 後のストリーム群を to_excel で .xlsx に書き出す。
成分ごとのモル流量、T/P/相、総モル・総質量が 1 シートにまとまる。

要 openpyxl: pip install chemflow2[excel]
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
    to_excel,
)

components = ["N2O4", "NO2"]
cond = StreamCondition(T=25, P="0.1MPaG", phase="gas")

Feed    = Stream(components, name="1. Feed",     order=1, flows={"N2O4": 10, "NO2": 0}, condition=cond)
Recycle = Stream(components, name="5. Recycle",  order=5, condition=cond)
Mixed   = Stream(components, name="2. Mixed",    order=2, condition=cond)
Rout    = Stream(components, name="3. ReactOut", order=3, condition=StreamCondition(T=120, P="0.1MPaG", phase="gas"))
Product = Stream(components, name="4. Product",  order=4, condition=cond)

rxn = Reaction({"N2O4": -1, "NO2": 2}, name="N2O4->2NO2")
M1  = Mixer([Feed, Recycle], Mixed, name="M1")
R1  = Reactor(Mixed, Rout, [rxn], key_component="N2O4", conversion=0.5, name="R1")
SP1 = Splitter(Rout, [Product, Recycle], ratios=[0.7, 0.3], name="SP1")

problem = Problem([Feed, Recycle, Mixed, Rout, Product], [M1, R1, SP1], name="Recycle")
problem.solve(bounds=(0, np.inf))

out_dir = os.path.join(os.path.dirname(__file__), "output")
os.makedirs(out_dir, exist_ok=True)
path = os.path.join(out_dir, "streams.xlsx")
to_excel(problem.streams, path, sheet="Streams")
print(f"Excel 出力: {path}")
