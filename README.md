# chemflow2

シンプルな定常状態プロセスシミュレーター。
各装置が「満たすべき残差（= 0）」を返し、`Problem` が全体を連立方程式として
`scipy.optimize` で解く（既存 chemflow の核だけを継承し、周辺機能を削ぎ落とした版）。

**➡ 全機能の使い方と実行例は [docs/GUIDE.md](docs/GUIDE.md) を参照。**

## 設計思想

- **状態はモル流量だけ** — 質量・体積・分率はすべて導出。
- **装置は残差を返すだけ** — 出口を直接計算しない。循環系も 1 つの連立方程式として解ける。
- **グローバル状態なし** — `Problem` にストリームとユニットを明示的に渡す宣言的スタイル。
- **循環は Stream 共有で表現** — テアストリームを未知として宣言し、複数ユニットで同じオブジェクトを使うだけ。

## フォルダ構成

```
chemflow2/
├── core/         不変の核（概念と求解）
│   ├── stream.py      Stream, StreamCondition
│   ├── reaction.py    Reaction
│   ├── unit.py        Unit 基底クラス（拡張ポイント）+ 物質収支ヘルパー
│   ├── problem.py     Problem, Constraint, Solution（求解）
│   ├── expr.py        Expr（遅延評価: constrain(S1.total_flow, 165) を可能にする）
│   ├── component.py   示性式 → 分子量（molmass）
│   └── errors.py
├── units/        差し替え・追加可能な装置
│   ├── mixer.py       Mixer
│   ├── reactor.py     Reactor（転化率 + 選択率 + 複数反応）
│   ├── separator.py   Separator（物質収支のみ。分配は制約で指定）
│   ├── splitter.py    Splitter（組成そのまま比率分割。分流・パージ）
│   └── gibbs.py       GibbsReactor（Cantera 平衡。オプション依存）
└── io/           出力
    ├── table.py       テキスト表 / CSV
    ├── excel.py       Excel(.xlsx) 出力（オプション依存）
    └── diagram.py     Mermaid フロー図（連続線 + 番号バッジ / 凡例つき HTML）
```

### Separator と Splitter の使い分け

| | 課す残差 | 組成 | 用途 |
|---|---|---|---|
| `Splitter(inlet, [out...], ratios=[...])` | 出口 = 入口 × 比率 | 入口と同一 | 分流・パージ・循環分岐 |
| `Separator(inlet, [out...])` | 入口 = Σ出口（収支のみ） | 変わる | 分離塔・気液分離（分配は制約で指定） |

## 使い方

```python
import numpy as np
from chemflow2 import Stream, StreamCondition, Reaction, Reactor, Separator, Problem

components = ["H2", "CO", "CO2", "CH3OCH3", "H2O"]

S1 = Stream(components, name="1. Feed", order=1,
            condition=StreamCondition(T=300, P=1, phase="gas"))
S3 = Stream(components, name="3. ReactOut", order=3)
S4 = Stream(components, name="4. Gas", order=4)
S5 = Stream(components, name="5. Liquid", order=5)

rxn1 = Reaction(stoich={"CO": -2, "H2": -4, "CH3OCH3": 1, "H2O": 1}, name="Tandem")
R1   = Reactor(inlet=S1, outlet=S3, reactions=[rxn1],
               key_component="CO", conversion=0.8, selectivities=[1], name="R1")
Sep1 = Separator(inlet=S3, outlet=[S4, S5], name="Sep1")

problem = Problem(streams=[S1, S3, S4, S5], units=[R1, Sep1])
problem.constrain_fracs(S1, {"H2": 0.48, "CO": 0.24})
problem.constrain(S1.total_flow, 165, name="S1 total flow specified")
# ... 分離の指定など

sol = problem.solve(bounds=(0, np.inf))
sol.print_report()
```

実行例:

```bash
PYTHONPATH=. python3 examples/example_declarative.py   # 宣言的な書き方
PYTHONPATH=. python3 examples/example_recycle.py       # 循環系（Separator + 制約）
PYTHONPATH=. python3 examples/example_splitter.py      # 循環系（Splitter で簡潔に）
PYTHONPATH=. python3 examples/example_diagram.py       # Mermaid フロー図の出力
PYTHONPATH=. python3 examples/example_recovery.py      # 成分回収率で分離を指定
PYTHONPATH=. python3 examples/example_gibbs.py         # Gibbs 平衡反応器（要 cantera）
PYTHONPATH=. python3 examples/example_excel.py         # Excel 出力（要 openpyxl）
PYTHONPATH=. python3 examples/example_multibasis.py    # 多 basis 出力
PYTHONPATH=. python3 -m pytest tests/ -q
```

### 多 basis 出力

`stream_table` / `to_csv` / `to_excel` は `basis` 引数を取り、mol/h・g/h・NL/h・各分率
（mol%/wt%/vol%）を出せる。複数指定で basis ごとにセクションを重ねる（既定 `"mol"`）。

```python
from chemflow2 import stream_table
print(stream_table(problem.streams, basis=["mol", "mass", "normal_volume", "mole_frac"]))
```

### Excel 出力（オプション）

`pip install chemflow2[excel]`。成分ごとの流量・T/P/相・合計を 1 シートに出力する。`basis` 対応。

```python
from chemflow2 import to_excel
to_excel(problem.streams, "streams.xlsx", sheet="Streams", basis=["mol", "mass"])
```

### 分離の指定（回収率）

`Separator` は物質収支だけを課すノード。分配は `constrain_recovery` で指定する。

```python
Sep1 = Separator(inlet=S3, outlet=[Gas, Liquid])
# H2O は全量 Liquid、それ以外は全量 Gas（= Liquid への回収率 0）
problem.constrain_recovery(S3, Liquid, {"H2O": 1.0, "H2": 0.0, "CO": 0.0})
```

### Gibbs 平衡反応器（オプション）

転化率ではなく、指定 T・P での平衡組成を Cantera が解く。`pip install chemflow2[gibbs]`

```python
from chemflow2 import GibbsReactor
G1 = GibbsReactor(inlet=Feed, outlet=Out, species=["CH4","H2O","CO","CO2","H2"],
                  T=850, P="0.1MPa")   # T[°C] / P は "2MPaG" 等も可
```

### フロー図の出力

ストリームを **連続線の中央に乗る番号バッジ**（丸数字 + 名前、番号 = `Stream.order`）、
装置を矩形で描く。フィード（青丸）/プロダクト（緑丸）を端点に置き、循環は 1 本の連続線になる。
番号↔名前の対応は HTML の凡例表と `stream_table` を参照。

> Mermaid の構造上、線の途中に図形ノード（ひし形）を「乗せる」ことはできないため、
> 番号はエッジラベルとして中央に表示する。真のひし形マーカーが必要なら SVG 直描画が別途必要。

```python
from chemflow2 import generate_mermaid, export_mermaid

print(generate_mermaid(problem))                       # Mermaid ソース文字列
export_mermaid(problem, "flow.html", title="My Flow")  # 凡例表つき HTML

# 線の中央に番号入りの白いひし形（SVG）を乗せる PFD 風スタイル
export_mermaid(problem, "flow.html", style="diamond")  # このHTML専用（securityLevel:'loose'）
```

`style="diamond"` は実 PFD のように番号をひし形マーカーで線上に描く。Mermaid のラベルに
インライン SVG を注入するため `securityLevel:'loose'` が要り、`export_mermaid` の HTML で
のみ描画される（生ソースを GitHub 等に貼る用途は既定の `"badge"` を使う）。

## 新しい装置を足す

`units/` にファイルを追加し、`Unit` を継承して `residuals()` を書くだけ。

```python
# chemflow2/units/purge.py
import numpy as np
from chemflow2.core.unit import Unit, component_union, flows_on

class Purge(Unit):
    def __init__(self, inlet, outlet, ratio, name=None):
        self._inlet, self._outlet, self.ratio, self.name = inlet, outlet, ratio, name

    @property
    def inlets(self):  return [self._inlet]
    @property
    def outlets(self): return [self._outlet]

    def residuals(self):
        f = component_union([self._inlet, self._outlet])
        return flows_on(self._outlet, f) - self.ratio * flows_on(self._inlet, f)
```

`units/__init__.py` に import を足せば `from chemflow2 import Purge` で使える。
`Problem(units=[...])` に渡すだけで求解に組み込まれる。

## 物質収支ツールとしての安全網

- **原子収支の検査**: `Reaction` は生成時に `Σ_i ν_i·(元素ベクトル)_i == 0` を検査し、
  崩れていれば `ReactionError`。係数タイポで原子が生成・消滅する「閉じない MB」を防ぐ
  （`Reaction(..., check=False)` で無効化可）。
- **選択率の和**: `Reactor` は `Σselectivities == 1` を検査。指定 conversion と実効 conversion のズレを防ぐ。
- **質量閉包の可視化**: レポート最下段に各ストリームの総質量行を出力。`in = out` を目視確認できる。
- **収束判定**: `least_squares` は solver の早期停止フラグではなく **最終残差ノルム** `‖residual‖ < tol` で判定。
- **初期推定**: 大規模・循環で all-ones が効かない場合は `Stream(..., guess={...})` で初期値を渡せる。

## 既存 chemflow / spec との関係

`.kiro/specs/chemflow-intuitive-api`（既存 chemflow の演算子オーバーロード版 API リデザイン）とは**別物**。
chemflow2 は演算子マジック・グローバル状態を持たず、明示的な `Problem` に寄せた再設計。
将来的に reaction_rate 側と連携する場合、「reaction_rate の PFR で出した転化率を `Reactor(conversion=...)` に渡す」か、
「reaction_rate を呼ぶ独自 Unit を `units/` に足す」のどちらもこの設計に自然に乗る。

## v1 の範囲（意図的に外したもの）

温度・圧力・相は `StreamCondition` にメタ情報として保持するだけで、物質収支には未使用
（Gibbs 反応器の入力条件としてのみ使う）。
エネルギー収支・相平衡・多段吸収(Henry)・reactflow 出力・JSON 復元は
将来の拡張ポイント（`units/` と `io/` に足す）として空けてある。
```
