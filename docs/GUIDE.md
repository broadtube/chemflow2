# chemflow2 使い方ガイド

全機能の使い方と最小の実行例をまとめる。各スニペットはそのまま実行できる
（`PYTHONPATH=.` でリポジトリ直下から実行するか、`pip install -e .` 後に import）。

- [1. 基本概念](#1-基本概念)
- [2. Stream（ストリーム）](#2-stream)
- [3. StreamCondition（状態量）](#3-streamcondition)
- [4. Reaction（反応）](#4-reaction)
- [5. 装置（Units）](#5-装置-units)
  - [Mixer](#mixer) / [Reactor](#reactor) / [Separator](#separator) / [Splitter](#splitter) / [GibbsReactor](#gibbsreactor)
- [6. Problem と制約](#6-problem-と制約)
- [7. Expr（制約 DSL）](#7-expr制約-dsl)
- [8. 出力（I/O）](#8-出力io)
- [9. 新しい装置を追加する](#9-新しい装置を追加する)
- [10. よくあるパターン](#10-よくあるパターン)

---

## 1. 基本概念

chemflow2 は **定常状態の物質収支** を解く。考え方は 3 つだけ。

1. **Stream の状態はモル流量ベクトルだけ**。質量・体積・分率はすべて導出。
2. **各装置は「満たすべき残差（= 0）」を返すだけ**。出口を直接計算しない。
3. **Problem が全変数をまとめ、連立方程式として解く**（`scipy.optimize`）。

ストリームには **固定**（`flows=` を与えた既知値）と **未知**（solve 対象）がある。

```python
from chemflow2 import Stream, Mixer, Problem

a = Stream(["N2", "H2"], flows={"N2": 1, "H2": 3})   # 固定（既知）
b = Stream(["N2", "H2"], flows={"N2": 2, "H2": 1})   # 固定
out = Stream(["N2", "H2"], name="out")               # 未知（solve 対象）

p = Problem(streams=[a, b, out], units=[Mixer([a, b], out)])
sol = p.solve()
print(out.flow_of("N2"), out.flow_of("H2"))   # 3.0 4.0
```

**自由度**（変数の数 = 方程式の数）が合わないと `SolveError`。事前確認:

```python
print(p.degrees_of_freedom())   # (変数, 方程式) 例: (2, 2)
```

---

## 2. Stream

```python
Stream(components, *, name=None, order=None, condition=None,
       flows=None, basis="mol", total=None, guess=None)
```

| 引数 | 説明 |
|------|------|
| `components` | 成分の示性式リスト（`molmass` が解決できる文字列） |
| `flows` | 与えると **固定ストリーム**。省略で **未知ストリーム** |
| `basis` | `flows` の単位。`"mol"`(既定)/`"mass"`/`"normal_volume"`/`"mole_frac"`/`"mass_frac"`/`"volume_frac"` |
| `total` | 分率系 basis のときの合計量（mol/mass/NL） |
| `guess` | 未知ストリームの初期推定（dict / list / ndarray） |
| `name` / `order` | 表示名・表示順 |
| `condition` | `StreamCondition`（T/P/相） |

```python
# いろいろな指定方法
Stream(["H2", "N2"], flows={"H2": 20, "N2": 60})                          # mol/h
Stream(["H2", "N2"], flows={"H2": 2.016, "N2": 28.0}, basis="mass")       # g/h → mol
Stream(["H2", "N2"], flows={"H2": 0.75, "N2": 0.25},
       basis="mole_frac", total=100)                                       # 組成 + 合計
Stream(["H2", "CO", "CO2"], name="Feed")                                   # 未知
Stream(["H2", "CO"], guess={"H2": 50, "CO": 20})                           # 未知 + 初期推定
```

**プロパティ / メソッド**（数値で欲しいものは `float()` か `flow_of`、制約に使うものは `*_expr` / `frac_of`）:

| 取得 | 数値（表示用） | Expr（制約用） |
|------|----------------|----------------|
| 成分流量 | `s.flow_of("CO")` | `s.flow_expr("CO")` |
| 成分モル分率 | — | `s.frac_of("CO")` |
| 総モル流量 | `float(s.total_flow)` | `s.total_flow` |
| 総質量流量 | `float(s.total_mass_flow)` | `s.total_mass_flow` |
| 総ノルマル体積流量 | `float(s.total_normal_volume_flow)` | `s.total_normal_volume_flow` |
| モル分率ベクトル | — | `s.mole_fractions` |

```python
s = Stream(["CO", "H2"], flows={"CO": 10, "H2": 30})
print(s.flow_of("CO"))            # 10.0
print(float(s.total_flow))        # 40.0
print(f"{s.total_mass_flow:.2f}") # 340.53   （Expr は f-string でも数値化できる）
```

---

## 3. StreamCondition

温度・圧力・相の **メタ情報**。v1 では物質収支には使わず、表示と GibbsReactor の入力条件に使う。

```python
from chemflow2 import StreamCondition

cond = StreamCondition(T=300, P="0.5MPaG", phase="gas")   # T[°C], P(数値Pa or 文字列), 相
Stream(["H2", "N2"], name="Feed", condition=cond)
```

圧力文字列は `"2MPaG"`(ゲージ)/`"2MPa"`(絶対)/`"3atm"`/`"50kPa"`/`"1bar"` に対応。単体でも使える:

```python
from chemflow2 import parse_pressure
parse_pressure("2MPaG")   # 2101325.0  (Pa 絶対圧)
```

---

## 4. Reaction

```python
Reaction(stoich, name=None, check=True, atol=1e-9)
```

`stoich` は係数の dict（反応物は負、生成物は正）。生成時に **原子収支を自動検査** する。

```python
from chemflow2 import Reaction

r = Reaction({"CO": -2, "H2": -4, "CH3OCH3": 1, "H2O": 1}, name="DME")

# 係数タイポ（原子が消える）は生成時にエラー
from chemflow2 import ReactionError
try:
    Reaction({"CO": -2, "H2": -4, "CH3OCH3": 1})   # H2O 抜け
except ReactionError as e:
    print(e)   # 原子収支が崩れています (O: -1, H: -2)...

# 検査を無効化したい場合
Reaction({"A2": -1, "A": 2}, check=False)

# 元素収支を確認（デバッグ用、0 が理想）
print(r.element_balance())   # {'C': 0.0, 'O': 0.0, 'H': 0.0}
```

---

## 5. 装置 (Units)

すべて共通で `inlet` / `outlet`（Stream オブジェクト）で結線する。循環は「同じ Stream を
複数装置で共有する」だけで表現できる。

### Mixer

複数入口 → 1 出口。残差: 成分ごとに `出口 - Σ入口 = 0`。

```python
from chemflow2 import Mixer
Mixer(inlet=[S1, S2], outlet=Mixed, name="M1")
```

### Reactor

転化率 + 選択率で複数反応を扱う。残差: `出口 - (入口 + Σ 進行度·係数) = 0`。

```python
Reactor(inlet, outlet, reactions, *, key_component, conversion,
        selectivities=None, name=None)
```

- `key_component` の総消費量 = `conversion × 入口の key 流量`
- それを `selectivities`（各反応への配分割合、**和 = 1**）で各反応に分ける
- 単一反応なら `selectivities` 省略で `[1.0]`

```python
from chemflow2 import Reaction, Reactor

# 単一反応
r = Reaction({"N2O4": -1, "NO2": 2})
Reactor(inlet, out, [r], key_component="N2O4", conversion=0.5)

# 複数反応（CO を key に、メタノール 70% / メタン化 30% へ配分）
r1 = Reaction({"CO": -1, "H2": -2, "CH3OH": 1}, name="methanol")
r2 = Reaction({"CO": -1, "H2": -3, "CH4": 1, "H2O": 1}, name="methanation")
Reactor(inlet, out, [r1, r2], key_component="CO",
        conversion=0.9, selectivities=[0.7, 0.3])
# → CO 転化率 90%、CH3OH = 0.7·(消費CO)、CH4 = 0.3·(消費CO)
```

### Separator

1 入口 → 複数出口の **一般分離ノード**。課すのは物質収支だけ（`入口 = Σ出口`）。
どの成分がどこへ行くかは **制約で指定**（[6 章](#6-problem-と制約)の `constrain_recovery`）。

```python
from chemflow2 import Separator
Sep1 = Separator(inlet=S3, outlet=[Gas, Liquid], name="Sep1")
problem.constrain_recovery(S3, Liquid, {"H2O": 1.0, "H2": 0.0, "CO": 0.0})  # 分配を指定
```

### Splitter

**組成そのまま比率分割**（分流器・パージ）。残差: 出口ごとに `出口 = 入口 × 比率`。
比率の和は 1（検査あり）。

```python
from chemflow2 import Splitter
Splitter(inlet=Rout, outlet=[Product, Recycle], ratios=[0.7, 0.3], name="SP1")
```

Separator との違い: Splitter は全出口が入口と **同一組成**。組成が変わる分離は Separator。

### GibbsReactor

指定 T・P での **平衡組成** を Cantera で解く（転化率では表現できない平衡計算）。
オプション依存: `pip install chemflow2[gibbs]`。

```python
GibbsReactor(inlet, outlet, *, species, T=None, P=None,
             mechanism="gri30.yaml", name=None)
```

- `species`: 平衡に含める化学種（`mechanism` に存在する名前）。出口の成分はこれに揃える
- `T`[°C] / `P`: 省略時は `inlet.condition` から取得

```python
from chemflow2 import GibbsReactor, Stream, StreamCondition, Problem

species = ["CH4", "H2O", "CO", "CO2", "H2"]
feed = Stream(species, flows={"CH4": 1, "H2O": 2, "CO": 0, "CO2": 0, "H2": 0},
              condition=StreamCondition(T=850, P="0.1MPa"))
out = Stream(species, name="eq")
G = GibbsReactor(feed, out, species=species)   # T/P は condition から
Problem([feed, out], [G]).solve()
print(out.flow_of("H2"))   # 平衡 H2 流量
```

---

## 6. Problem と制約

```python
Problem(streams, units, name=None)
```

すべてのストリームと装置を明示的に渡す。制約は 3 種類。

```python
# (a) 任意の等式: lhs == rhs（残差 = lhs - rhs）
problem.constrain(S1.total_flow, 165, name="S1 total flow")
problem.constrain(S1.total_mass_flow, S2.total_mass_flow)   # 両辺 Expr でも可

# (b) モル分率の指定
problem.constrain_fracs(S1, {"H2": 0.48, "CO": 0.24})       # CO2 は残り

# (c) 成分回収率（inlet の成分 c の frac が outlet へ）
problem.constrain_recovery(S3, Liquid, {"H2O": 1.0, "H2": 0.0})
```

**求解**:

```python
sol = problem.solve()                      # root（既定、複数手法を自動試行）
sol = problem.solve(bounds=(0, np.inf))    # 非負制約付き（least_squares）
sol = problem.solve(tol=1e-10)             # 収束判定の残差ノルム閾値
sol = problem.solve(method="hybr")         # 手法を明示（フォールバックしない）
```

**Solution**:

```python
print(sol.success)          # 収束したか
print(sol.message)          # ソルバのメッセージ（‖residual‖ を含む）
print(sol.nfev)             # 残差評価回数
sol.print_report()          # ストリーム表を表示
text = sol.report()         # 表を文字列で取得
```

`solve()` が収束しない場合や自由度が合わない場合は `SolveError`。

---

## 7. Expr（制約 DSL）

`S1.total_flow` などは即座に float を返さず **遅延評価の `Expr`** を返す。これにより
「その値の参照」を制約に渡せる。同時に `float()` / f-string で数値化もできる。

```python
problem.constrain(S1.total_flow, 165)      # solve 中に S1 の総流量を評価
print(f"{S1.total_flow:.2f}")              # solve 後は数値として表示

# 用意されたプロパティで足りないときは Expr を自作（成分流量ベクトルなど）
from chemflow2 import Expr
import numpy as np
problem.constrain(Expr(lambda: Recycle.molar_flows - 0.3 * Rout.molar_flows))
```

`Expr` は `- + * /` を最小限サポートするので `rout.total_flow * 0.3` のような式も書ける。

---

## 8. 出力 (I/O)

```python
from chemflow2 import stream_table, to_csv, to_excel, generate_mermaid, export_mermaid

print(stream_table(problem.streams))          # テキスト表（total[g/h] 行で質量閉包を確認）
to_csv(problem.streams, "streams.csv")        # CSV
to_excel(problem.streams, "streams.xlsx")     # Excel（要 openpyxl: pip install chemflow2[excel]）

print(generate_mermaid(problem))              # Mermaid ソース文字列（Markdown に貼れる）
export_mermaid(problem, "flow.html", title="My Flowsheet")  # 自己完結 HTML
```

**多 basis 出力**: `stream_table` / `to_csv` / `to_excel` は `basis` 引数を取る。
単一 or 複数を指定でき、複数指定すると basis ごとにセクションを重ねる。

| basis | 単位 |
|-------|------|
| `"mol"`（既定） | mol/h |
| `"mass"` | g/h |
| `"normal_volume"` | NL/h |
| `"mole_frac"` / `"mass_frac"` / `"volume_frac"` | mol% / wt% / vol% |

```python
# モル・質量・ノルマル体積・モル分率を一度に
print(stream_table(problem.streams, basis=["mol", "mass", "normal_volume", "mole_frac"]))
to_excel(problem.streams, "streams.xlsx", basis=["mol", "mass"])
```

`stream_table` / 出力はすべて `Stream.order` の昇順で並ぶ。

Mermaid はストリームを **連続線の中央に乗る番号バッジ**（丸数字 + 名前、番号 = `Stream.order`）、
装置を矩形で描く。フィード（青丸）/プロダクト（緑丸）を端点に置き、循環は 1 本の連続線になる。
番号と名前の対応は `export_mermaid` が HTML に出力する凡例表、または `stream_table` を参照。
`order` を付けておくと図の番号が意味を持つ。

```python
# 既定: 丸数字 + 名前（ポータブル）
export_mermaid(problem, "flow.html")
# PFD 風: 線の中央に番号入りの白いひし形（SVG をラベルに注入。このHTML専用）
export_mermaid(problem, "flow.html", style="diamond")
```

`style="diamond"` は Mermaid のエッジラベルにインライン SVG を入れて、線の中央に本物の
ひし形マーカーを乗せる。`securityLevel:'loose'` が必要なため `export_mermaid` の HTML で
のみ描画される。生ソースを GitHub 等に貼るなら既定の `"badge"`（丸数字）を使う。

---

## 9. 新しい装置を追加する

`Unit` を継承し `residuals()` を実装するだけ。`Problem(units=[...])` に渡せば求解に組み込まれる。

```python
# chemflow2/units/purge.py
import numpy as np
from chemflow2.core.unit import Unit, component_union, flows_on

class Purge(Unit):
    """出口 = 入口 × 比率（1 出口版の分流）。"""
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

ヘルパ:
- `component_union(streams)` … 複数ストリームの成分を出現順で統合
- `flows_on(stream, formulas)` … その成分順にモル流量ベクトルを整列（無い成分は 0）

`units/__init__.py` に import を足すと `from chemflow2 import Purge` で使える。

---

## 10. よくあるパターン

### 循環（リサイクル）

テアストリームを未知として宣言し、Mixer 入口と Splitter 出口で共有する。

```python
Recycle = Stream(comps, name="Recycle")            # テア（未知）
Mixer([Feed, Recycle], Mixed)                      # 入口で使う
Reactor(Mixed, Rout, [rxn], key_component=..., conversion=...)
Splitter(Rout, [Product, Recycle], ratios=[0.7, 0.3])   # 出口で同じオブジェクトを使う
Problem([Feed, Recycle, Mixed, Rout, Product], [...]).solve(bounds=(0, np.inf))
```

### フィードの組成 + 総流量を指定

```python
S1 = Stream(comps, name="Feed")                     # 未知
problem.constrain_fracs(S1, {"H2": 0.48, "CO": 0.24})
problem.constrain(S1.total_flow, 165)
```

### 完全分離（成分ごとに行き先を固定）

```python
Sep = Separator(S3, [Gas, Liquid])
problem.constrain_recovery(S3, Liquid, {"H2O": 1.0, "H2": 0.0, "CO": 0.0, "CO2": 0.0})
```

### 平衡反応器

```python
G = GibbsReactor(Feed, Out, species=["CH4","H2O","CO","CO2","H2"], T=850, P="0.1MPa")
```

---

## 実行できる完全な例

`examples/` にすべて揃っている:

| ファイル | 内容 |
|----------|------|
| `example_declarative.py` | 宣言的な書き方（Feed 組成指定 + 反応 + 分離） |
| `example_recycle.py` | 循環系（Separator + 生の Expr 制約） |
| `example_splitter.py` | 循環系（Splitter で簡潔に） |
| `example_recovery.py` | 成分回収率で分離を指定 |
| `example_gibbs.py` | Gibbs 平衡反応器（水蒸気メタン改質、要 cantera） |
| `example_diagram.py` | Mermaid フロー図の出力 |
| `example_multibasis.py` | 多 basis 出力（mol/h・g/h・NL/h・分率） |
| `example_excel.py` | Excel 出力（要 openpyxl） |
