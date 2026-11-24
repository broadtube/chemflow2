# 改質器 + CO2 除去比較スクリプトの初期条件・拘束条件

対象:

- `examples/example_reformer_co2removal.py`（後工程 plant3）
- `examples/example_reformer_co2removal_plant4.py`（後工程 plant4）

どちらも **2 段階で解く**。

1. 改質器フローシート（`example_reformer.py`）を単独で解く
2. その出口ガス（PlantFeed = DryGas、`a_upstream` では CO2 除去後）を後工程の新鮮供給
   `1. Steam1` として**固定**し、後工程の循環ループを解く

苛性ソーダ吸収塔の薬剤収支（`example_co2_removal.caustic_scrubber`）は、主フローシートの
解を入力にした**後処理**で、循環計算には影響しない。

---

## 0. コマンドライン既定値

| 引数 | plant3 版 | plant4 版 | 意味 |
|---|---|---|---|
| `--reform` | `A` | `A` | 改質条件（A〜D） |
| `--cases` | `baseline,a_upstream,b_recycle` | `baseline,b_recycle,c_interstage` | CO2 除去位置 |
| `--eta` | `0.90` | `0.90` | CO2 除去率（カンマ区切りで複数値を順に計算） |
| `--purge` | `0.05` | `0.05` | 凝縮器ガスのパージ率 |
| `--sv-mode` | （なし） | **`plant3`** | 触媒体積の決め方（§3.3） |
| `--synthesis-model` | `VBF` | `VBF` | R1 のメタノール合成速度論（`KOGAS` に切替可） |
| `--solve-tol` | `1e-6` | `1e-4` | 内側求解の残差ノルム合格判定（絶対値） |
| `--retry-tol` | `1e-5` | `1e-3` | 1 回目が不合格のとき 1 回だけ緩めて再試行 |
| `--solver-tols` | `1e-12` | `1e-12` | least_squares の ftol / xtol / gtol |
| `--x-scale` | `jac` | `jac` | least_squares の変数スケール |
| `--outdir` | （なし。`examples/output/reformer/` 固定） | `examples/output/reformer/` | 出力先 |

> plant4 の `--solve-tol` が緩いのは、PFR の積分誤差（rtol 1e-8）で残差ノルムが約 3e-6 で
> 頭打ちになるため。
>
> `example_plant4.py` 単体の既定は `SV_MODE_DEFAULT = "per-bed"` だが、
> **このスクリプトのコマンドライン既定は `plant3`**。

### CO2 除去位置（ケース）

| ケース | 塔 | 位置 | plant3 版 | plant4 版 |
|---|---|---|---|---|
| `baseline` | — | 除去なし | ○ | ○ |
| `a_upstream` | T-100 | 改質器 DryGas → PlantFeed（後工程の外） | ○ | ○（既定外） |
| `b_recycle` | T-101 | 凝縮器ガス → パージ分岐の手前 | ○ | ○（Cond2Gas の後） |
| `c_interstage` | T-102 | Cond1Gas → R2（床の間） | — | ○ |

塔は主フローシートでは **Separator 1 個 + 「CO2 の回収率 = η」の制約 1 本**で表す。
CO2Vent は CO2 だけの 1 成分ストリームなので、他成分の回収率 0 は書かない
（書くと変数の無い式が増えて自由度が崩れる）。

---

## 1. 改質器（`example_reformer.py`）

### 1.1 固定値

| 項目 | 値 |
|---|---|
| 平衡種（GibbsReactor, Cantera） | CO2, CH4, H2O, CO, H2 |
| 凝縮器温度 | 25 ℃（改質器と同圧） |
| 標準状態モル体積 | 22.414 NL/mol |
| 炭素活量の設計値 | a_C = 0.90 |
| 切替圧力 `P_CROSSOVER` | 0.417 MPaG |

| ケース | T [℃] | P [MPaG] | 供給 |
|---|---|---|---|
| A | 835 | 1.41 | **実測値で全固定** CH4 23.2 / CO2 52.0 / H2O 32.7 / H2 5.1 NL/h（H2/CO は結果 1.067） |
| B | 875 | 1.41 | 未知（下の拘束で決まる） |
| C | 875 | 0.417 | 同上 |
| D | 875 | 0.1 | 同上 |

### 1.2 拘束条件

全ケース共通:

- Mixer（M0）、GibbsReactor（G1, 平衡）、凝縮器 Separator の物質収支
- DryGas の H2O は 25 ℃ 飽和（Antoine）
- 凝縮水への H2 / CO / CO2 / CH4 の溶解は Henry 則
- `a_upstream` のみ: DryGas の CO2 を η だけ T-100 で除去

B〜D のみ（4 本追加して自由度を閉じる）:

1. 供給合計 = 113 NL/h
2. 供給 CH4 : H2 = 78.9 : 21.1
3. DryGas の H2/CO = 4/3（3CO + 4H2 → MA + H2O の量論）
4. 水蒸気量を決める 1 本。**圧力から自動で選ぶ**（`binding_mode`）
   - P > 0.417 MPaG: CH4 転化率 = 0.88（A の 87.93% をわずかに上回る）
   - P ≤ 0.417 MPaG: 炭素活量 a_C = 0.90

自由度は A が変数 19 / 方程式 19、B〜D が 23 / 23（塔を付けると +6）。

### 1.3 初期推定と解き方

| ストリーム | 初期推定 [mol/h] |
|---|---|
| Feed（B〜D） | CH4 1.9 / CO2 1.8 / H2O 0.8 / H2 0.5 |
| ReactOut | CO 2.5 / H2 3.3 / CO2 0.8 / CH4 0.2 / H2O 0.6 |
| DryGas | CO 2.5 / H2 3.3 / CO2 0.8 / CH4 0.2 / H2O 0.02 |
| PlantFeed（塔あり） | CO 2.5 / H2 3.3 / CO2 0.1 / CH4 0.2 / H2O 0.02 |
| その他 | 全成分 1.0（Stream の既定） |

- **B〜D は連続法**: 1.41 → 1.0 → 0.7 → 0.5 → 0.417 → 0.3 → 0.2 → 0.15 → 0.1 MPaG と
  高圧から降り、各段の解を次段の初期値にする。各段の拘束 4 はその段の圧力で選ぶ。
  CH4 < 0.1 mol/h の縮退解は初期値に使わない。途中段の失敗は無視して進む。
- **CO2 除去塔は最後に足す**: 塔なしの解（改質器部分の厳密解）を初期値にして解き直す。
- ソルバー:
  - A（塔なし）: `root`、tol 1e-7
  - B〜D と塔ありの A: `least_squares`（流量 0〜1000、tol 1e-7、ftol/xtol/gtol 1e-14、
    x_scale `jac`）。失敗したら `root` で解き直し、負の流量（< −1e-9）が無ければ採用。

---

## 2. 苛性ソーダ吸収塔（後処理サブフローシート）

`example_co2_removal.caustic_scrubber`。両スクリプト共通。

| 項目 | 値 |
|---|---|
| 構成 | Mixer（気液接触）→ Reactor（化学吸収）→ Separator（分液） |
| 反応 | CO2 + 2NaOH → Na2CO3 + H2O（原子収支は Reaction が検査） |
| 条件 | 40 ℃・5 MPaG |
| 供給ガス | 主フローシートで解いた塔入口ガスを**固定値**として渡す |
| NaOH 供給 | 2 × η × CO2 × 過剰率 1.2 |
| 苛性液濃度 | 20 wt% NaOH aq |
| 拘束 | Reactor の CO2 転化率 = η、H2O は全量液側 |

---

## 3. 後工程

### 3.1 共通の固定値

| 項目 | 値 |
|---|---|
| 新鮮供給 `1. Steam1` | 改質器の PlantFeed（成分ごとに mol/h で固定） |
| ガス条件 | 250 ℃・5 MPaG（= 51.0 bar）、等温 |
| 液条件 | 40 ℃・5 MPaG |
| R1 ハイブリッド床 | synthesis（Cu/ZnO, VBF or KOGAS）+ dehydration（ZSM-5, Ortega）、`k_eq3="thermo"` |
| 体積分率 | synthesis 側 `F_SYN` = 1/1.9 |
| R2 カルボニル化床 | H-MOR、`DTU-Cheung2007-2`（K3 は逆数形 = reaction_rate cecc1da 以降） |
| 充填密度 | Cu/ZnO 1300 / ZSM-5 800 / H-MOR 837 kg/m³ |
| SV | 5000 /h（標準状態 22.414 L/mol 基準） |
| パージ | SP-gas で Purge : RecycleGas = 0.05 : 0.95 |
| 反応器 | KineticReactor: 出口 = `reaction_rate.pfr(入口)` を残差として課す |

成分: H2, CO, CO2, CH4, H2O, CH3OH, CH3OCH3, CH3COOCH3, CH3CHO, CH3COOH, N2

### 3.2 構成と分離の拘束（第 2 出口への回収率）

**plant3**（R1 → R2 直結、凝縮器 1 つ）

```
M1 → R1(合成+脱水) → R2(カルボニル化) → Condenser ┬ CondGas ─[T-101]→ SP-gas ┬ Purge
                                                  │                          └ RecycleGas → M1
                                                  └ CondLiq → Column1 ┬ MA 製品
                                                                      └ Sep1Liq → Column2 ┬ RecycleMeOH → M1
                                                                                          └ Water
```

| 分離器 → 第 2 出口 | 回収率 1 | 回収率 0 |
|---|---|---|
| Condenser → CondLiq | H2O, CH3OH, MA, CH3CHO, AcOH | H2, CO, CO2, CH4, N2, DME |
| Column1 → Sep1Liq | H2O, CH3OH, AcOH | MA, CH3CHO（→ MA 製品）、軽質ガス |
| Column2 → Water | H2O, AcOH | CH3OH（全量リサイクル）、その他 |

**plant4**（床の間で H2O / MeOH を抜く、凝縮器 2 つ）

```
M1 → R1 → Cond1 ┬ Cond1Gas ─[T-102]→ R2 → Cond2 ┬ Cond2Gas ─[T-101]→ SP-gas ┬ Purge
                │                               │                            └ RecycleGas → M1
                │                               └ Cond2Liq → Column3 ┬ MA 製品
                │                                                    └ Col3Bottoms（系外）
                └ Cond1Liq → Column1 ┬ Col1Top（実質 0）
                                     └ Sep1Liq → Column2 ┬ RecycleMeOH → M1
                                                         └ Water
```

| 分離器 → 第 2 出口 | 回収率 1 | 回収率 0 |
|---|---|---|
| Cond1 → Cond1Liq | H2O, CH3OH, MA, CH3CHO, AcOH | H2, CO, CO2, DME, CH4, N2（→ R2） |
| Column1 → Sep1Liq | H2O, CH3OH, AcOH | MA, CH3CHO, 軽質ガス |
| Column2 → Water | H2O, AcOH | CH3OH（全量リサイクル） |
| Cond2 → Cond2Liq | MA, H2O, CH3OH, CH3CHO, AcOH | **DME はガス側**に残して M1 へ、軽質ガス |
| Column3 → Col3Bottoms | MA 以外すべて | MA（製品として留出） |

CO2 除去塔: 塔入口ガスの CO2 回収率 = η（§0 の表の位置）。

### 3.3 触媒体積の外側反復

循環があるので反応器入口流量は解くまで分からない。そこで
「触媒体積を仮定 → フローシートを解く → 入口流量から体積を再計算」を繰り返す。

V_SV(n) = n [mol/h] × 22.414e-3 / 5000 [m³]

| | plant3 | plant4 `plant3` モード（既定） | plant4 `per-bed` モード |
|---|---|---|---|
| 未知 | V_tot（スカラー） | (V1, V2) | (V1, V2) |
| 目標 | V_tot = V_SV(ReactorIn) | V1 = V_SV(ReactorIn)/3、V2 = 2·V1 | V1 = V_SV(ReactorIn)、V2 = V_SV(R2 入口) |
| 配分 | hybrid : MA = 1 : 2 | 1 : 2（毎回 V2 = 2·V1 に射影） | 各床独立 |
| 初期値 | V_SV(6.5 × 新鮮供給) | V_SV(6.5 × 新鮮供給) を 1:2 | V1 = V_SV(6.5 × 供給)、V2 = V_SV(0.95 × 6.5 × 供給) |
| 更新 | 割線法（1 歩目は固定点ステップ） | 成分ごとの割線法 | 同左 |
| 1 回の更新幅 | 0.2〜5 倍 | 0.2〜5 倍 | 0.2〜5 倍 |
| 最大反復 | 8 | 10 | 10 |
| 収束判定 | \|g\| ≤ 1e-4 × V | 各成分 \|g\| ≤ 1e-4 × V | 同左 |

`6.5` は `RECYCLE_GUESS`（反応器入口 / 新鮮供給）、`0.95` は `GAS_FRACTION_GUESS`
（Cond1 のガス側に残る割合）。

> `example_plant3.py` の docstring には「外側の固定点反復・減衰0.5」とあるが、
> 現在のコードは割線法（`solve_with_sv`）。

### 3.4 内側の求解

| 項目 | 値 |
|---|---|
| 未知ストリームの初期推定 | 外側 1 回目は全成分 1.0 mol/h、2 回目以降は前回の解 |
| ソルバー | `least_squares`、流量の範囲 0〜∞ |
| 合格判定 | 残差ノルム < `--solve-tol` |
| 不合格時 | `--retry-tol` に緩めて 1 回だけ再試行、それでも駄目なら失敗として記録し次のケースへ |

---

## 4. 出力

`--outdir`（plant3 版は `examples/output/reformer/` 固定）にケースごとに書く。
`examples/output/` は `.gitignore` 対象なのでリポジトリには入らない。

| ファイル | 内容 |
|---|---|
| `co2removal_<label>_reformer.xlsx / .html` | 改質器のストリーム表 / Mermaid フロー図 |
| `co2removal_<label>_plant3.xlsx / .html`（plant4 版は `_plant4`） | 後工程のストリーム表 / Mermaid フロー図 |
| `co2removal_<label>_caustic.xlsx / .html` | 吸収塔サブフローシート（除去ありのケースのみ） |
| `co2removal_<label>_metrics.json` | 指標（成功時のみ。実行前に古いものを削除） |

ラベルは plant3 版が `R<改質>_<ケース>[_eta<η>]`、plant4 版が `p4R<改質>_<ケース>[_eta<η>]`。
