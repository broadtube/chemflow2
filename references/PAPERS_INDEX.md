# 文献インデックス（CO2 改質・炭素析出・炭素活量）

調査日: 2026-08-01 / 更新: 2026-08-02
フォルダ名は sibling プロジェクト `../reaction_rate/references/` の慣例に合わせた。

炭素活量 a_C ≤ 1 を pattern1（改質器）の運転制約として課すための根拠文献。

**→ 内容のまとめは [`carbon_activity.html`](carbon_activity.html)**（定義・反応式・平衡定数・
析出限界線図・金属比較・原論文との照合・実装。`../reaction_rate/references/rate_equations.html` と同形式）

> ⚠️ **論文 PDF はこのリポジトリに含まれない**（`.gitignore` で `references/*.pdf` を除外）。
> 本リポジトリは公開で、各 PDF の著作権は出版社にあるため。**「✅取得済」は調査時に
> 手元で確認したという意味**で、リポジトリ内にあるという意味ではない。書誌・DOI・入手先は
> すべて記載してあるので、必要なら各自で取得すること。

---

## 0. 炭素活量 — 一次出典と主根拠 ★2026-08-02 更新

### 0-1. ★取得済 — §5（金属比較・ウィスカーカーボン）の主根拠

- **J. R. Rostrup-Nielsen, J.-H. Bak Hansen. "CO2-Reforming of Methane over Transition
  Metals." *J. Catal.* 144(1), 38–49 (1993). DOI 10.1006/jcat.1993.1312**
- → `rostrupnielsen1993.pdf`（全 12 頁・**スキャン、テキスト層なし**）
- abstract: "All catalysts show **smaller equilibrium constants for methane decomposition
  than that based on graphite, the effect being largest for the noble metals**.
  Ru and Rh show high selectivity for carbon-free operation…"
- **Table 2**: 改質活性 TOF。序列 **Ru, Rh > Ir > Ni, Pt, Pd**
- **Table 3**: 炭素形態。**ウィスカーを作るのは Ni と Pd のみ**（Rh/Pt のものは Fe 汚染由来）
- **Fig. 4**: CH4 分解平衡定数。Graphite > Ni-a(300nm) > Ni-b(<10nm) > Ni-1 > Rh > Ru > Ir > Pt > Pd
  → 粒子径が小さいほど、また貴金属ほど析出しにくい。**a_C ≤ 1 は Ni で 1.5〜4 倍、
  貴金属で 5〜7 倍の保守側**（値は図の目視読み取り・±30% 程度）
- AIChE 2008 の参考文献 1

### 0-2. ★取得済 — 千代田自身の定義式（唯一の一次資料）

- **坂口 順一（千代田化工建設）「ＣＯ2を排出しない合成ガスの製造方法」WO2012140994A1**
  優先日 2011-04-12 / 公開 2012-10-18 / 32 pp.
- → `WO2012140994A1_Chiyoda_syngas.pdf`（**スキャン、テキスト層なし**）
- 「**カーボン活性＝Ｋ×（Ｐco)²／（Ｐco2)**」「この**カーボン活性の値が１を超えると**
  触媒上にカーボンが析出しやすくなる。」＝ **Boudouard 経路そのもの**。
  AIChE 2008 は式を書いていないので、**千代田の定義が確認できるのはこの特許だけ**。
- ⚠️ 引用は Google Patents のテキスト層から。保存 PDF からは直接読めない。

### 0-3. 概念の原著（いずれも未入手・DOI 解決確認済）

- J. R. Rostrup-Nielsen. "Equilibria of decomposition reactions of carbon monoxide
  and methane over nickel catalysts." *J. Catal.* 27(3), 343–356 (1972).
  DOI 10.1016/0021-9517(72)90170-4
  〔**概念の原著**。黒鉛平衡からのずれを最初に実測。1993 年論文の土台〕
- J. Rostrup-Nielsen, D. L. Trimm. "Mechanisms of carbon formation on nickel-containing
  catalysts." *J. Catal.* 48(1), 155–165 (1977). DOI 10.1016/0021-9517(77)90087-2
  〔※Crossref は著者 1 名しか返さないため共著者 Trimm は書誌検索ベース〕
- J. R. Rostrup-Nielsen. "Sulfur-passivated nickel catalysts for carbon-free steam
  reforming of methane." *J. Catal.* 85(1), 31–43 (1984). DOI 10.1016/0021-9517(84)90107-6
  〔SPARG の基礎〕
- I. Alstrup. "A new model explaining carbon filament growth on nickel, iron, and Ni–Cu
  alloy catalysts." *J. Catal.* 109(2), 241–251 (1988). DOI 10.1016/0021-9517(88)90207-2
  〔**粒子径補正 ΔG_c,dev = 2.6 + 93/d(nm) の一次候補。当該式がこの論文にあるかは未確認**。
  式自体は二次資料（*Appl. Catal. A* 2018 の whisker carbon 論文、有料）経由。
  **1993 Fig. 4 が実測で同じ物理を与えるので、そちらを主根拠にした**〕

### 0-4. 二次資料（無料・要注意）

- **R. A. Dagle 他. "Review of Novel Catalysts for Biomass Tar Cracking and Methane
  Reforming." PNNL-16950, Pacific Northwest National Laboratory (2007).**
- → `PNNL-16950_tar_reforming_catalyst_review.pdf`（テキスト層あり）
- **価値**: 米国政府技術レポート＝**無料で誰でも入手でき、引用も再配布も自由**。
- **性質**: バイオマスガス化のタール改質**触媒**レビュー。目次は
  1 序論 / 2 貴金属触媒によるメタン改質 / 3 代替 Ni 系触媒 / 4 代替触媒によるタール分解 / 5 結論。
- **どこまで**: 1993 年論文は **§2 の数段落で実験条件と結果を要約するのみ**。
  TOF 比（Ni 1.4–1.6 / Rh 4.3 / Ru 3.1 / Ir 10.2 / Pd 8.9 / Pt 5.6）は原著 Table 2 と一致。
- **⚠️ 誤り**: 活性序列を「Rh, Ru > Ir > Ni, Pt, Pd」と書くが、**原著は「Ru, Rh > …」で
  Ru と Rh が逆**（原著 Table 2 は Ru 8.9 > Rh 8.1）。
- **⚠️ 扱っていない**: 全文検索で `carbon activity` 0 件 / `carbon limit` 0 件 /
  `thermodynamic` 0 件 / `Alstrup` 0 件。**炭素活量・熱力学的判定には一切踏み込まない**。
- **使いどころ**: 原著が入手できない人向けの代替。序列の誤りに注意して補助的に使う。

---

## 1. 千代田化工建設 CT-CO2AR 系列

### 1-1. ★取得済み — 炭素活量の図が載っている一次資料

- Mikuriya, T.; Yagi, F.; Shimura, M.（千代田化工建設）
  "**Development of a new Synthesis Gas Production Catalyst and Process with
  CO2 and H2O Reforming**"
  *AIChE Annual Meeting 2008*, Philadelphia, PA, Nov. 16–21, 2008. 全8頁.
- <https://skoge.folk.ntnu.no/prost/proceedings/aiche-2008/data/papers/P140551.pdf>
- 入手状況: **取得済み** → `Mikuriya_Yagi_Shimura_2008_AIChE_CO2_H2O_reforming.pdf`（317 KB）

**確認した内容（全文検索済み）**

- "Carbon Activity" の語は本文に **1 回だけ**登場し、それは
  **Fig.4 "Feed Gas Conditions for Process Studies" 中の等高線 "Carbon Activity=1"**。
- Fig.4 の軸: 縦 = **H/C 比**、横 = **O/C 比**（いずれも原子比）。
  等温線 600 / 700 / 800 / 900 / 950 ℃。平衡組成は
  **リフォーマー出口 950 ℃・1.5 MPa** で計算したもの。
- ⚠️ **本論文は a_C の定義式を書いていない。** 図は "carbon limit diagram⁴⁾" として
  **参考文献 4** を引いている（→ §1-3）。
- 本文の要点（原文引用）:
  > operations on the **left side of the carbon limit curve are more economical**
  > conditions because lower steam and carbon dioxide addition in feed to produce
  > the synthesis gas for a given H2/CO ratio. But in these conditions, product
  > synthesis gases have **high thermodynamic potential for carbon formation**.
  > In our process, feed gas conditions on the left side of carbon limit line
  > **could be selected because our catalyst has high resistance against carbon
  > formation**.

**⚠️ 設計方針への含意**

千代田の技術は「a_C ≤ 1 を守る技術」ではなく、**a_C > 1 の領域でも触媒の耐炭素析出性で
運転できる**ことが売り。したがって本プロジェクトが課す **a_C ≤ 1 は、千代田触媒を
前提にしない保守的な条件**にあたる。

- パイロット実績: H2/CO = 2.0 で約 7,000 時間、H2/CO = 1.0 で約 1,000 時間の安定運転。
- 触媒設計: 貴金属をエッグシェル状に外表面へ選択担持。強ルイス塩基性の金属酸化物担体。

### 1-2. ★未取得 — 同じ内容の日本語版（著者 八木冬樹）

いずれも千代田化工建設・八木冬樹による解説。**3 件とも J-STAGE には本文が無く、
オープンアクセスではない**（日本エネルギー学会誌 91(4) の J-STAGE 目次は p.303 から
始まっており、解説欄の p.274 は収録されていない）。

| # | 誌名 | 巻号・頁 | 年 | タイトル |
|---|---|---|---|---|
| J1 | **触媒** | 51(5), 331–335 | 2009 | CO2 を原料として合成ガスを製造する CO2/H2O リフォーミングプロセス |
| J2 | **ペトロテック** | 34(6), 378–383 | 2011 | CO2 を原料とする合成ガス製造プロセス |
| J3 | **日本エネルギー学会誌** | 91(4), 274–278 | 2012 | CO2 リフォーミングプロセスによる合成ガス製造 |

- J3 の CiNii CRID: `1520290882653584768`
  NDL: <https://ndlsearch.ndl.go.jp/books/R000000004-I023744959>（参考文献 10 件収録）
- 入手経路の候補: 国立国会図書館の遠隔複写、大学図書館の相互貸借、各学会の会員向け
  アーカイブ（触媒学会・石油学会・日本エネルギー学会）。
- 優先度: **J1（触媒 2009）> J3（エネルギー学会誌 2012）> J2（ペトロテック 2011）**。
  J1 が AIChE 2008 の直後で内容が最も近いと見込まれる。

### 1-3. ★未取得・最重要 — 炭素活量の定義の一次出典

AIChE 2008 論文の参考文献 4。carbon limit diagram の出典。

- Udengaard, N. R.; Bak Hansen, J.-H.; Hanson, D. C.; Stal, J. A.（Haldor Topsøe）
  "**Sulfur passivated reforming process lowers syngas H2/CO ratio**"
  *Oil & Gas Journal* **1992**, 90(10), 62–67（Mar. 9, 1992）.
- OSTI 書誌: <https://www.osti.gov/biblio/5650478>
- 内容: Topsøe の **SPARG プロセス**（硫黄不動態化改質。1987 年に Sterling 社
  Texas City プラントで商業化）。硫黄を活性点に選択吸着させて炭素核生成を阻害する。
- **a_C の定義式を確認するならこの文献。** ただし商業誌（Oil & Gas Journal）なので
  J-STAGE 等では入手できない。OSTI 経由か図書館経由。

---

## 2. 炭素活量 a_C の定義（現行実装の根拠）

**§0-2 の特許 WO2012140994A1 で千代田自身の定義（Boudouard 経路）が確認できた**ので、
実装の根拠は「教科書的定義」から「一次資料で確認済み」に格上げされている。
§1-3 の Udengaard（線図としての一次出典）は依然未入手だが、優先度は下がった。

炭素析出を起こしうる 3 反応それぞれについて、固体炭素（黒鉛）の活量を求める。
千代田は B のみで定義しているが、**平衡気相では 3 経路が一致する**ので
改質器出口に適用する限り同値（→ `carbon_activity.html` §2）。

| # | 反応 | a_C |
|---|---|---|
| B | Boudouard: 2CO ⇌ C(gr) + CO2 | a_C = K_B · p_CO² / p_CO2 |
| M | メタン分解: CH4 ⇌ C(gr) + 2H2 | a_C = K_M · p_CH4 / p_H2² |
| R | CO 還元: CO + H2 ⇌ C(gr) + H2O | a_C = K_R · p_CO·p_H2 / p_H2O |

- **a_C > 1 なら炭素析出が熱力学的に有利**、a_C ≤ 1 なら析出しない。
- 分圧は **atm 基準**（Cantera の `reference_pressure` が 101325 Pa のため）。
- 3 つのうち**最大値**で判定する（どれか 1 つでも 1 を超えれば析出しうる）。

平衡定数は Cantera の標準生成ギブズエネルギーから計算する:
気相は `gri30.yaml`、固体炭素は `graphite.yaml`（相 `C(gr)`）。

### 検算（2026-08-01 実施）

| T [K] | K_B (Boudouard) | K_M (CH4分解) | K_R (CO還元) |
|---|---|---|---|
| 773.15 (500℃) | 2.4296e+02 | 4.7347e-01 | 4.7471e+01 |
| 973.15 (700℃) | 9.9978e-01 | 7.7947e+00 | 6.2037e-01 |
| 1123.15 (850℃) | 6.0235e-02 | 3.4111e+01 | 6.5862e-02 |
| 1173.15 (900℃) | 2.7833e-02 | 5.1449e+01 | 3.5432e-02 |
| 1273.15 (1000℃) | 7.1738e-03 | 1.0646e+02 | 1.1877e-02 |

妥当性の確認:
- Boudouard は**発熱**反応なので高温で K が下がる（500℃ 2.43e+2 → 1000℃ 7.17e-3）✓
- メタン分解は**吸熱**なので高温で K が上がる（500℃ 4.73e-1 → 1000℃ 1.06e+2）✓
- **K_B ≈ 1.000 が 973.15 K（700℃）** で得られる。これは Boudouard 平衡の
  古典的な基準温度と一致する ✓

---

## 3. 炭素析出の平衡計算を扱った日本語論文（検証用・オープンアクセス）

千代田系ではないが、a_C ≤ 1 制約の妥当性を実験値で確認するのに使える。

### ★★★ 未取得・最優先
- 村長 潔（東北肥料）「**メタン-水蒸気-炭酸ガス反応における炭素析出とモル比**」
  *工業化学雑誌* **1961**, 64(11), 1933–1938.
- DOI: 10.1246/nikkashi1898.64.11_1933 ／ **無料 PDF（1577 KB）**
- <https://www.jstage.jst.go.jp/article/nikkashi1898/64/11/64_11_1933/_article/-char/ja/>
- なぜ要るか: 595–982 ℃ で炭素析出曲線の**実験値と平衡計算値を直接比較**している。
  抄録より「温度816℃以上では炭素析出曲線の実験値と理論値はほぼ一致するが,706℃以下では
  全般に実験値の方が低く」→ **816 ℃以上なら平衡計算で析出境界を予測してよい**根拠。
  本プロジェクトの運転温度 850–900 ℃ はこの範囲に入る。

### ★★☆ 未取得
- 冨重 圭一, 藤元 薫（東京大学）「メタンの炭酸ガスリフォーミング反応による合成ガス
  製造における触媒特性と炭素析出挙動解析 — NiO-MgO 固溶体触媒の開発」
  *石油学会誌* **2001**, 44(2), 65–79. DOI: 10.1627/jpi1958.44.65 ／ **無料 PDF（3612 KB）**
- <https://www.jstage.jst.go.jp/article/jpi1958/44/2/44_2_65/_article/-char/ja/>
- 常圧では耐炭素析出だが**高圧では対策が必要**と報告。本件は 1.04 MPaG 運転なので関連。

### ✗ 取得不要と判断
- 森田 義郎「石炭ガス化反応の基礎」*日本エネルギー学会誌* **1978**, 58(2), 141–.
  → PDF を取得して全文検索したが「炭素活量」「活量」「Boudouard」「炭素析出」
    いずれも **0 件**。ガス化速度論の総説で用途が合わない。
