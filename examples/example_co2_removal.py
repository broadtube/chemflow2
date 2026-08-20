r"""サンプル: 循環系の CO2 希釈を解消する — 除去位置 (a)/(b) の比較。

example_pattern1_plant3.py の結論「不活性（CO2+CH4）最小化は MA 生産では裏目」を受けて、
**改質器への CO2 供給は A_base のまま（= ドライリフォーミングの取り分は捨てない）**、
CO2 を**どこで抜くか**だけを変えて比較する。

    baseline  除去なし（= A_base をそのまま plant3 へ）
    (a) 前段  改質器の下流・plant3 へ渡す手前で抜く   T-100（solve_reformer の中）
    (b) 循環  凝縮器出口ガスとパージ分岐の"間"で抜く  T-101（p3.build の中）

    (a)  RG+CO2+H2O → G1 → 凝縮 → DryGas ─[T-100]→ PlantFeed → plant3（塔なし）
    (b)  ... → DryGas → plant3: 反応器 → 凝縮器 → CondGas ─[T-101]→ ScrubbedGas
                                                                 → パージ / 循環

**(b) をパージより手前に置く理由**: CO2 を先に抜けばパージすべき蓄積成分が CH4 だけに
なるので、同じ CH4 濃度をパージ率を下げて達成できる（= CO パージ損失が減る）。
この効果は下の WGS の話とは独立に効く。本スクリプトは第一段としてパージ率を 5% に
固定して比べる（--purge で上書き可）。

────────────────────────────────────────────────────────────────────────
⚠️ 事前の予測: **CO2 除去は MA を減らす可能性が高い。**
────────────────────────────────────────────────────────────────────────
example_pattern1_plant3 の出力を物質収支で読み直すと、CO2 は全ケースで循環系の中で
**正味"生成"側**だった（供給 < パージ）。つまり動いているのは RWGS ではなく順方向の
水性ガスシフト CO + H2O → CO2 + H2 で、これは **CO を食う**反応。

    ケース      供給CO2   パージCO2   正味生成    MA
    A_base       1.675      1.767     +0.092    1.0416
    C_03MPaG     0.174      0.517     +0.344    0.9653

A→C で増えた正味生成 0.252 mol/h は、減った MA 0.076 × 3 = 0.229 mol の CO とほぼ一致
する。**ループ内の CO2 は WGS の生成物側阻害として CO を守っていた**、というのがデータの
読み。だとすれば CO2 を抜くほど WGS が進んで CO が失われる。

一方で希釈が解消して CO 分圧が上がる分、MA 合成速度は上がり触媒量は減る。差し引きは
計算しないと決まらない。**このスクリプトはその差し引きを測るためのもの**であって、
CO2 除去が良いという前提で書かれてはいない。

────────────────────────────────────────────────────────────────────────
✅ 検証結果（総当たり 28 runs、2026-08-20 改訂）
────────────────────────────────────────────────────────────────────────
予測（MA は減る）は当たった。**ただし「だから CO2 除去は駄目」という結論は誤り**
だった。2 つの前提を置き換えると評価が反転する。

**(1) 指標を投入炭素基準に取り直す。** 新鮮 CO 基準は分母が DryGas の CO だけで、
改質器に投入した CO2 と CH4 を勘定に入れていない。分母を投入炭素
（CH4 2.137 + 改質器 CO2 供給）にすると、改質条件の順位が完全に逆転する
（A_base 57.12% が最下位、D_01MPaG 75.06% が最良）。

**(2) 回収 CO2 を改質器へ戻す。** 改質器の CO2 供給は H2/CO=1.3837 の制約で
**総量が決まっている**ので、その一部を回収 CO2 で置き換えても Mixed 組成は同一
＝改質器・DryGas・下流の解はすべて不変。**新鮮 CO2 の勘定だけが変わる**ので、
再計算なしで正確に評価できる（紙計算ではない）。

    条件        ケース            MA      回収CO2  新鮮CO2  C効率(捨)  C効率(戻)  触媒
    D_01MPaG  b_recycle η0.9   0.9020    0.852    0.994    67.94%   **86.43%**  276mL
    D_01MPaG  b_recycle η0.5   0.9267    0.758    1.088    69.80%     86.20%   244mL
    C_03MPaG  b_recycle η0.9   0.8871    0.856    1.011    66.46%     84.53%   259mL
    A_base    b_recycle η0.99  0.8553    2.232    1.102    46.90%     79.21%   216mL
    D_01MPaG  baseline         0.9966    0       1.846    75.06%     75.06%   217mL
    A_base    baseline         1.0416    0       3.334    57.12%     57.12%   305mL

**回収 CO2 を捨てる前提なら CO2 除去は常に損（C効率が下がる）。戻す前提なら
常に得**（D_01MPaG で 75.06% → 86.43%）。MA は 9.5% 落ちるが、新鮮 CO2 の消費を
1.846 → 0.994 mol/h と 46% 削るので、炭素あたりでは大きく勝つ。

    ⚠ 代償: 触媒量 +27%（217 → 276 mL）。MA 減と合わせてトレードオフの評価が要る。
    ⚠ **NaOH では実現できない**。苛性ソーダ吸収は不可逆（Na2CO3 になって戻らない）
      なので回収 CO2 を改質器へ戻せない。アミン（MDEA）や PSA など再生可能な方式が
      前提になる。NaOH 吸収塔は原理確認用の実験装置であって、循環構成のモデルでは
      ない（references/naoh_absorber_sizing.html §7）。
    ⚠ 未評価: エネルギー（改質器熱負荷・循環コンプレッサ動力）、経済性。

**後工程を plant4 にすると成立しない。** plant4 は未反応 DME を循環させる設計なので
ループ CO2 を抜くと DME が濃縮してパージで失われ、MA が −57% 落ちる。回収 CO2 を
戻す前提でも最良は除去なし（63.15%）。詳細は example_co2_removal_plant4.py。

────────────────────────────────────────────────────────────────────────
苛性ソーダ吸収塔の扱い（4. の設計）
────────────────────────────────────────────────────────────────────────
吸収塔がガスループに及ぼす影響は「CO2 除去率 η」ただ 1 つに縮退する（苛性液側の情報は
ガスに戻らない）。そこで 2 つに分ける:

  主フローシート  Separator 1 個 + 回収率制約 1 本。循環求解（121変数・約10分）を
                  重くしない。上の T-100 / T-101 がこれ。
  吸収塔サブ      Mixer + Reactor + Separator。NaOH 消費・Na2CO3 廃液・原子収支を出す。
                  主解を入力に一瞬で解ける。caustic_scrubber() がこれ。

**新しい Unit は作らない**。Separator の設計思想（ユニットは収支だけを課し、分配は制約で
決める）にそのまま乗るうえ、Reactor を使えば Reaction の原子収支検査が自動でかかるため。

    CO2 + 2NaOH → Na2CO3 + H2O     C:1=1  O:4=4  H:2=2  Na:2=2 で閉じる

（液中では CO2 + OH⁻ → HCO3⁻ → CO3²⁻ の 2 段だが正味は同じ。過剰 NaOH がある限り
炭酸水素塩まで戻らないので 1 式で足りる。EXCESS がその余裕。）

**塔の寸法（塔高・充填量・塔径）は chemflow2 の外**。η → 寸法の手計算は
`references/naoh_absorber_sizing.html` にまとめてある（NTU = ln(1/(1-η))、
V = NTU·G/(k_G a·P)、Hatta 数による律速域の確認、実験規模の薬剤量と濡れ性の制約）。
結論だけ言うと、本件の規模では必要充填体積は数十 mL のオーダーで、**塔の寸法を実際に
決めるのは物質移動ではなく濡れ性・液分散・耐圧**のほう。

────────────────────────────────────────────────────────────────────────
実行方法
────────────────────────────────────────────────────────────────────────
要 Cantera（改質器）と reaction_rate（速度論反応器）。
PYTHONPATH の "." は examples.example_plant3 を絶対 import するため、
"../reaction_rate/src" は速度論 PFR のため。**どちらも必須**。

■ 総当たり（改質条件 4 種 × 除去位置 2 種 × η 3 点 = 28 runs, 約 5 時間）

    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal.py --reform A_base   --cases baseline,a_upstream,b_recycle --eta 0.5,0.9,0.99 --solve-tol 1e-4
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal.py --reform B_900C   --cases baseline,a_upstream,b_recycle --eta 0.5,0.9,0.99 --solve-tol 1e-4
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal.py --reform C_03MPaG --cases baseline,a_upstream,b_recycle --eta 0.5,0.9,0.99 --solve-tol 1e-4
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal.py --reform D_01MPaG --cases baseline,a_upstream,b_recycle --eta 0.5,0.9,0.99 --solve-tol 1e-4

    python3 examples/merge_xlsx.py --all     # 最後に 3 つの xlsx を 1 つに結合

Windows (cmd) はこう:

    set PYTHONPATH=.;..\reaction_rate\src
    python -u examples\example_co2_removal.py --reform C_03MPaG --cases baseline,a_upstream,b_recycle --eta 0.5,0.9,0.99 --solve-tol 1e-4

**`-u` と `--solve-tol 1e-4` を付ける理由:**

  -u             付けないと stdout がブロックバッファになり、リダイレクト時に
                 進捗が何も見えなくなる（実際に 52 分間無音になった事例あり）
  --solve-tol    既定 1e-6 は残差ノルムの**絶対値**なので、流量 50 mol/h 規模の
                 本問題では相対 2e-8 に相当し厳しすぎる。1e-4 なら相対 2e-6 で
                 精度は十分（plant4 版の既定と同じ）

⚠ **solve_tol は求解そのものを速くしない。** これは合否判定にしか使われず
（core/problem.py の `ok = resid_norm < tol`）、least_squares の停止条件は
solve_with_sv がハードコードしている ftol/xtol/gtol=1e-12 のほう。つまり:

    solve_tol を緩める効果 = 「判定を外して自動再試行に入る」のを防ぐ（時間が倍に
                             なるのを避ける）だけ。1 回目の求解時間は変わらない。

したがって 1 回目が長いケースでは、solve_tol を緩めても 1 回目は同じだけかかる。
効くのは「残差が 1e-6 に届かず 1e-4 には収まる」ケースで、そこで再試行を回避できる。
1 回目自体を速くしたいなら ftol/xtol/gtol を緩める必要があるが、それは解の精度に
直接効くので既定は変えていない。

**検算（A_base baseline）:** solve_tol を 1e-6 → 1e-4 にしても解は一致する。
所要時間も変わらない（上のとおり合否判定にしか使われないため）。

    指標              記録(1e-6)   今回(1e-4)
    酢酸メチル [mol/h]   1.0416     1.0416226491
    全触媒体積 [mL]     305.1059     305.1057
    反応器入口 [mol/h]   68.0614      68.0614
    所要時間               338 s        332 s

**改質条件ごとに 4 回に分けるのは仕様**。--reform は 1 条件しか取らない。
1 回の比較の中で改質条件を固定しないと「改質条件を変えた効果」と「除去位置の
効果」が混ざって読めなくなるため（BASE_CASE のコメント参照）。

■ 手早く 1 条件だけ（baseline + (a) + (b) を η=0.9 で = 3 runs, 約 30 分）

    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal.py --eta 0.9 --solve-tol 1e-4

■ 引数

    --reform     A_base / B_900C / C_03MPaG / D_01MPaG（既定 A_base）
                 example_pattern1_plant3.CASES のキーそのまま。T と圧力が名前に入る。
                 誤った名前は有効値を列挙して落ちるので、黙って変な結果にはならない。
    --cases      baseline / a_upstream / b_recycle をカンマ区切り
                 （既定は 3 つとも。baseline は η に依らないので 1 回だけ回る）
    --eta        CO2 除去率。カンマ区切りで複数指定すると η スイープ
    --purge      パージ率（既定 0.05）
    --progress-every  残差評価 N 回ごとに進捗を表示（既定 200、0 で無効）
    --solve-tol  Problem.solve の残差ノルム判定（既定 1e-6）
    --retry-tol  1 回目が収束判定を外したときの緩和 solve_tol（既定 1e-5）

⚠ solve_tol は残差ノルムの**絶対値**なのでスケール依存。流量 50 mol/h 規模の
本問題では既定の 1e-6 は相対 2e-8 に相当し、かなり厳しい（plant4 版は既定 1e-4 で、
「PFR の積分誤差で残差は 3e-6 で頭打ちになる」というコメント付き）。
**特定のケースだけ何十分も返ってこないときは、まず --solve-tol 1e-4 を試すこと。**
収束判定を外して自動再試行に入ると時間が倍になるが、最初から緩めておけば 1 回で済む。

■ 出力（examples/output/）

    co2removal_{改質条件}_{ケース}_eta{η}_reformer.xlsx / .html
    co2removal_{改質条件}_{ケース}_eta{η}_plant3.xlsx   / .html
    co2removal_{改質条件}_{ケース}_eta{η}_caustic.xlsx  / .html
    co2removal_{改質条件}_{ケース}_eta{η}_merged.xlsx   ← merge_xlsx.py が作る

xlsx は 3 つとも MERGE_BASIS / MERGE_COMPONENTS で行構造を揃えてあるので、
merge_xlsx.py が列を連結するだけで 1 枚のストリーム表になる。html は Mermaid 図。

■ 所要時間の目安（実測）

plant3 は 1 run 6〜18 分（収束の渋いケースで最長 50 分超）。28 runs で約 5 時間。
改質器と吸収塔サブフローシートは一瞬なので、時間はすべて plant3 の循環求解。
1 ケースが収束判定を外しても自動で 1 回だけ緩和再試行するので、途中で止まらない。

求解中は --progress-every の刻みでこう出る（無音にならない）:

      [外側1] nfev=   200  ‖residual‖=5.954e+00  経過 1.2 分
      [外側1] nfev=   400  ‖residual‖=2.117e+00  経過 2.5 分
      [1] V_cat=... → 反応器入口 ... → V_cat(SV=5000/h)=...   (64.4s, PFR積分 1228回, nfev=62)

‖residual‖ が階段状（しばらく一定 → 急に下がる）なのは正常。least_squares は
ヤコビアンを差分で作るので、変数 133 個なら約 133 回ごとに 1 歩進む。

⚠ Windows で `python -u` を付けずにリダイレクトすると stdout がブロックバッファに
なり、この進捗も含めて何も見えなくなる。**必ず -u を付けること。**
"""

import argparse
import os
import time

import examples.example_pattern1_plant3 as pp
import examples.example_plant3 as p3
from chemflow2 import (
    Mixer,
    Problem,
    Reaction,
    Reactor,
    Separator,
    Stream,
    StreamCondition,
    export_mermaid,
    stream_table,
    to_excel,
)

OUT = os.path.join(os.path.dirname(__file__), "output")

#: 比較の土台にする改質条件（既定は A_base = 現状条件）。--reform で切り替える。
#: 1 回の比較の中では固定すること。動かすと「改質条件を変えた効果」と
#: 「除去位置の効果」が混ざって読めなくなる。
#:
#: なお改質条件と CO2 除去は**同じレバー**（どちらもループ CO2 を下げる）なので、
#: 4 条件 × 2 位置 × η の総当たりは半分近くが縮退する:
#:   - (a) 前段除去は C/D では抜くものがほとんど無い（DryGas CO2 が 0.17/0.09 mol/h）
#:   - (b) 循環除去は 4 条件すべてで効く（ループ内で WGS が CO2 を作り続けるため）
#: 総当たりを回す前に、まず 1 条件で η を振って向きを決めること。
BASE_CASE = dict(pp.CASES["A_base"])

#: 苛性ソーダの化学量論に対する過剰率。過剰 NaOH が無いと炭酸水素塩側に寄るので余裕を持つ。
EXCESS = 1.2
#: 苛性ソーダ水溶液の濃度 [wt%]。工業的な希釈苛性の常用域。
NAOH_WT = 0.20

CAUSTIC = Reaction({"CO2": -1, "NaOH": -2, "Na2CO3": 1, "H2O": 1},
                   name="CO2 + 2NaOH -> Na2CO3 + H2O")


#: 3 つのフローシート（改質器 / 吸収塔 / plant3）の出力で**共通に使う成分行**。
#: 全フローシートの成分の和集合（13）を、軽質ガス → 含酸素 → 液/塩 の順に並べたもの。
#: これを to_excel(components=...) に渡すと **集合も順序も**固定されるので、
#: 3 ファイルの行構造が完全に一致し、merge_xlsx.py が列を連結するだけで済む。
#:
#: ⚠ Stream 側には持たせない。chemflow2 では流量ゼロの成分も変数 1 個 + 式 1 本を
#: 生むので、計算に寄与しないまま問題サイズだけ膨らむ（plant3 は 121 変数・約10分）。
#: あくまで**出力層だけ**の話。
MERGE_COMPONENTS = [
    "H2", "CO", "CO2", "CH4", "N2", "H2O",
    "CH3OH", "CH3OCH3", "CH3COOCH3", "CH3CHO", "CH3COOH",
    "NaOH", "Na2CO3",
]

#: 同じ理由で basis も 3 ファイル共通にする。to_excel は "mass" が bases に無いときだけ
#: 質量閉包行を足すので、揃えておけばその行の有無も自動で揃う。
#: 注: normal_volume は液ストリームでは物理的に無意味だが、既存の plant3 出力でも
#: 同じなので新たな問題ではない。
MERGE_BASIS = ["mol", "mole_frac", "mass", "normal_volume"]


def sheet_name(name: str) -> str:
    """Excel のシート名は 31 文字まで。超える分は頭から詰める。

    ラベルは「改質条件_ケース_eta_種別」の順で、種別（plant3/reformer/caustic）は
    ファイル名でも区別できるので、切るなら末尾から。openpyxl は 31 文字超でも
    書き込むが警告を出し、アプリによっては読めない。
    """
    return name[:31]


# --------------------------------------------------------------------------- #
# 吸収塔サブフローシート（既存ユニットの組み合わせ・新 Unit は作らない）
# --------------------------------------------------------------------------- #
def caustic_scrubber(gas_in: dict[str, float], eta: float, *,
                     excess: float = EXCESS, naoh_wt: float = NAOH_WT,
                     T: float = 40.0, P: str = "5MPaG", name: str = "T-x"):
    """苛性ソーダ吸収塔を Mixer + Reactor + Separator で組んで解く。

    gas_in  : 塔に入るガスのモル流量 [mol/h]（主フローシートの解をそのまま渡す）
    eta     : CO2 除去率（= Reactor の CO2 基準転化率そのもの）

    自由度: 変数 3·n_g + 7 / 方程式 3·(n_g+2) + 1 で閉じる（n_g = ガス成分数）。
    ガス側に残る水は「全量が液に落ちる」と近似した（51 bar・40 ℃ の飽和水分は
    y_H2O ≈ 0.15% で、NaOH 収支にも下流にも効かない）。

    戻り値 (problem, ScrubGas, Spent, Caustic)。
    """
    gas_comps = list(gas_in)
    abs_comps = gas_comps + [f for f in ("NaOH", "Na2CO3") if f not in gas_comps]
    liq_comps = ["H2O", "NaOH", "Na2CO3"]

    gas_c = StreamCondition(T=T, P=P, phase="gas")
    liq_c = StreamCondition(T=T, P=P, phase="liquid")

    # 苛性ソーダ供給は「除去する CO2 の 2 倍 × 過剰率」。ガスが既知なので前もって決まる。
    co2_removed = eta * gas_in.get("CO2", 0.0)
    naoh = 2.0 * excess * co2_removed
    # naoh_wt [wt%] の水溶液にするための水: m_w = m_NaOH·(1-wt)/wt
    water = naoh * 40.0 * (1.0 - naoh_wt) / naoh_wt / 18.015

    GasIn   = Stream(gas_comps, name="1. GasIn",   order=1, condition=gas_c, flows=gas_in)
    Caustic = Stream(["NaOH", "H2O"], name="2. Caustic", order=2, condition=liq_c,
                     flows={"NaOH": naoh, "H2O": water})
    AbsIn   = Stream(abs_comps, name="3. AbsIn",  order=3, condition=gas_c)
    AbsOut  = Stream(abs_comps, name="4. AbsOut", order=4, condition=gas_c)
    ScrubGas = Stream(gas_comps, name="5. ScrubbedGas", order=5, condition=gas_c)
    Spent    = Stream(liq_comps, name="6. SpentLiquor", order=6, condition=liq_c)

    problem = Problem(
        streams=[GasIn, Caustic, AbsIn, AbsOut, ScrubGas, Spent],
        units=[Mixer([GasIn, Caustic], AbsIn, name=f"{name} 気液接触"),
               Reactor(inlet=AbsIn, outlet=AbsOut, reactions=[CAUSTIC],
                       key_component="CO2", conversion=eta, name=f"{name} 化学吸収"),
               Separator(AbsOut, [ScrubGas, Spent], name=f"{name} 分液")],
        name=f"{name} NaOH scrubber (eta={eta})")
    problem.constrain_recovery(AbsOut, Spent, {"H2O": 1.0}, name="水は全量 液側へ")

    sol = problem.solve()
    if not sol.success:
        raise RuntimeError(f"吸収塔が収束せず: {sol}")
    return problem, ScrubGas, Spent, Caustic


def report_caustic(label: str, gas_in: dict[str, float], eta: float, *, name: str):
    """吸収塔を解いて薬剤収支を出し、xlsx と Mermaid を書き出す。"""
    problem, ScrubGas, Spent, Caustic = caustic_scrubber(gas_in, eta, name=name)
    co2_removed = eta * gas_in.get("CO2", 0.0)
    naoh = Caustic.flow_of("NaOH")
    print(f"\n--- {name} 苛性ソーダ吸収塔（{label}, η={eta:.0%}）---")
    print(f"塔入口 CO2        : {gas_in.get('CO2', 0.0):.4f} mol/h")
    print(f"除去 CO2          : {co2_removed:.4f} mol/h")
    print(f"NaOH 供給         : {naoh:.4f} mol/h = {naoh * 40.0:.1f} g/h"
          f"（量論の {EXCESS:.1f} 倍）")
    print(f"苛性液 供給       : {Caustic.flow_of('H2O') * 18.015 + naoh * 40.0:.0f} g/h"
          f"（{NAOH_WT:.0%} NaOH aq）")
    print(f"Na2CO3 生成       : {Spent.flow_of('Na2CO3'):.4f} mol/h = "
          f"{Spent.flow_of('Na2CO3') * 106.0:.1f} g/h")
    print(f"未反応 NaOH       : {Spent.flow_of('NaOH'):.4f} mol/h（0 未満なら過剰率不足）")
    print(f"廃液 計           : {float(Spent.total_flow.eval()):.4f} mol/h")
    print(f"塔出口ガス CO2    : {ScrubGas.flow_of('CO2'):.4f} mol/h")
    print(f"原子収支          : Reaction が生成時に検査済 "
          f"({CAUSTIC.element_balance()})")

    to_excel(problem.streams, os.path.join(OUT, f"co2removal_{label}_caustic.xlsx"),
             sheet=sheet_name(f"{label}_caustic"),
             basis=MERGE_BASIS, components=MERGE_COMPONENTS)
    export_mermaid(problem, os.path.join(OUT, f"co2removal_{label}_caustic.html"),
                   title=f"{name} NaOH scrubber ({label}, eta={eta:.0%})", style="diamond")
    return Spent


# --------------------------------------------------------------------------- #
# 3 ケース
# --------------------------------------------------------------------------- #
def run_case(tag: str, eta: float, purge: float, reform: dict, label: str,
             solve_tol: float = 1e-6, progress_every: int = 0) -> dict:
    """1 ケース解いて、改質側 + plant3 側の指標を返す。

    tag は "baseline" / "a_upstream" / "b_recycle"、label は出力ファイル名の接頭辞。
    solve_tol は plant3 の残差ノルム判定（p3.solve_with_sv 参照）。
    """
    print(f"\n{'=' * 78}\n=== {label}（η={eta:.0%}, パージ率={purge:.0%}）\n{'=' * 78}")

    # --- 前工程（改質器）。(a) のときだけ改質器側に T-100 を置く ---
    ref_removal = eta if tag == "a_upstream" else None
    rp, DryGas, CO2f, H2Of, PlantFeed = pp.solve_reformer(
        verbose=False, co2_removal=ref_removal, **reform)
    print(f"改質器: H2O={H2Of.flow_of('H2O'):.4f} CO2供給={CO2f.flow_of('CO2'):.4f} mol/h")
    print(stream_table(rp.streams, basis=["mol", "mole_frac"]))
    to_excel(rp.streams, os.path.join(OUT, f"co2removal_{label}_reformer.xlsx"),
             sheet=sheet_name(f"{label}_reformer"),
             basis=MERGE_BASIS, components=MERGE_COMPONENTS)
    export_mermaid(rp, os.path.join(OUT, f"co2removal_{label}_reformer.html"),
                   title=f"Reformer + CO2 removal ({label})", style="diamond")

    # (a) の吸収塔は改質器出口ガスを処理する
    if tag == "a_upstream":
        report_caustic(label, {f: DryGas.flow_of(f) for f in DryGas.formulas}, eta,
                       name="T-100")

    # --- 後工程（plant3）。(b) のときだけ plant3 の中に T-101 を置く ---
    feed = {f: PlantFeed.flow_of(f) for f in PlantFeed.formulas}
    plant_removal = eta if tag == "b_recycle" else None
    problem, streams, v_tot = p3.solve_with_sv(feed=feed, co2_removal=plant_removal,
                                               purge=purge, solve_tol=solve_tol,
                                               progress_every=progress_every)
    by_name = {s.name: s for s in problem.streams}
    Steam1, ReactorIn = by_name["1. Steam1"], by_name["2. ReactorIn"]
    Purge, MA = by_name["7. Purge"], by_name["9. MethylAcetate"]
    print(stream_table(problem.streams, basis=["mol", "mole_frac"]))
    to_excel(problem.streams, os.path.join(OUT, f"co2removal_{label}_plant3.xlsx"),
             sheet=sheet_name(f"{label}_plant3"),
             basis=MERGE_BASIS, components=MERGE_COMPONENTS)
    export_mermaid(problem, os.path.join(OUT, f"co2removal_{label}_plant3.html"),
                   title=f"plant3 + CO2 removal ({label})", style="diamond")

    # (b) の吸収塔は凝縮器出口ガスを処理する
    if tag == "b_recycle":
        cond_gas = by_name["5. CondGas"]
        report_caustic(label, {f: cond_gas.flow_of(f) for f in cond_gas.formulas}, eta,
                       name="T-101")

    n_in = float(ReactorIn.total_flow.eval())
    # 抜いた CO2 は (a) なら改質器側の、(b) なら plant3 側のベントに出る
    ref_vent = next((s for s in rp.streams if s.name == "9. CO2Vent"), None)
    up_removed = ref_vent.flow_of("CO2") if ref_vent is not None else 0.0
    loop_removed = (by_name["14. CO2Vent"].flow_of("CO2")
                    if "14. CO2Vent" in by_name else 0.0)
    removed = up_removed + loop_removed
    return {
        "供給 CO2 [mol/h]": Steam1.flow_of("CO2"),
        "除去 CO2 [mol/h]": removed,
        "反応器入口 [mol/h]": n_in,
        "反応器入口 CO2 [mol%]": ReactorIn.flow_of("CO2") / n_in * 100,
        "反応器入口 CO  [mol%]": ReactorIn.flow_of("CO") / n_in * 100,
        "反応器入口 CH4 [mol%]": ReactorIn.flow_of("CH4") / n_in * 100,
        "全触媒体積 [mL]": v_tot * 1e6,
        "酢酸メチル [mol/h]": MA.flow_of("CH3COOCH3"),
        "CO パージ損失 [mol/h]": Purge.flow_of("CO"),
        "CO2 パージ [mol/h]": Purge.flow_of("CO2"),
        "新鮮CO基準 炭素収率 [%]": MA.flow_of("CH3COOCH3") * 3 / Steam1.flow_of("CO") * 100,
        # plant3 ループの CO2 収支（出 − 入）。正 = ループ内で CO が CO2 に食われている量
        # （順方向 WGS の正味進行度）。(a) の前段除去は plant3 の外なので入れない。
        "CO2 正味生成 [mol/h]": (Purge.flow_of("CO2") + loop_removed
                                 - Steam1.flow_of("CO2")),
    }


def main():
    ap = argparse.ArgumentParser(description="CO2 除去位置 (a)/(b) の比較")
    ap.add_argument("--eta", default="0.90",
                    help="CO2 除去率。カンマ区切りで複数指定すると η スイープになる"
                         "（例 --eta 0.5,0.9,0.99）")
    ap.add_argument("--cases", default="baseline,a_upstream,b_recycle",
                    help="実行するケースをカンマ区切りで（baseline/a_upstream/b_recycle）")
    ap.add_argument("--reform", default="A_base",
                    help=f"土台にする改質条件（{'/'.join(pp.CASES)}）")
    ap.add_argument("--purge", type=float, default=p3.PURGE,
                    help=f"パージ率（既定 {p3.PURGE}）")
    ap.add_argument("--progress-every", type=int, default=200,
                    help="残差評価 N 回ごとに進捗を表示（既定 200、0 で無効）。"
                         "1 回の求解が数十分かかることがあり、無音だと計算中か"
                         "固まったか判別できないため")
    ap.add_argument("--solve-tol", type=float, default=1e-6,
                    help="Problem.solve の残差ノルム判定（既定 1e-6）。**絶対値**なので"
                         "スケール依存で、流量 50 mol/h 規模の本問題では 1e-6 は相対 2e-8 に"
                         "相当し厳しい。収束に時間がかかるケースは 1e-5〜1e-4 で回すとよい"
                         "（plant4 版は既定 1e-4）")
    ap.add_argument("--retry-tol", type=float, default=1e-5,
                    help="1 回目が収束判定を外したときに使う緩和 solve_tol（既定 1e-5）")
    args = ap.parse_args()

    if args.reform not in pp.CASES:
        ap.error(f"--reform は {list(pp.CASES)} のいずれか")
    reform = dict(pp.CASES[args.reform])
    etas = [float(e) for e in args.eta.split(",") if e.strip()]

    os.makedirs(OUT, exist_ok=True)
    tags = [t.strip() for t in args.cases.split(",") if t.strip()]

    # baseline は η に依らないので 1 回だけ。(a)/(b) は η ごとに回す。
    plan: list[tuple[str, float, str]] = []
    for tag in tags:
        if tag == "baseline":
            plan.append((tag, 0.0, f"{args.reform}_baseline"))
        else:
            for eta in etas:
                plan.append((tag, eta, f"{args.reform}_{tag}_eta{eta:g}"))

    print(f"実行計画: {len(plan)} runs（plant3 は 1 run 約 10 分 → 約 "
          f"{len(plan) * 10} 分, solve_tol={args.solve_tol:g}）")
    for _tag, eta, label in plan:
        print(f"  - {label}  η={eta:.0%}")

    # 総当たりは数時間走るので、1 ケースが収束しなくても残りを続ける。
    # 落ちたケースは results に入れず、最後にまとめて報告する。
    #
    # 1 度目に失敗したら判定を 1 桁緩めて 1 回だけ再試行する。solve_tol は残差ノルムの
    # **絶対値**なのでスケール依存で、既定の 1e-6 は流量 50 mol/h 規模の本問題では
    # 相対 2e-8 に相当し厳しすぎる。実際 1.8e-6 で頭打ちになって落ちる点があった
    # （物理的には収束済み）。緩めた解は results に retried 印をつけて区別する。
    results, failed, retried = {}, [], []
    for i, (tag, eta, label) in enumerate(plan, 1):
        t0 = time.time()
        try:
            results[label] = run_case(tag, eta, args.purge, reform, label,
                                      solve_tol=args.solve_tol,
                                      progress_every=args.progress_every)
            print(f"[{i}/{len(plan)}] {label} 完了 ({time.time() - t0:.0f}s)")
        except Exception as e:                      # noqa: BLE001 — 打ち切らず続行
            print(f"[{i}/{len(plan)}] {label} 1回目失敗: {type(e).__name__}: {e}")
            print(f"    → solve_tol を {args.retry_tol:g} に緩めて再試行")
            try:
                results[label] = run_case(tag, eta, args.purge, reform, label,
                                          solve_tol=args.retry_tol,
                                          progress_every=args.progress_every)
                retried.append(label)
                print(f"[{i}/{len(plan)}] {label} 完了（緩和判定） "
                      f"({time.time() - t0:.0f}s)")
            except Exception as e2:                 # noqa: BLE001
                failed.append((label, f"{type(e2).__name__}: {e2}"))
                print(f"[{i}/{len(plan)}] {label} 失敗 ({time.time() - t0:.0f}s): "
                      f"{type(e2).__name__}: {e2}")

    if len(results) > 1:
        print(f"\n{'=' * 78}\n=== 除去位置の比較（改質条件 {args.reform} 固定・"
              f"パージ率 {args.purge:.0%}）\n{'=' * 78}")
        metrics = list(next(iter(results.values())))
        # 列見出しは改質条件名（全列共通）を落として短くする
        heads = {t: t[len(args.reform) + 1:] or "baseline" for t in results}
        print(f"{'指標':>24s}" + "".join(f"{heads[t]:>18s}" for t in results))
        for m in metrics:
            print(f"{m:>24s}" + "".join(f"{results[t][m]:18.4f}" for t in results))
        print("\n読み方: 「CO2 正味生成」が baseline より増えていれば、CO2 を抜いたぶん\n"
              "順方向 WGS が進んで CO が食われている。それが MA の減少と 1:3 で対応する\n"
              "なら、docstring 冒頭の予測どおり。対応しないなら別の機構を探す必要がある。")

    if retried:
        print(f"\n※ 緩和判定（solve_tol={args.retry_tol:g}）で解いたケース {len(retried)} 件: "
              f"{', '.join(retried)}")
    if failed:
        print(f"\n⚠ 収束しなかったケース {len(failed)} 件:")
        for label, msg in failed:
            print(f"  - {label}: {msg}")


if __name__ == "__main__":
    main()
