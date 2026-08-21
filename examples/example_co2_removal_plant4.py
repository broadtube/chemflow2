r"""サンプル: CO2 除去位置の比較 — 後工程を **plant4（2段反応器）** にした版。

example_co2_removal.py（後工程 plant3）の plant4 版。共通部分（苛性ソーダ吸収塔の
サブフローシート・成分行・basis）は例の相互 import で再利用し、**plant4 固有の部分
だけ**を書いてある。plant3 版とは別ファイルにしてあるのは:

  - plant4 には plant3 に無い除去位置 (c) 段間 がある → --cases の選択肢自体が違う
  - 触媒体積が plant3 はスカラー、plant4 は (V1, V2) の 2 つ → 指標の構造が違う
  - 仮説が違う（下記）ので docstring が別物になる

リポジトリの流儀（example_plant → plant2 → plant3 → plant4 がコピー＆改変の系列）
にも合わせている。3 つ目のプラントが出てきたら共通化を検討する。

────────────────────────────────────────────────────────────────────────
plant4 の構成と、除去位置が 3 つある理由
────────────────────────────────────────────────────────────────────────
    plant3:  M1 → R1(合成+脱水) → R2(カルボニル化) → Condenser → SP-gas
    plant4:  M1 → R1(合成+脱水) → Cond1 ┬ ガス → R2 → Cond2 ┬ ガス → SP-gas ┬ Purge
                                        └ 液(H2O,MeOH)      └ 液 → Col3 → MA

plant4 は **R1 と R2 の間で H2O と MeOH を抜く**（H-MOR のカルボニル化は水で強く
阻害されるため）。この Cond1 の存在が、除去位置と機構の両方に効いてくる。

    (a) 前段   改質器の下流・plant4 へ渡す手前        T-100  ← plant3 版と同じ
    (b) 循環   Cond2Gas とパージ分岐の間              T-101  ← plant3 の (b) に相当
    (c) 段間   Cond1Gas と R2 の間                    T-102  ← **plant4 でのみ可能**

(c) は Cond1 が既に水を抜いている場所なので、そこで CO2 も抜けば**カルボニル化床に
入る CO 分圧を直接上げられる**。H-MOR は水で阻害されるうえ CO2 で希釈されると不利
なので、plant3 には作れなかったこの位置に見込みがある。

────────────────────────────────────────────────────────────────────────
仮説: plant4 なら CO2 除去が効くかもしれない
────────────────────────────────────────────────────────────────────────
plant3 では **CO2 除去は 28 runs すべてで酢酸メチルを減らした**（例外なし）。
新鮮 CO の行き先を数えると、除去なしと (b) η=0.9 の差は完全に説明できた:

    酢酸メチル（3 CO/MA）   3.1248 → 2.5575   −0.5673
    順方向 WGS で CO2 へ    0.0922 → 0.5484   +0.4562  ← 損失の 80%
    DME としてパージ        0.0652 → 0.3196   +0.2544  ← 損失の 45%
    CO のままパージ         0.2675 → 0.1245   −0.1430  ← 25% 戻る

主因は **順方向 WGS（CO + H2O → CO2 + H2）が CO を食う**こと。ループ CO2 を抜くほど
生成物側阻害が外れて WGS が進む。

**plant4 は Cond1 で H2O を抜くので、この機構が弱まるはず**:
  - WGS が進める場所が R1 に限られ、R2 は乾いたガスで運転される
  - CO を食う経路が短くなる
→ plant4 は CO2 除去が実際に効くかもしれない唯一の構成。

────────────────────────────────────────────────────────────────────────
❌ 検証結果（A_base・η=0.9・3 runs）: **仮説は外れた。plant4 のほうがはるかに悪い。**
────────────────────────────────────────────────────────────────────────

    指標                    除去なし   (b)循環    (c)段間
    酢酸メチル [mol/h]       1.0651    0.4584    0.4582   ← **−57%**
    新鮮CO基準 炭素収率 [%]   90.01     38.74     38.72
    反応器入口 [mol/h]       62.95     47.38     47.38
    反応器入口 CO2 [mol%]    55.94      4.06      4.06   ← 狙いどおり下がる
    反応器入口 CO  [mol%]    10.39      8.97      8.97   ← **なのに CO も下がる**
    全触媒体積 [mL]         537.99    398.44    387.90
    CO パージ損失 [mol/h]    0.1573    0.0369    0.0369
    DME パージ [mol/h]       0.0536    0.7196    0.7198   ← **13 倍**
    CO2 正味生成 [mol/h]     0.0901    0.6987    0.6988

**WGS は確かに弱まったが、別の漏れ口が開いた。** 新鮮 CO の行き先で差を数えると
（供給 3.5498 に対し 3 ケースとも合計 3.5498 で厳密に閉じる）:

    酢酸メチル（3 CO/MA）  −1.8201
      順方向 WGS で CO2 へ  +0.6086   損失の 33%（plant3 では 80%）
      DME としてパージ      +1.3320   損失の 73%（plant3 では 45%）← **主因が入れ替わった**
      CO のままパージ       −0.1204   7% 戻る

**なぜ plant4 のほうが DME に弱いか:** plant4 は Cond2 で**未反応 DME をあえてガス側に
残して M1 へ戻す**設計（docstring 冒頭の構成図）。CO2 を抜くとループ全体の流量が
62.95 → 47.38 mol/h に縮み、CO の絶対量も 6.54 → 4.25 mol/h に落ちる。すると
カルボニル化が追いつかず DME が溜まり、パージ率 5% がその濃い DME を大量に系外へ
持ち出す。**CO2 を抜いたのに反応器入口の CO 濃度まで下がっている**のがその現れ。

────────────────────────────────────────────────────────────────────────
総当たり 40 runs の結果（4 改質条件 × 3 位置 × η 3 点、9.4 時間、40/40 成功）
────────────────────────────────────────────────────────────────────────

**(1) (b) 循環と (c) 段間は完全に等価だった。** 12 組すべてで差は 4 桁目以下:

    A_base η0.9 (b)0.4584/(c)0.4582   B_900C η0.9 (b)0.4068/(c)0.4068
    C_03MPaG η0.9 (b)0.4000/(c)0.4000  D_01MPaG η0.9 (b)0.4161/(c)0.4160

Cond1Gas と Cond2Gas の CO2 は絶対量が同一（R2 は CO2 を作りも消しもしない）なので、
どちらで抜いてもループへの効果は原理的に同じ。**カルボニル化床の CO 分圧を上げる
という (c) の狙いは実測でゼロだった。**ループ全体が痩せる影響に飲み込まれている。
差が出るのは塔の処理量だけ（Cond1Gas のほうが 4% 多い）。

**(2) CO2 を改質器へ戻す前提でも、plant4 では除去が効かない。**

改質器の CO2 供給は H2/CO=1.3837 の制約で総量が決まっているので、回収 CO2 を
戻しても改質器の解は不変（新鮮 CO2 の勘定だけが変わる）。その前提で計算しても:

    構成                        MA      C効率(捨)  C効率(戻)  触媒
    B_900C  baseline          0.9547     63.15%    63.15%   415 mL  ← plant4 の最良
    A_base  a_upstream η0.5   0.9748     53.45%    63.11%   424 mL
    A_base  b_recycle η0.9    0.4584     25.14%    44.21%   398 mL

最良は**除去なし**。plant3 では戻す前提で評価が反転する（D_01MPaG b_recycle η0.9 が
75.06% → 86.43%）が、plant4 では反転しない。MA の落ち込みが大きすぎて新鮮 CO2 の
節約で取り返せない。

**(3) plant4 が plant3 に勝つのは A_base だけ。**

    baseline の MA    plant3    plant4
    A_base            1.0416    1.0651   ← plant4 の勝ち
    B_900C            1.0036    0.9547
    C_03MPaG          0.9653    0.8005   ← 大きく負け
    D_01MPaG          0.9966    0.7963   ← 大きく負け

**真因は DME パージ。** plant4 は Cond2 で未反応 DME をあえてガス側に残して M1 へ
戻す設計なので、パージガスが DME に富む。ループ CO2 が少ないほど DME が濃縮され、
パージ率 5% がそれを大量に持ち出す。baseline の DME パージは A_base 0.054 に対し
D_01MPaG 0.439 と 8 倍。**plant4 の弱点は「CO2 希釈剤を失うこと」全般**であって、
CO2 除去に限らない。低 CO2 の改質条件でも同じ機構で悪化する。

新鮮 CO の行き先で数えると（供給 3.5498 に対し合計 3.5498 で厳密に閉じる）、
A_base 除去なし → (b) η0.9 の差は:

    酢酸メチル（3 CO/MA）  −1.8201
      DME としてパージ      +1.3320   損失の 73%  ← 主因
      順方向 WGS で CO2 へ  +0.6086   損失の 33%
      CO のままパージ       −0.1204   7% 戻る

仮説どおり WGS の寄与は plant3 の 80% から 33% に弱まった（Cond1 で H2O を抜く効果は
出ている）。だが DME という別の漏れ口が開き、差し引きで大幅に悪化した。

**結論: plant4 は CO2 除去の受け皿として不適。** 全構成の総合順位（投入炭素基準・
CO2 を改質器へ戻す前提）は plant3 + 循環除去が上位を占める:

    86.43%  plant3  D_01MPaG  b_recycle η0.9
    75.06%  plant3  D_01MPaG  baseline
    63.15%  plant4  B_900C    baseline        ← plant4 の最良

**次に効きそうなのは CO2 除去ではなく、パージからの DME 回収かパージ率の引き下げ。**
plant4 は除去なしでも DME を 0.05〜0.44 mol/h 捨てている（CO 換算で最大 0.88）。

────────────────────────────────────────────────────────────────────────
実行方法
────────────────────────────────────────────────────────────────────────
要 Cantera（改質器）と reaction_rate（速度論反応器）。
PYTHONPATH の "." は examples.* を絶対 import するため、"../reaction_rate/src" は
速度論 PFR のため。**どちらも必須**。

■ まず向きを見る（除去なし + (b) + (c) を η=0.9 で = 3 runs）

    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal_plant4.py \
        --cases baseline,b_recycle,c_interstage --eta 0.9 --solve-tol 1e-4 --stop-at-tol

■ 総当たり（改質条件 4 種 × 除去位置 3 種 × η 3 点 = 40 runs）

    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal_plant4.py --reform A_base   --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99 --solve-tol 1e-4 --stop-at-tol
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal_plant4.py --reform B_900C   --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99 --solve-tol 1e-4 --stop-at-tol
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal_plant4.py --reform C_03MPaG --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99 --solve-tol 1e-4 --stop-at-tol
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal_plant4.py --reform D_01MPaG --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99 --solve-tol 1e-4 --stop-at-tol

    python3 examples/merge_xlsx.py --all     # 最後に 3 つの xlsx を 1 つに結合

Windows (cmd) はこう:

    set PYTHONPATH=.;..\reaction_rate\src
    set OMP_NUM_THREADS=1
    set OPENBLAS_NUM_THREADS=1
    python -u examples\example_co2_removal_plant4.py --reform A_base --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99 --solve-tol 1e-4 --stop-at-tol

**スレッド数の固定について:** BLAS のスレッド数で浮動小数の加算順序が変わり、
悪条件な本問題では求解の経路が分岐する。ただし実測では総時間は 5% しか変わらず、
**遅さの対策としては --stop-at-tol と --x-scale jac（既定）のほうが桁違いに効く**。
詳細は plant3 版の docstring。

改質条件ごとに分けるのは仕様（1 回の比較の中で固定しないと「改質条件の効果」と
「除去位置の効果」が混ざる）。plant3 版と同じ理由。

**`-u` を付ける理由:** 付けないと stdout がブロックバッファになり、リダイレクト時に
進捗が何も見えなくなる（plant3 版で 52 分間無音になった事例あり）。
**`--solve-tol` について:** plant4 版はもともと既定 1e-4 なので明示は不要だが、
plant3 版と揃えて書いてある。1e-6 のような厳しい値にすると収束が渋いケースで
時間が倍になる（PFR の積分誤差で残差は 3e-6 で頭打ちになるため到達できない）。

■ 出力（examples/output/）

    co2removal_p4_{改質条件}_{ケース}_eta{η}_reformer.xlsx / .html
    co2removal_p4_{改質条件}_{ケース}_eta{η}_plant4.xlsx   / .html
    co2removal_p4_{改質条件}_{ケース}_eta{η}_caustic.xlsx  / .html

ラベルに "p4_" を付けて plant3 版の出力と混ざらないようにしてある。
xlsx は 3 つとも plant3 版と同じ MERGE_BASIS / MERGE_COMPONENTS で行構造を揃えて
あるので、merge_xlsx.py が列を連結するだけで 1 枚のストリーム表になる。

■ 所要時間（実測）

plant4 は plant3 より重い（変数 165〜177、触媒体積が 2 次元なので外側反復も 2 次元）。

    1 run        平均 13.5 分（最短 5.8 / 最長 24.7）
    総当たり 40 runs   **約 9 時間**
    まず向きを見る 3 runs  約 40 分

    参考: plant3 版は 28 runs で約 5 時間（1 run 平均 11 分）

**総当たりは半日仕事**なので、まず 3 runs で向きを見てから広げること。実際その
段取りで、plant4 の CO2 除去が −57% と桁違いに悪いことが 40 分で分かった
（docstring 冒頭の検証結果）。

1 ケースが収束判定を外しても自動で 1 回だけ緩和再試行するので、途中で止まらない。
失敗したケースは最後に一覧される。
"""

import argparse
import os
import time

import numpy as np

import examples.example_co2_removal as cr      # 吸収塔サブフローシートと出力設定を共有
import examples.example_pattern1_plant3 as pp  # 前工程（改質器）
import examples.example_plant4 as p4
from chemflow2 import export_mermaid, stream_table, to_excel

OUT = cr.OUT

#: 比較の土台にする改質条件（既定は A_base = 現状条件）。--reform で切り替える。
#: plant3 版と同じく、1 回の比較の中では固定すること。
BASE_CASE = dict(pp.CASES["A_base"])

#: plant4 のストリーム名。plant3 とは番号も名前も違うのでここに集約しておく。
#: （plant3: 7.Purge / 9.MethylAcetate / 5.CondGas / 14.CO2Vent）
S_FEED = "1. Steam1"
S_REACTOR_IN = "2. ReactorIn"
S_COND1_GAS = "4. Cond1Gas"          # R2 の入口（段間除去を入れると 17 に変わる）
S_COND2_GAS = "11. Cond2Gas"         # 循環側（パージ分岐の手前）
S_PURGE = "12. Purge"
S_MA = "15. MethylAcetate"
S_VENT = "18. CO2Vent"
S_SCRUBBED = "17. ScrubbedGas"

#: ケース定義: tag → (改質器側の除去率を使うか, plant4 の co2_position, 塔の名前, 説明)
CASES = {
    "baseline":      (False, "recycle",    None,    "除去なし"),
    "a_upstream":    (True,  "recycle",    "T-100", "(a) 前段: 改質器の下流"),
    "b_recycle":     (False, "recycle",    "T-101", "(b) 循環: Cond2Gas → パージ分岐"),
    "c_interstage":  (False, "interstage", "T-102", "(c) 段間: Cond1Gas → R2"),
}


def run_case(tag: str, eta: float, purge: float, reform: dict, label: str,
             solve_tol: float = 1e-4, progress_every: int = 0,
             stop_at_tol: bool = False, solver_tols: float = 1e-12,
             x_scale: str | float = "jac") -> dict:
    """1 ケース解いて、改質側 + plant4 側の指標を返す。"""
    use_upstream, position, tower, title = CASES[tag]
    print(f"\n{'=' * 78}\n=== {label}  {title}（η={eta:.0%}, パージ率={purge:.0%}）"
          f"\n{'=' * 78}")

    # --- 前工程（改質器）。(a) のときだけ改質器側に T-100 を置く ---
    ref_removal = eta if use_upstream else None
    rp, DryGas, CO2f, H2Of, PlantFeed = pp.solve_reformer(
        verbose=False, co2_removal=ref_removal, **reform)
    print(f"改質器: H2O={H2Of.flow_of('H2O'):.4f} CO2供給={CO2f.flow_of('CO2'):.4f} mol/h")
    print(stream_table(rp.streams, basis=["mol", "mole_frac"]))
    to_excel(rp.streams, os.path.join(OUT, f"co2removal_{label}_reformer.xlsx"),
             sheet=cr.sheet_name(f"{label}_reformer"),
             basis=cr.MERGE_BASIS, components=cr.MERGE_COMPONENTS)
    export_mermaid(rp, os.path.join(OUT, f"co2removal_{label}_reformer.html"),
                   title=f"Reformer + CO2 removal ({label})", style="diamond")

    if use_upstream:
        cr.report_caustic(label, {f: DryGas.flow_of(f) for f in DryGas.formulas},
                          eta, name=tower)

    # --- 後工程（plant4）。(b)/(c) のときだけ plant4 の中に塔を置く ---
    feed = {f: PlantFeed.flow_of(f) for f in PlantFeed.formulas}
    plant_removal = None if tag in ("baseline", "a_upstream") else eta
    problem, streams, v = p4.solve_with_sv(
        feed=feed, co2_removal=plant_removal, co2_position=position,
        purge=purge, solve_tol=solve_tol, progress_every=progress_every,
        stop_at_tol=stop_at_tol, solver_tols=solver_tols,
        solver_kwargs={"x_scale": x_scale})

    by = {s.name: s for s in problem.streams}
    Steam1, ReactorIn = by[S_FEED], by[S_REACTOR_IN]
    Purge, MA = by[S_PURGE], by[S_MA]
    # 段間に塔を入れると R2 の入口が Cond1Gas → ScrubbedGas に差し替わる
    r2_in = next(u.inlets[0] for u in problem.units if u.name == "R2 MA bed")

    print(stream_table(problem.streams, basis=["mol", "mole_frac"]))
    to_excel(problem.streams, os.path.join(OUT, f"co2removal_{label}_plant4.xlsx"),
             sheet=cr.sheet_name(f"{label}_plant4"),
             basis=cr.MERGE_BASIS, components=cr.MERGE_COMPONENTS)
    export_mermaid(problem, os.path.join(OUT, f"co2removal_{label}_plant4.html"),
                   title=f"plant4 + CO2 removal ({label})", style="diamond")

    # (b)/(c) の吸収塔は、それぞれ Cond2Gas / Cond1Gas を処理する
    if plant_removal is not None:
        src = by[S_COND1_GAS] if position == "interstage" else by[S_COND2_GAS]
        cr.report_caustic(label, {f: src.flow_of(f) for f in src.formulas},
                          eta, name=tower)

    n_in = float(ReactorIn.total_flow.eval())
    n_r2 = float(r2_in.total_flow.eval())
    ref_vent = next((s for s in rp.streams if s.name == "9. CO2Vent"), None)
    up_removed = ref_vent.flow_of("CO2") if ref_vent is not None else 0.0
    loop_removed = by[S_VENT].flow_of("CO2") if S_VENT in by else 0.0
    ma = MA.flow_of("CH3COOCH3")
    # plant4 ループの CO2 収支（出 − 入）。(a) の前段除去は plant4 の外なので入れない
    net_co2 = Purge.flow_of("CO2") + loop_removed - Steam1.flow_of("CO2")
    return {
        "供給 CO2 [mol/h]": Steam1.flow_of("CO2"),
        "除去 CO2 [mol/h]": up_removed + loop_removed,
        "反応器入口 [mol/h]": n_in,
        "反応器入口 CO2 [mol%]": ReactorIn.flow_of("CO2") / n_in * 100,
        "反応器入口 CO  [mol%]": ReactorIn.flow_of("CO") / n_in * 100,
        # (c) 段間除去の狙いはここ。カルボニル化床に入る CO 分圧を上げること
        "R2入口 CO  [mol%]": r2_in.flow_of("CO") / n_r2 * 100,
        "R2入口 CO2 [mol%]": r2_in.flow_of("CO2") / n_r2 * 100,
        "R2入口 H2O [mol%]": r2_in.flow_of("H2O") / n_r2 * 100,
        "触媒 V1 hybrid [mL]": v[0] * 1e6,
        "触媒 V2 MA床 [mL]": v[1] * 1e6,
        "触媒 合計 [mL]": (v[0] + v[1]) * 1e6,
        "酢酸メチル [mol/h]": ma,
        "CO パージ損失 [mol/h]": Purge.flow_of("CO"),
        "DME パージ [mol/h]": Purge.flow_of("CH3OCH3"),
        "CO2 正味生成 [mol/h]": net_co2,
        "新鮮CO基準 炭素収率 [%]": ma * 3 / Steam1.flow_of("CO") * 100,
        # 供給 CO = MA(3C) + CO パージ + 正味 CO2 + DME パージ(2C) で閉じるはず。
        # 閉じなければ指標の取り違えか、想定外の CO の行き先がある。
        "CO 収支計 [mol/h]": (3 * ma + Purge.flow_of("CO") + net_co2
                              + 2 * Purge.flow_of("CH3OCH3")),
        "供給 CO [mol/h]": Steam1.flow_of("CO"),
    }


def main():
    ap = argparse.ArgumentParser(
        description="CO2 除去位置 (a)/(b)/(c) の比較 — 後工程 plant4")
    ap.add_argument("--eta", default="0.90",
                    help="CO2 除去率。カンマ区切りで複数指定すると η スイープ")
    ap.add_argument("--cases", default="baseline,b_recycle,c_interstage",
                    help=f"実行するケース（{'/'.join(CASES)}）をカンマ区切りで")
    ap.add_argument("--reform", default="A_base",
                    help=f"土台にする改質条件（{'/'.join(pp.CASES)}）")
    ap.add_argument("--purge", type=float, default=p4.PURGE,
                    help=f"パージ率（既定 {p4.PURGE}）")
    ap.add_argument("--progress-every", type=int, default=200,
                    help="残差評価 N 回ごとに進捗を表示（既定 200、0 で無効）。"
                         "1 回の求解が数十分かかることがあり、無音だと計算中か"
                         "固まったか判別できないため")
    ap.add_argument("--stop-at-tol", action="store_true",
                    help="‖residual‖ が --solve-tol を下回った時点で打ち切る。"
                         "合否判定を満たした瞬間に止めるので、合格ラインをとうに"
                         "下回ってから延々と粘る無駄が無くなる（実測 7〜8 倍速）")
    ap.add_argument("--x-scale", default="jac",
                    help="least_squares の x_scale（既定 'jac'）。'jac' はヤコビアン各列の"
                         "ノルムの逆数から変数ごとのスケールを推定する（Moré 1978）。"
                         "本問題は流量が 25 桁にわたるため、既定の 1.0（全変数を同じ"
                         "物差しで測る）では trust region がうまく働かない。"
                         "'1.0' を渡すと従来の挙動")
    ap.add_argument("--solver-tols", type=float, default=1e-12,
                    help="least_squares の ftol/xtol/gtol に一括で入れる値（既定 1e-12）。"
                         "**実行時間を実際に決めているのはこれ**")
    ap.add_argument("--solve-tol", type=float, default=1e-4,
                    help="Problem.solve の残差ノルム判定（既定 1e-4）。**絶対値**なので"
                         "スケール依存。plant4 は PFR の積分誤差で残差が 3e-6 で"
                         "頭打ちになるため 1e-6 のような厳しい値には到達できない")
    ap.add_argument("--retry-tol", type=float, default=1e-3,
                    help="1 回目が収束判定を外したときの緩和 solve_tol（既定 1e-3）")
    args = ap.parse_args()

    if args.reform not in pp.CASES:
        ap.error(f"--reform は {list(pp.CASES)} のいずれか")
    reform = dict(pp.CASES[args.reform])
    etas = [float(e) for e in args.eta.split(",") if e.strip()]
    tags = [t.strip() for t in args.cases.split(",") if t.strip()]
    for t in tags:
        if t not in CASES:
            ap.error(f"--cases に未知のケース {t!r}（{list(CASES)} のいずれか）")

    # --x-scale は 'jac' か数値。数値なら float に直して渡す
    xs = args.x_scale if args.x_scale == "jac" else float(args.x_scale)

    os.makedirs(OUT, exist_ok=True)

    # baseline は η に依らないので 1 回だけ。その他は η ごとに回す。
    plan: list[tuple[str, float, str]] = []
    for tag in tags:
        if tag == "baseline":
            plan.append((tag, 0.0, f"p4_{args.reform}_baseline"))
        else:
            for eta in etas:
                plan.append((tag, eta, f"p4_{args.reform}_{tag}_eta{eta:g}"))

    print(f"実行計画: {len(plan)} runs（plant4 は 1 run 15〜30 分）")
    for _tag, eta, label in plan:
        print(f"  - {label}  η={eta:.0%}")

    # 1 ケースが収束判定を外しても 1 桁緩めて 1 回だけ再試行する（plant3 版と同じ）。
    # plant4 は既定 solve_tol が 1e-4（PFR の積分誤差で残差が 3e-6 で頭打ちになるため）
    # なので、緩和先は 1e-3。
    results, failed, retried = {}, [], []
    for i, (tag, eta, label) in enumerate(plan, 1):
        t0 = time.time()
        try:
            results[label] = run_case(tag, eta, args.purge, reform, label,
                                      solve_tol=args.solve_tol,
                                      progress_every=args.progress_every,
                                      stop_at_tol=args.stop_at_tol,
                                      solver_tols=args.solver_tols,
                                      x_scale=xs)
            print(f"[{i}/{len(plan)}] {label} 完了 ({time.time() - t0:.0f}s)")
        except Exception as e:                      # noqa: BLE001 — 打ち切らず続行
            print(f"[{i}/{len(plan)}] {label} 1回目失敗: {type(e).__name__}: {e}")
            print(f"    → solve_tol を {args.retry_tol:g} に緩めて再試行")
            try:
                results[label] = run_case(tag, eta, args.purge, reform, label,
                                          solve_tol=args.retry_tol,
                                          progress_every=args.progress_every,
                                          stop_at_tol=args.stop_at_tol,
                                          solver_tols=args.solver_tols,
                                          x_scale=xs)
                retried.append(label)
                print(f"[{i}/{len(plan)}] {label} 完了（緩和判定） "
                      f"({time.time() - t0:.0f}s)")
            except Exception as e2:                 # noqa: BLE001
                failed.append((label, f"{type(e2).__name__}: {e2}"))
                print(f"[{i}/{len(plan)}] {label} 失敗 ({time.time() - t0:.0f}s): "
                      f"{type(e2).__name__}: {e2}")

    if len(results) > 1:
        print(f"\n{'=' * 78}\n=== 除去位置の比較（plant4・改質条件 {args.reform} 固定・"
              f"パージ率 {args.purge:.0%}）\n{'=' * 78}")
        metrics = list(next(iter(results.values())))
        heads = {t: t[len(f"p4_{args.reform}") + 1:] or "baseline" for t in results}
        print(f"{'指標':>24s}" + "".join(f"{heads[t]:>18s}" for t in results))
        for m in metrics:
            print(f"{m:>24s}" + "".join(f"{results[t][m]:18.4f}" for t in results))
        print("\n判定: 酢酸メチルが「除去なし」より増えていれば、plant4 では CO2 除去が効く。\n"
              "      plant3 では 28 runs すべてで減った。plant4 は Cond1 で H2O を抜くぶん\n"
              "      順方向 WGS が弱まるはず、というのが仮説（docstring 参照）。\n"
              "      「CO2 正味生成」の増分が plant3 より小さければ、その仮説どおり。")

    if retried:
        print(f"\n※ 緩和判定（solve_tol={args.retry_tol:g}）で解いたケース {len(retried)} 件: "
              f"{', '.join(retried)}")
    if failed:
        print(f"\n⚠ 収束しなかったケース {len(failed)} 件:")
        for label, msg in failed:
            print(f"  - {label}: {msg}")


if __name__ == "__main__":
    main()
