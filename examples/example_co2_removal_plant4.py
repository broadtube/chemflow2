"""サンプル: CO2 除去位置の比較 — 後工程を **plant4（2段反応器）** にした版。

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
plant3 では **CO2 除去は 26 runs すべてで酢酸メチルを減らした**（例外なし）。
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

**(c) 段間は (b) と実質同じだった**（MA 0.4582 対 0.4584）。カルボニル化床に入る
CO 分圧を上げる狙いだったが、ループ全体が痩せる効果に打ち消されて差が出ない。
plant4 固有の位置を作った意味は、少なくともこの条件では無かった。

**得られた知見:** plant4 は除去なしなら plant3 より良い（MA 1.0651 対 1.0416、
炭素収率 90.01% 対 88.03%）。水を抜いてカルボニル化を守る設計は効いている。
**しかし CO2 除去に対しては plant3 より脆い。** 次に効きそうなのは CO2 除去ではなく
**パージから DME を回収する**か**パージ率を下げる**方向。

⚠ 総当たりに広げる価値は薄い。上の 3 runs で −57% と桁違いに悪いので、η や改質条件を
振っても向きは変わらない見込み。回すとしても知見の確認用と割り切ること。

────────────────────────────────────────────────────────────────────────
実行方法
────────────────────────────────────────────────────────────────────────
要 Cantera（改質器）と reaction_rate（速度論反応器）。
PYTHONPATH の "." は examples.* を絶対 import するため、"../reaction_rate/src" は
速度論 PFR のため。**どちらも必須**。

■ まず向きを見る（除去なし + (b) + (c) を η=0.9 で = 3 runs）

    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal_plant4.py \
        --cases baseline,b_recycle,c_interstage --eta 0.9

■ 総当たり（改質条件 4 種 × 除去位置 3 種 × η 3 点 = 40 runs）

    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal_plant4.py --reform A_base   --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal_plant4.py --reform B_900C   --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal_plant4.py --reform C_03MPaG --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_co2_removal_plant4.py --reform D_01MPaG --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99

    python3 examples/merge_xlsx.py --all     # 最後に 3 つの xlsx を 1 つに結合

改質条件ごとに分けるのは仕様（1 回の比較の中で固定しないと「改質条件の効果」と
「除去位置の効果」が混ざる）。plant3 版と同じ理由。

■ 出力（examples/output/）

    co2removal_p4_{改質条件}_{ケース}_eta{η}_reformer.xlsx / .html
    co2removal_p4_{改質条件}_{ケース}_eta{η}_plant4.xlsx   / .html
    co2removal_p4_{改質条件}_{ケース}_eta{η}_caustic.xlsx  / .html

ラベルに "p4_" を付けて plant3 版の出力と混ざらないようにしてある。
xlsx は 3 つとも plant3 版と同じ MERGE_BASIS / MERGE_COMPONENTS で行構造を揃えて
あるので、merge_xlsx.py が列を連結するだけで 1 枚のストリーム表になる。

■ 所要時間

plant4 は plant3 より重い（変数 165〜177、触媒体積が 2 次元なので外側反復も 2 次元）。
実測で外側反復 1 回あたり約 240 秒。1 run あたり 15〜30 分を見込むこと。
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
             solve_tol: float = 1e-4) -> dict:
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
        purge=purge, solve_tol=solve_tol)

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
            results[label] = run_case(tag, eta, args.purge, reform, label)
            print(f"[{i}/{len(plan)}] {label} 完了 ({time.time() - t0:.0f}s)")
        except Exception as e:                      # noqa: BLE001 — 打ち切らず続行
            print(f"[{i}/{len(plan)}] {label} 1回目失敗: {type(e).__name__}: {e}")
            print(f"    → solve_tol を {args.retry_tol:g} に緩めて再試行")
            try:
                results[label] = run_case(tag, eta, args.purge, reform, label,
                                          solve_tol=args.retry_tol)
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
              "      plant3 では 26 runs すべてで減った。plant4 は Cond1 で H2O を抜くぶん\n"
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
