r"""サンプル: 改質器 A〜D の DryGas を後工程 plant4（2段反応器）に流す。

example_reformer.py の 4 ケースを、そのまま plant3（速度論タンデム反応器）の新鮮供給に
する。**改質器へ戻る流れが無い一方向結合**なので逐次に解けば厳密。

    [reformer]  Feed(113 NL/h) → Gibbs(T,P) → 凝縮(25℃) → DryGas
    [plant4]    DryGas → M1 → R1(合成+脱水) → Cond1 ┬ ガス → R2 → Cond2 ┬ ガス → 循環
                                                     └ 液(H2O,MeOH)      └ 液 → 酢酸メチル

example_plant4.py は無改造で import して使う。plant4 は R1 と R2 の間で H2O/MeOH を
抜くので、H-MOR のカルボニル化を水の阻害から守れる。出力は
examples/output/reformer/ に分けてあり、既存の検討と混ざらない。

────────────────────────────────────────────────────────────────────────
何を見たいか
────────────────────────────────────────────────────────────────────────
改質器としては **D が圧倒的に良い**（改質 C 効率 55.8 → 93.7%、DryGas 中 CO2 が
25.3 → 1.5 mol%）。だが下流では逆に働く可能性がある。

過去の検討（example_co2_removal_plant4.py）で、**plant4 は「CO2 希釈剤を失うこと」に
弱い**ことが分かっている。plant4 は未反応 DME をあえて循環させる設計なので、ループ
CO2 が減ると DME が濃縮してパージで大量に失われた（損失の 73%）。D は DryGas の CO2 が
1.5 mol% しかないので、この機構で強く不利になる可能性が高い。

さらに plant4 が plant3 に勝ったのは A_base（CO2 が多い条件）だけで、低 CO2 の C/D 条件
では大きく負けていた（MA 0.80 対 0.97）。**今回の D は当時の C/D よりさらに CO2 が
少ない**ので、その傾向が強く出ると予想される。

一方で CO の絶対量は 41.9 → 82.6 NL/h と倍増する。**差し引きは解かないと決まらない**、
というのがこのスクリプトを回す理由。

指標は**投入炭素基準**で見ること。新鮮 CO 基準は改質器に入れた CO2 と CH4 を分母に
含めないので、改質条件を比較する用途には使えない（過去にこれで順位を取り違えた）。

────────────────────────────────────────────────────────────────────────
実行方法
────────────────────────────────────────────────────────────────────────
要 Cantera（改質器）と reaction_rate（速度論 PFR）。
**実測: 1 ケース平均 4.6 分（2.5〜7.2）、4 ケース計 23 分。**

    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_reformer_plant4.py
    PYTHONPATH=.:../reaction_rate/src python3 -u examples/example_reformer_plant4.py --cases A,D

Windows (cmd):

    set PYTHONPATH=.;..\reaction_rate\src
    python -u examples\example_reformer_plant4.py

**-u を付けること**（付けないと stdout がブロックバッファになり進捗が見えない）。
--stop-at-tol と x_scale='jac' は既定で有効にしてある（過去の検証で 691s → 204s）。
"""

import argparse
import os
import time

import examples.example_plant4 as p4
import examples.example_reformer as rf
from chemflow2 import export_mermaid, stream_table, to_excel

NV = rf.NV
OUT = rf.OUT

#: plant4 の 11 成分（plant3 と同一）
PLANT_C = p4.C
BASIS = rf.BASIS


def run_case(tag, *, solve_tol=1e-4, progress_every=200):
    """改質器 → plant3 を 1 ケース通し、指標を返す。"""
    print(f"\n{'=' * 78}\n=== {tag}: {rf.CASES[tag]['note']}\n{'=' * 78}", flush=True)

    # --- 前工程（改質器）---
    rpb, rst = rf.solve_case(tag)
    DryGas = rst[3]
    rsum = rf.summary(tag, rpb, rst)
    print(stream_table(rpb.streams, basis=["normal_volume", "mole_frac"]))
    to_excel(rpb.streams, os.path.join(OUT, f"p4_{tag}_reformer.xlsx"),
             sheet=f"p4_{tag}_reformer", basis=BASIS, components=PLANT_C)
    export_mermaid(rpb, os.path.join(OUT, f"p4_{tag}_reformer.html"),
                   title=f"Reformer {tag}", style="diamond")

    # --- 後工程（plant3）---
    t0 = time.time()
    problem, streams, v = p4.solve_with_sv(
        feed=rf.plant_feed(DryGas), solve_tol=solve_tol, stop_at_tol=True,
        progress_every=progress_every, solver_kwargs={"x_scale": "jac"})
    by = {s.name: s for s in problem.streams}
    Steam1, ReactorIn = by["1. Steam1"], by["2. ReactorIn"]
    # plant4 はストリーム名が plant3 と違う（12. Purge / 15. MethylAcetate）
    Purge, MA = by["12. Purge"], by["15. MethylAcetate"]
    r2_in = next(u.inlets[0] for u in problem.units if u.name == "R2 MA bed")
    print(stream_table(problem.streams, basis=["mol", "mole_frac"]))
    to_excel(problem.streams, os.path.join(OUT, f"p4_{tag}_plant4.xlsx"),
             sheet=f"p4_{tag}_plant3", basis=BASIS, components=PLANT_C)
    export_mermaid(problem, os.path.join(OUT, f"p4_{tag}_plant4.html"),
                   title=f"plant4 (feed = reformer {tag})", style="diamond")
    print(f"  → plant4 完了 ({time.time() - t0:.0f}s)", flush=True)

    n_in = float(ReactorIn.total_flow.eval())
    ma = MA.flow_of("CH3COOCH3")
    # 投入炭素 = 改質器に入れた CH4 + CO2（新鮮 CO 基準では改質器の消費が見えない）
    c_in = (rsum["供給 CH4 [NL/h]"] + rsum["供給 CO2 [NL/h]"]) / NV
    net_co2 = Purge.flow_of("CO2") - Steam1.flow_of("CO2")
    return {
        "改質 CH4 転化率 [%]": rsum["CH4 転化率 [%]"],
        "改質 CO2 転化率 [%]": rsum["CO2 転化率 [%]"],
        "改質 C 効率 [%]": rsum["改質 C 効率 [%]"],
        "DryGas H2/CO [-]": rsum["DryGas H2/CO [-]"],
        "DryGas [NL/h]": rsum["DryGas 合計 [NL/h]"],
        "DryGas CO2 [mol%]": rsum["DryGas CO2 [mol%]"],
        "投入炭素 [mol/h]": c_in,
        "反応器入口 [mol/h]": n_in,
        "反応器入口 CO2 [mol%]": ReactorIn.flow_of("CO2") / n_in * 100,
        "反応器入口 CO  [mol%]": ReactorIn.flow_of("CO") / n_in * 100,
        "R2入口 CO  [mol%]": r2_in.flow_of("CO") / float(r2_in.total_flow.eval()) * 100,
        # plant4 は触媒体積が 2 つ（ハイブリッド床とカルボニル化床）
        "触媒 V1 hybrid [mL]": v[0] * 1e6,
        "触媒 V2 MA床 [mL]": v[1] * 1e6,
        "触媒 合計 [mL]": (v[0] + v[1]) * 1e6,
        "酢酸メチル [mol/h]": ma,
        "CO パージ損失 [mol/h]": Purge.flow_of("CO"),
        "DME パージ [mol/h]": Purge.flow_of("CH3OCH3"),
        "CO2 正味生成 [mol/h]": net_co2,
        "新鮮CO基準 収率 [%]": ma * 3 / Steam1.flow_of("CO") * 100,
        "投入炭素基準 収率 [%]": ma * 3 / c_in * 100,
    }


def main():
    ap = argparse.ArgumentParser(description="改質器 A〜D → plant4")
    ap.add_argument("--cases", default=",".join(rf.CASES),
                    help=f"実行するケース（{'/'.join(rf.CASES)}）をカンマ区切りで")
    ap.add_argument("--solve-tol", type=float, default=1e-4)
    ap.add_argument("--progress-every", type=int, default=200)
    args = ap.parse_args()

    tags = [t.strip() for t in args.cases.split(",") if t.strip()]
    for t in tags:
        if t not in rf.CASES:
            ap.error(f"未知のケース {t!r}（{list(rf.CASES)} のいずれか）")
    os.makedirs(OUT, exist_ok=True)

    results, failed = {}, []
    for i, tag in enumerate(tags, 1):
        t0 = time.time()
        try:
            results[tag] = run_case(tag, solve_tol=args.solve_tol,
                                    progress_every=args.progress_every)
            print(f"[{i}/{len(tags)}] {tag} 完了 ({time.time() - t0:.0f}s)", flush=True)
        except Exception as e:                      # noqa: BLE001
            failed.append((tag, f"{type(e).__name__}: {e}"))
            print(f"[{i}/{len(tags)}] {tag} 失敗: {type(e).__name__}: {e}", flush=True)

    if len(results) > 1:
        print(f"\n{'=' * 78}\n=== 改質条件が plant4 に効く度合い\n{'=' * 78}")
        print(f"{'指標':>24s}" + "".join(f"{t:>14s}" for t in results))
        for k in next(iter(results.values())):
            print(f"{k:>24s}" + "".join(f"{results[t][k]:14.4f}" for t in results))
        print("\n読み方: **投入炭素基準 収率**で比較すること。新鮮 CO 基準は改質器に"
              "入れた\n        CH4・CO2 を分母に含めないので、改質条件の比較には使えない。")
    if failed:
        print(f"\n⚠ 失敗 {len(failed)} 件:")
        for tag, msg in failed:
            print(f"  - {tag}: {msg}")
    print(f"\n出力: {OUT}")


if __name__ == "__main__":
    main()
