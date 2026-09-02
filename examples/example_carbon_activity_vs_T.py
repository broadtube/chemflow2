r"""サンプル: 改質器 A〜D の炭素活量 a_C の温度依存性を描く。

example_reformer.py の 4 ケースについて、**横軸 = 温度、縦軸 = 炭素活量**の図を
1 ケース 1 枚の PNG で出力する。温度は 300 ℃ から各ケースの運転温度まで。

────────────────────────────────────────────────────────────────────────
3 本の曲線 — 同じ「温度」でも見ている場所が違う
────────────────────────────────────────────────────────────────────────
a_C は分圧と温度の両方で決まる。温度を下げるとき組成をどう扱うかで話が変わるので、
両方を 1 枚に描く。

  (1) **frozen（組成固定・実線）** — 改質器出口の組成をそのまま持って冷やす。
      廃熱ボイラや配管など、**触媒が無くて組成が凍結している下流**がこれ。
      K(T) だけが動くので 3 経路はばらける（気相が平衡していないため）。
      3 経路それぞれと、その**最大値**（保守側の判定値）を描く。

  (2) **equilibrium（各温度で再平衡・破線）** — その温度で改質器を運転したら、
      という仮定。触媒が効いている改質管の中がこれに当たる。5 種
      （CH4/CO2/H2O/CO/H2）のギブズ平衡なので **3 経路は厳密に一致**し 1 本になる。

  (3) **inlet feed（入口の供給ガス・点線）** — 予熱中の未反応の原料。
      **CO をまだ含まないので Boudouard と CO 還元は厳密に 0**、CH4 分解だけが
      効く。CH4 分解は吸熱なので **温度を上げるほど a_C は上がる**（frozen とは
      逆向き）。予熱温度を上げても入口の析出リスクは下がらない、が要点。

運転温度では出口組成そのものが平衡組成なので、**右端で (1) と (2) は必ず一致する**。

────────────────────────────────────────────────────────────────────────
読み方
────────────────────────────────────────────────────────────────────────
- Boudouard (2CO ⇌ C + CO2) は**発熱**なので低温ほど K が大きい。組成を固定して
  冷やすと a_C は跳ね上がり、frozen の最大値は低温側で桁で 1 を超える。
  → **出口ガスは「冷やす途中」が最も危ない**。ただしこれは熱力学的な可能性で、
    触媒の無い領域では速度が伴わない（村長 1961: 706 ℃以下では実験値は理論値より
    析出しにくい。references/PAPERS_INDEX.md §3）。
- equilibrium 側は逆に低温ほど a_C が下がる。低温では平衡がメタン化側に寄って
  CO が減るため。改質器そのものの析出リスクは高温・低酸化剤側で厳しい。

⚠️ **この図を「入口をここまで予熱せよ」と読んではいけない。** 実線は出口ガスの話で、
   入口の原料は点線のほう。入口は CO が無いので B・R 経路が存在せず、残る CH4 分解は
   低温ほど弱い（400〜500 ℃ で K_M = 0.066〜0.47、900 ℃ の 1/100 以下）。
   実際、千代田 WO2012140994A1 の実施例はいずれも **入口 500 ℃ / 出口 850 ℃**で、
   出口のカーボン活性を 0.90 に設計している（references/PAPERS_INDEX.md §0-2）。
   低温で入れて管内で昇温させるのが正しく、予熱を上げると点線が示すとおり逆効果。
   a_C ≤ 1 の判定は **平衡に達した出口ガス（equilibrated gas）に対する規約**であって、
   管内の各点に課すものではない。

要 Cantera + matplotlib:
    pip install chemflow2[gibbs] matplotlib
実行:
    PYTHONPATH=. python3 -u examples/example_carbon_activity_vs_T.py
    PYTHONPATH=. python3 -u examples/example_carbon_activity_vs_T.py --cases A,D
出力: examples/output/reformer/carbon_activity_vs_T_{A..D}.png

⚠️ 図中のラベルは英語（環境に CJK フォントが無いと豆腐になるため。
   example_carbon_limit.py と同じ方針）。
"""

import argparse
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from chemflow2.core.carbon_activity import (  # noqa: E402
    CARBON_ROUTES,
    carbon_activity,
    equilibrium_constants,
)
from chemflow2.core.pressure import parse_pressure  # noqa: E402
from examples.example_reformer import (  # noqa: E402
    CARBON_ACTIVITY_DESIGN,
    CASES,
    OUT,
    SPECIES,
    solve_case,
)

#: 描画する温度範囲 [℃]。上端は各ケースの運転温度。
T_MIN = 300.0
T_STEP = 5.0

# --- 配色: example_carbon_limit.py と揃える ---
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
ACCENT = "#eb6834"          # 運転点・a_C = 1 の注記色
ROUTE_COLOR = {"Boudouard": "#2a78d6", "CH4_cracking": "#7fb069", "CO_reduction": "#c9a227"}
ROUTE_LABEL = {"Boudouard": "B: 2CO $\\rightleftharpoons$ C + CO$_2$",
               "CH4_cracking": "M: CH$_4$ $\\rightleftharpoons$ C + 2H$_2$",
               "CO_reduction": "R: CO + H$_2$ $\\rightleftharpoons$ C + H$_2$O"}
MAXCOL = "#0b0b0b"          # frozen の最大値（判定に使う値）
EQCOL = "#8b2f8b"           # 各温度で再平衡
FEEDCOL = "#b5495b"         # 入口の供給ガス（未反応）


def _mole_fractions(stream):
    """Stream から 5 成分のモル分率を返す。持っていない成分は 0。

    供給ストリームは CH4/CO2/H2O/H2 の 4 成分で **CO を持たない**。この 0 が
    そのまま「入口では Boudouard と CO 還元が起こりようがない」を意味する。
    """
    total = float(stream.molar_flows.sum())
    have = set(stream.formulas)
    return {s: (stream.flow_of(s) / total if s in have else 0.0) for s in SPECIES}


def frozen_curves(x_frac, P_pascal, temps):
    """組成固定のまま温度を振ったときの経路別 a_C。{経路: 配列} を返す。"""
    P_atm = P_pascal / 101325.0
    p_atm = {s: x * P_atm for s, x in x_frac.items()}
    out = {r: np.empty(len(temps)) for r in CARBON_ROUTES}
    for i, T in enumerate(temps):
        routes = carbon_activity(p_atm, float(T), per_route=True)
        for r in CARBON_ROUTES:
            out[r][i] = routes[r]
    return out


def equilibrium_curve(x_frac, P_pascal, temps):
    """各温度で 5 種のギブズ平衡を取り直したときの a_C（3 経路は一致するので 1 本）。

    平衡させる元素組成は改質器出口と同じなので、供給を引き直さず出口組成から始める。
    """
    import cantera as ct

    selected = [s for s in ct.Species.list_from_file("gri30.yaml") if s.name in SPECIES]
    gas = ct.Solution(thermo="ideal-gas", species=selected)
    P_atm = P_pascal / 101325.0
    a = np.empty(len(temps))
    spread = np.empty(len(temps))
    for i, T in enumerate(temps):
        gas.TPX = float(T) + 273.15, P_pascal, x_frac
        gas.equilibrate("TP")
        p = {s: gas.X[gas.species_index(s)] * P_atm for s in SPECIES}
        routes = carbon_activity(p, float(T), per_route=True)
        v = np.array([routes[r] for r in CARBON_ROUTES])
        a[i] = v.max()
        spread[i] = (v.max() - v.min()) / v.max()      # 一致しているかの検算
    return a, spread


#: 管内の途中を表す「平衡への接近度」f。0 = 未反応のまま、1 = その温度で局所平衡。
APPROACH = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]


def approach_curves(feed_flows, P_pascal, temps, fractions=APPROACH):
    """管内の途中の a_C。{f: 配列} と、a_C = 1 を横切る温度 {f: T} を返す。

    **なぜ線形内挿でよいか**: 組成は反応進行度 ξ の 1 次式（n = n_feed + ν·ξ）なので、
    n_feed と局所平衡 n_eq(T) を f : 1−f で混ぜた点は「平衡の f 倍だけ反応が進んだ状態」
    そのものであり、元素収支も自動的に合う。f は接近度（approach to equilibrium）で、
    **f = 0 が「触媒が全く効いていない」、f = 1 が「触媒が十分速い」の両極**。
    実際の管内は温度が上がるにつれて f が 0 から 1 へ動く経路をたどる。
    """
    import cantera as ct

    selected = [s for s in ct.Species.list_from_file("gri30.yaml") if s.name in SPECIES]
    gas = ct.Solution(thermo="ideal-gas", species=selected)
    idx = {s: gas.species_index(s) for s in SPECIES}
    n_feed = np.array([feed_flows.get(s, 0.0) for s in SPECIES])
    total_feed = n_feed.sum()

    out = {f: np.empty(len(temps)) for f in fractions}
    for i, T in enumerate(temps):
        gas.TPX = float(T) + 273.15, P_pascal, {s: n_feed[k] / total_feed
                                                for k, s in enumerate(SPECIES)}
        gas.equilibrate("TP")
        # モル分率 → モル数（元素が保存するよう、全モル数は水素原子数などから決まる）。
        # ここは比だけ要るので、平衡のモル分率を feed の元素量に合わせて規格化する。
        x_eq = np.array([gas.X[idx[s]] for s in SPECIES])
        # C 原子基準で規格化すると全モル数が決まる（C は CH4 + CO + CO2 に入る）
        c_feed = n_feed[SPECIES.index("CH4")] + n_feed[SPECIES.index("CO")] \
            + n_feed[SPECIES.index("CO2")]
        c_frac = x_eq[SPECIES.index("CH4")] + x_eq[SPECIES.index("CO")] \
            + x_eq[SPECIES.index("CO2")]
        n_eq = x_eq * (c_feed / c_frac)
        for f in fractions:
            n = n_feed + f * (n_eq - n_feed)
            x = n / n.sum()
            p_atm = {s: x[k] * P_pascal / 101325.0 for k, s in enumerate(SPECIES)}
            out[f] = out[f]
            out[f][i] = carbon_activity(p_atm, float(T))

    # a_C > 1 になる温度の窓。曲線は f によって増加にも減少にもなるので、
    # 「何 ℃ 以上なら安全」ではなく**窓の下端と上端**で報告する。
    risk = {}
    for f in fractions:
        above = np.where(out[f] > 1.0)[0]
        risk[f] = (float(temps[above.min()]), float(temps[above.max()])) if len(above) else None
    return out, risk


def draw_tube(tag, temps, paths, risk, T_op, P_label, a_op, path):
    """管内の途中（接近度 f を振った束）を描いて保存する。"""
    fig, ax = plt.subplots(figsize=(8.4, 6.4), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.axhspan(1.0, 1e12, color=ACCENT, alpha=0.055, zorder=1, lw=0)
    ax.axhline(1.0, color=ACCENT, lw=1.4, zorder=3, label="$a_C = 1$  (carbon forms above)")

    fracs = sorted(paths)
    ramp = plt.get_cmap("viridis")
    for k, f in enumerate(fracs):
        col = FEEDCOL if f == 0.0 else (EQCOL if f == 1.0 else ramp(0.15 + 0.6 * k / len(fracs)))
        lw = 2.4 if f in (0.0, 1.0) else 1.4
        ls = (0, (1, 1.6)) if f == 0.0 else ("--" if f == 1.0 else "-")
        lab = ("f = 0.00  (feed, no reaction yet)" if f == 0.0 else
               "f = 1.00  (local equilibrium)" if f == 1.0 else f"f = {f:.2f}")
        ax.plot(temps, paths[f], color=col, lw=lw, ls=ls, zorder=5, label=lab)

    ax.plot([T_op], [a_op], marker="o", ms=9, color=ACCENT,
            markeredgecolor=SURFACE, markeredgewidth=2.0, zorder=7)
    ax.annotate(f"outlet {T_op:.0f} °C,  $a_C$ = {a_op:.4f}",
                xy=(0.985, 0.955), xycoords="axes fraction", color=INK,
                fontsize=10, fontweight="bold", ha="right", va="top", zorder=7)

    ax.set_yscale("log")
    ax.set_xlim(T_MIN, T_op)
    ax.set_xlabel("Local gas temperature along the tube  [°C]", color=INK2, fontsize=10.5)
    ax.set_ylabel("Carbon activity  $a_C$  [-]", color=INK2, fontsize=10.5)
    ax.set_title(f"Case {tag} — inside the tube: approach to equilibrium   "
                 f"({T_op:.0f} °C, {P_label})", color=INK, fontsize=12, pad=14, loc="left")
    ax.grid(True, which="major", color=GRID, lw=0.8, zorder=0)
    ax.grid(True, which="minor", color=GRID, lw=0.4, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9.5)
    leg = ax.legend(loc="lower left", frameon=True, fontsize=9, labelcolor=INK2,
                    borderaxespad=0.8, facecolor=SURFACE, edgecolor="none", framealpha=0.92,
                    title="approach to equilibrium", title_fontsize=9)
    leg.get_title().set_color(MUTED)
    leg.set_zorder(8)

    ax.text(0.0, -0.125,
            "f is how far the local gas has gone toward equilibrium at that temperature "
            "(0 = untouched feed, 1 = fully equilibrated).\n"
            "The real gas starts near f = 0 at the inlet and reaches f = 1 at the outlet, "
            "so the tube travels through this band.",
            transform=ax.transAxes, color=INK2, fontsize=8.6, ha="left", va="top",
            linespacing=1.5)

    fig.tight_layout(rect=(0, 0.085, 1, 1))
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)
    return path


def draw(tag, temps, frozen, a_max, a_eq, a_feed, T_op, P_label, a_op, path):
    """1 ケース分の図を描いて保存する。"""
    fig, ax = plt.subplots(figsize=(8.4, 6.4), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # --- a_C = 1（析出の熱力学的境界）と設計値 0.90 ---
    # この 2 本は**凡例で説明する**。y 軸が 12 桁の対数なので 1.0 と 0.90 は
    # ほぼ同じ高さに来てしまい、線の傍に文字を置くと必ずどれかの曲線と重なる。
    ax.axhline(1.0, color=ACCENT, lw=1.4, ls="-", zorder=3,
               label="$a_C = 1$  (carbon forms above)")
    ax.axhline(CARBON_ACTIVITY_DESIGN, color=MUTED, lw=1.0, ls=":", zorder=3,
               label=f"design value {CARBON_ACTIVITY_DESIGN:.2f}")

    # --- frozen: 最大値は**帯**で下に敷く ---
    # 最大値は必ずどれか 1 経路と重なる（この系では常に Boudouard）ので、細線で
    # 描くと下の経路線を完全に隠してしまう。太い半透明の帯にして経路線を透かす。
    ax.plot(temps, a_max, color=MAXCOL, lw=6.5, alpha=0.17, zorder=3.5,
            solid_capstyle="round",
            label="frozen — max of 3 (the value used for the check)")
    for r in CARBON_ROUTES:
        ax.plot(temps, frozen[r], color=ROUTE_COLOR[r], lw=1.8, ls="-",
                alpha=0.95, zorder=5, label=f"frozen — {ROUTE_LABEL[r]}")
    ax.plot(temps, a_eq, color=EQCOL, lw=2.2, ls="--", zorder=5,
            label="re-equilibrated at each T (3 routes coincide)")
    # 入口の供給ガス（未反応）。CO が無いので M 経路だけが効き、**高温ほど上がる**。
    ax.plot(temps, a_feed, color=FEEDCOL, lw=2.2, ls=(0, (1, 1.6)), zorder=5,
            label="inlet feed, unconverted (only CH$_4$ cracking is possible)")

    # --- 運転点。注記は曲線と重ならない右上の空きに置く ---
    ax.plot([T_op], [a_op], marker="o", ms=9, color=ACCENT,
            markeredgecolor=SURFACE, markeredgewidth=2.0, zorder=7)
    ax.annotate(f"outlet {T_op:.0f} °C,  $a_C$ = {a_op:.4f}",
                xy=(0.985, 0.955), xycoords="axes fraction", color=INK,
                fontsize=10, fontweight="bold", ha="right", va="top", zorder=7)

    ax.set_yscale("log")
    ax.set_xlim(T_MIN, T_op)
    ax.set_xlabel("Temperature  [°C]", color=INK2, fontsize=10.5)
    ax.set_ylabel("Carbon activity  $a_C$  [-]", color=INK2, fontsize=10.5)
    ax.set_title(f"Case {tag} — carbon activity vs. temperature   "
                 f"({T_op:.0f} °C, {P_label}, reformer outlet gas)",
                 color=INK, fontsize=12, pad=14, loc="left")
    ax.grid(True, which="major", color=GRID, lw=0.8, zorder=0)
    ax.grid(True, which="minor", color=GRID, lw=0.4, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9.5)
    leg = ax.legend(loc="lower left", frameon=True, fontsize=9,
                    labelcolor=INK2, borderaxespad=0.8,
                    facecolor=SURFACE, edgecolor="none", framealpha=0.92)
    leg.set_zorder(8)

    ax.text(0.0, -0.125,
            "Solid  = OUTLET gas, composition frozen, then cooled (no catalyst downstream).\n"
            "Dashed = gas re-equilibrated at each temperature.\n"
            "Dotted = INLET feed as charged \u2014 no CO yet, so only CH$_4$ cracking; "
            "preheating makes it worse, not better.",
            transform=ax.transAxes, color=INK2, fontsize=8.6, ha="left", va="top",
            linespacing=1.5)

    fig.tight_layout(rect=(0, 0.085, 1, 1))
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)
    return path


def run_case(tag):
    """1 ケース解いて図を出す。数値も print する。"""
    cfg = CASES[tag]
    print(f"\n=== ケース {tag}: {cfg['T']:.0f} ℃ / {cfg['P']} — {cfg['note']} ===",
          flush=True)
    _pb, streams = solve_case(tag)
    ReactOut = streams[2]
    x_frac = _mole_fractions(ReactOut)
    P_pascal = parse_pressure(cfg["P"])
    T_op = float(cfg["T"])

    n = int(round((T_op - T_MIN) / T_STEP)) + 1
    temps = np.linspace(T_MIN, T_op, n)

    frozen = frozen_curves(x_frac, P_pascal, temps)
    a_max = np.maximum.reduce([frozen[r] for r in CARBON_ROUTES])
    a_eq, spread = equilibrium_curve(x_frac, P_pascal, temps)

    # 入口の供給ガス（未反応・組成固定）。CO が無いので B と R は厳密に 0 で、
    # M（CH4 分解）だけが残る。M は吸熱なので**予熱するほど a_C は上がる**。
    x_feed = _mole_fractions(streams[0])
    feed_routes = frozen_curves(x_feed, P_pascal, temps)
    a_feed = np.maximum.reduce([feed_routes[r] for r in CARBON_ROUTES])

    print("  出口組成 [mol%]: " +
          "  ".join(f"{s} {x*100:6.3f}" for s, x in x_frac.items()))
    print(f"  平衡側の 3 経路の相対ばらつき: 最大 {spread.max():.2e}（0 なら厳密一致）")
    print("  供給組成 [mol%]: " +
          "  ".join(f"{s} {x*100:6.3f}" for s, x in x_feed.items() if x > 0))
    print(f"  {'T[℃]':>7s}{'frozen max':>12s}{'Boudouard':>12s}"
          f"{'CH4分解':>11s}{'CO還元':>11s}{'再平衡':>12s}{'入口供給':>12s}")
    for T in (300, 400, 500, 600, 700, 800, T_op):
        if T < T_MIN or T > T_op:
            continue
        i = int(np.argmin(np.abs(temps - T)))
        print(f"  {temps[i]:7.0f}{a_max[i]:12.4g}{frozen['Boudouard'][i]:12.4g}"
              f"{frozen['CH4_cracking'][i]:11.4g}{frozen['CO_reduction'][i]:11.4g}"
              f"{a_eq[i]:12.4g}{a_feed[i]:12.4g}")

    # frozen 側が a_C = 1 を切る温度（これより下は熱力学的に析出域）
    below = np.where(a_max <= 1.0)[0]
    if len(below) and len(below) < len(temps):
        T_cross = temps[below.min()]
        print(f"  frozen の a_C が 1 を下回るのは {T_cross:.0f} ℃ 以上"
              f"（それ未満は熱力学的に析出域）")
    elif not len(below):
        print("  frozen の a_C はこの範囲で常に 1 を超える")

    # --- 管内の途中（接近度 f を振る）---
    feed_flows = {c: streams[0].flow_of(c) for c in streams[0].formulas}
    paths, risk = approach_curves(feed_flows, P_pascal, temps)
    print("  管内の途中 — a_C > 1 になる温度の窓（熱力学的に析出しうる範囲）:")
    for f in sorted(paths):
        w = risk[f]
        where = f"{w[0]:.0f} 〜 {w[1]:.0f} ℃" if w is not None else "なし（常に 1 以下）"
        print(f"    f = {f:4.2f}  →  {where}"
              f"   （最大 {paths[f].max():.3g} @ {temps[int(np.argmax(paths[f]))]:.0f} ℃）")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"carbon_activity_vs_T_{tag}.png")
    draw(tag, temps, frozen, a_max, a_eq, a_feed, T_op, cfg["P"], float(a_eq[-1]), path)
    print(f"  → {path}")
    path2 = os.path.join(OUT, f"carbon_activity_in_tube_{tag}.png")
    draw_tube(tag, temps, paths, risk, T_op, cfg["P"], float(a_eq[-1]), path2)
    print(f"  → {path2}")
    return [path, path2]


def main():
    ap = argparse.ArgumentParser(description="改質器 A〜D の a_C 温度依存グラフ")
    ap.add_argument("--cases", default=",".join(CASES),
                    help=f"実行するケース（{'/'.join(CASES)}）をカンマ区切りで")
    args = ap.parse_args()
    paths = [q for t in args.cases.split(",") if t.strip() for q in run_case(t.strip())]
    print("\n出力:")
    for p in paths:
        print(f"  {p}")


if __name__ == "__main__":
    main()
