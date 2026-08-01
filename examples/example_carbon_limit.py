"""サンプル: 炭素析出限界線図（carbon limit diagram）を自前の定義で描く。

千代田化工建設 Mikuriya, Yagi, Shimura, AIChE Annual Meeting 2008,
"Development of a new Synthesis Gas Production Catalyst and Process with CO2 and
H2O Reforming" の Fig.4 と同じ図を、こちらの炭素活量の定義で計算して描画する。
（`references/PAPERS_INDEX.md` §1-1 参照。原論文は等高線 "Carbon Activity=1" を
示すだけで定義式を書いていないため、下記の教科書的定義で再現できるかを確かめる。）

炭素活量 a_C: その気相と平衡にある固体炭素（黒鉛）の活量。基準状態は純黒鉛で a=1。
**a_C > 1 なら気相が黒鉛に対して過飽和 = 炭素析出が熱力学的に有利。**

    B  2CO  ⇌ C(gr) + CO2    a_C = K_B · p_CO² / p_CO2
    M  CH4  ⇌ C(gr) + 2H2    a_C = K_M · p_CH4 / p_H2²
    R  CO+H2 ⇌ C(gr) + H2O   a_C = K_R · p_CO·p_H2 / p_H2O

3 式の差は水性ガスシフト（B−R）と水蒸気改質（M−R）なので、**気相が平衡していれば
3 つの a_C は一致する**。改質器出口は Gibbs 平衡なので一意に決まる。非平衡の気相に
当てる場合に備えて最大値を取る（保守側）。

分圧は atm 基準（Cantera の reference_pressure = 101325 Pa）。K は標準生成ギブズ
エネルギーから K = exp(−ΔG°/RT) で計算する。気相 gri30.yaml、固体炭素 graphite.yaml。

線図の座標: 原子比 O/C（横軸）と H/C（縦軸）。CH4 + CO2 + H2O の供給に対して
    C = CH4 + CO2,  H = 4·CH4 + 2·H2O,  O = 2·CO2 + H2O
なので、CH4+CO2 のみなら (0,4)→(2,0) の直線、CH4+H2O のみなら (0,4) から傾き 2 の
直線になり、この 2 本が物理的に到達可能な領域の左側境界（V 字）を成す。

要 Cantera + matplotlib:
    pip install chemflow2[gibbs] matplotlib
実行: python3 examples/example_carbon_limit.py

⚠️ 図中のラベルは英語。環境に CJK フォントが無いと日本語が豆腐になるため
   （原論文の Fig.4 も英語なので比較にも都合が良い）。
"""

import os

import cantera as ct
import matplotlib
import numpy as np
from scipy.optimize import brentq

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SPECIES = ["CH4", "CO", "CO2", "H2", "H2O"]

# --- 配色: 温度は順序量なので単一色相の逐次ランプ（青 250→650）を使う ---
TEMPS = [600, 700, 800, 900, 950]
TEMP_COLOR = {600: "#86b6ef", 700: "#5598e7", 800: "#2a78d6",
              900: "#1c5cab", 950: "#104281"}
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
ACCENT = "#eb6834"      # 運転点マーカー（系列色ではなく注記色）

_gas_th = ct.Solution("gri30.yaml")
_graphite = ct.Solution("graphite.yaml")
_gas_eq = ct.Solution("gri30.yaml")


def equilibrium_constants(T_kelvin: float) -> tuple[float, float, float]:
    """(K_B, K_M, K_R) を標準生成ギブズエネルギーから返す（1 atm 基準）。"""
    _gas_th.TP = T_kelvin, _gas_th.reference_pressure
    _graphite.TP = T_kelvin, _graphite.reference_pressure
    g = {s: _gas_th.standard_gibbs_RT[_gas_th.species_index(s)] for s in SPECIES}
    g["C"] = _graphite.standard_gibbs_RT[0]
    return (np.exp(-((g["C"] + g["CO2"]) - 2 * g["CO"])),
            np.exp(-((g["C"] + 2 * g["H2"]) - g["CH4"])),
            np.exp(-((g["C"] + g["H2O"]) - g["CO"] - g["H2"])))


def carbon_activity(p_atm: dict[str, float], T_kelvin: float) -> float:
    """分圧 [atm] の気相に対する炭素活量（3 経路の最大値）。"""
    K_B, K_M, K_R = equilibrium_constants(T_kelvin)
    eps = 1e-300
    return max(
        K_B * p_atm["CO"] ** 2 / max(p_atm["CO2"], eps),
        K_M * p_atm["CH4"] / max(p_atm["H2"] ** 2, eps),
        K_R * p_atm["CO"] * p_atm["H2"] / max(p_atm["H2O"], eps),
    )


def feed_from_ratios(oc: float, hc: float) -> tuple[float, float, float] | None:
    """元素比 (O/C, H/C) → (CH4, CO2, H2O)（C=1 に規格化）。到達不能なら None。"""
    ch4 = (hc - 2 * oc + 4) / 8.0
    co2 = 1.0 - ch4
    h2o = oc - 2.0 + 2 * ch4
    if min(ch4, co2, h2o) < -1e-10:
        return None
    return max(ch4, 0.0), max(co2, 0.0), max(h2o, 0.0)


def equilibrium_carbon_activity(oc: float, hc: float, T_celsius: float,
                                P_pascal: float) -> float | None:
    """元素比 (O/C, H/C) の気体を T, P で平衡化し、その平衡気相の a_C を返す。"""
    feed = feed_from_ratios(oc, hc)
    if feed is None:
        return None
    ch4, co2, h2o = feed
    _gas_eq.TPX = T_celsius + 273.15, P_pascal, {"CH4": ch4, "CO2": co2, "H2O": h2o}
    _gas_eq.equilibrate("TP")
    P_atm = P_pascal / 101325.0
    p = {s: _gas_eq.X[_gas_eq.species_index(s)] * P_atm for s in SPECIES}
    return carbon_activity(p, T_celsius + 273.15)


def carbon_limit_oc(hc: float, T_celsius: float, P_pascal: float) -> float | None:
    """その H/C・T・P で a_C = 1 になる O/C。析出しないなら None。"""
    lo = abs(4.0 - hc) / 2.0 + 1e-4            # 組成が非負になる下限（V 字の境界）
    hi = min(3.0, (hc + 4.0) / 2.0 - 1e-4)     # CH4 = 0 の上限、図の右端で切る
    if lo >= hi:
        return None
    a_lo = equilibrium_carbon_activity(lo, hc, T_celsius, P_pascal)
    if a_lo is None or a_lo < 1.0:
        return None                            # 最も酸素の少ない端でも析出しない
    try:
        return brentq(lambda x: equilibrium_carbon_activity(x, hc, T_celsius, P_pascal) - 1.0,
                      lo, hi, xtol=1e-7)
    except ValueError:
        return None


def limit_curve(T_celsius: float, P_pascal: float, hc_max: float = 9.0, n: int = 130):
    """1 本の等温線（a_C = 1 の軌跡）を (O/C 配列, H/C 配列) で返す。"""
    oc, hc = [], []
    for h in np.linspace(0.05, hc_max, n):
        x = carbon_limit_oc(float(h), T_celsius, P_pascal)
        if x is not None:
            oc.append(x)
            hc.append(float(h))
    return np.array(oc), np.array(hc)


def draw(P_pascal: float, path: str, title: str, operating_point=None):
    """炭素析出限界線図を描いて PNG に保存する。"""
    fig, ax = plt.subplots(figsize=(8.4, 6.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # --- 到達可能領域の境界（V 字）: 純 CH4+CO2 と 純 CH4+H2O ---
    ax.plot([0, 2], [4, 0], color=MUTED, lw=1.2, zorder=2)
    ax.plot([0, 3], [4, 10], color=MUTED, lw=1.2, zorder=2)
    ax.text(0.95, 0.55, "CH$_4$ + CO$_2$ only", color=MUTED, fontsize=8.5,
            rotation=-31, ha="center", va="center")
    ax.text(1.35, 6.95, "CH$_4$ + H$_2$O only", color=MUTED, fontsize=8.5,
            rotation=48, ha="center", va="center")

    # --- 等温線（a_C = 1）---
    # 直接ラベルは各曲線の**最も右に張り出す点**に置く。終端に置くと 5 本とも
    # CH4+H2O 境界の一点に集まって重なるため。
    for T in TEMPS:
        oc, hc = limit_curve(T, P_pascal)
        if len(oc) == 0:
            continue
        ax.plot(oc, hc, color=TEMP_COLOR[T], lw=2.0, zorder=4,
                label=f"{T} °C", solid_capstyle="round")
        i = int(np.argmax(oc))
        ax.annotate(f"{T} °C", xy=(oc[i], hc[i]), xytext=(9, 0),
                    textcoords="offset points", color=TEMP_COLOR[T],
                    fontsize=9, fontweight="bold", ha="left", va="center", zorder=6)

    # --- 領域の意味（境界は温度ごとに違うので、位置固定の領域ラベルは置かない）---
    ax.text(0.02, -0.135,
            "Left of a curve: $a_C > 1$, carbon forms.   Right: $a_C \\leq 1$, carbon-free.",
            transform=ax.transAxes, color=INK2, fontsize=9.5, ha="left", va="top")

    # --- 運転点 ---
    if operating_point is not None:
        oc0, hc0, lbl = operating_point
        ax.plot([oc0], [hc0], marker="o", ms=9, color=ACCENT,
                markeredgecolor=SURFACE, markeredgewidth=2.0, zorder=7)
        ax.annotate(lbl, xy=(oc0, hc0), xytext=(12, -14), textcoords="offset points",
                    color=INK, fontsize=9.5, fontweight="bold", zorder=7)

    ax.set_xlim(0, 3)
    ax.set_ylim(0, 10)
    ax.set_xlabel("O / C  atomic ratio", color=INK2, fontsize=10.5)
    ax.set_ylabel("H / C  atomic ratio", color=INK2, fontsize=10.5)
    ax.set_title(title, color=INK, fontsize=12, pad=14, loc="left")
    ax.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9.5)
    leg = ax.legend(title="Outlet temperature", loc="upper left", frameon=False,
                    fontsize=9.5, title_fontsize=9.5, labelcolor=INK2)
    leg.get_title().set_color(MUTED)

    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)
    return path


def main():
    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)

    # --- (1) 原論文 Fig.4 と同条件（1.5 MPa）— 照合用 ---
    p1 = draw(1.5e6, os.path.join(out, "carbon_limit_1.5MPa.png"),
              "Carbon limit diagram  ($a_C = 1$)   —   1.5 MPa")
    print(f"図(1) 原論文 Fig.4 と同条件: {p1}")

    # --- (2) pattern1 の実運転圧力 ---
    CH4, CO2, H2O, H2 = 2.137257, 3.333756, 2.776607, 0.885957
    C = CH4 + CO2
    oc0 = (2 * CO2 + H2O) / C
    hc0 = (4 * CH4 + 2 * H2O + 2 * H2) / C
    P_op = 1.04e6 + 101325.0
    p2 = draw(P_op, os.path.join(out, "carbon_limit_1.04MPaG.png"),
              "Carbon limit diagram  ($a_C = 1$)   —   1.04 MPaG (1.1413 MPa abs)",
              operating_point=(oc0, hc0, "pattern1"))
    print(f"図(2) pattern1 の運転圧力  : {p2}")

    # --- 数値も出す ---
    print(f"\npattern1 の運転点: O/C = {oc0:.4f}, H/C = {hc0:.4f}")
    print(f"{'T[°C]':>7s}{'a_C':>10s}{'限界O/C':>10s}  判定")
    for T in (800, 850, 900):
        a = equilibrium_carbon_activity(oc0, hc0, T, P_op)
        lim = carbon_limit_oc(hc0, T, P_op)
        judge = "炭素フリー" if a <= 1.0 else "析出域"
        print(f"{T:7d}{a:10.4f}{lim:10.3f}  {judge}（限界より "
              f"{oc0 - lim:+.3f} 右）")

    # --- 原論文 Fig.4 の読み取り値との照合（1.5 MPa）---
    print(f"\n--- 原論文 Fig.4 との照合（1.5 MPa, a_C=1 の O/C）---")
    read = {4.0: {950: 1.007, 900: 1.131, 800: 1.366, 700: 1.710, 600: 2.000},
            2.0: {950: 1.159, 900: 1.255, 800: 1.600, 700: 1.959, 600: 2.207}}
    print(f"{'H/C':>5s}{'T[°C]':>7s}{'本実装':>9s}{'図読取':>9s}{'差':>8s}")
    worst = 0.0
    for hc, row in read.items():
        for T, fig_val in row.items():
            calc = carbon_limit_oc(hc, T, 1.5e6)
            worst = max(worst, abs(calc - fig_val))
            print(f"{hc:5.1f}{T:7d}{calc:9.3f}{fig_val:9.3f}{calc - fig_val:+8.3f}")
    print(f"\n最大差 {worst:.3f}（図の線幅だけで約 0.03 相当。読み取り誤差の範囲）")

    # --- CSV ---
    csv_path = os.path.join(out, "carbon_limit.csv")
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write("P_MPa_abs,T_celsius,H_over_C,O_over_C_at_aC1\n")
        for P, tag in ((1.5e6, "1.5"), (P_op, f"{P_op/1e6:.4f}")):
            for T in TEMPS:
                oc, hc = limit_curve(T, P)
                for x, y in zip(oc, hc):
                    fh.write(f"{tag},{T},{y:.4f},{x:.6f}\n")
    print(f"\nCSV 出力: {csv_path}")


if __name__ == "__main__":
    main()
