r"""検討: 千代田の特許の Ac と本実装の差を「平衡定数 K だけの違い」として説明できるか。

`example_chiyoda_patent_check.py`（§6 の 21 点照合）の続き。WO2012140994A1 表5 の
実施例 1〜4（全部 Ac = 0.90 ちょうど）に対して本実装は 2〜4% 高い。これが K の違い
だけで説明できるかを、モデル選択の形で確かめる。

────────────────────────────────────────────────────────────────────────
考え方
────────────────────────────────────────────────────────────────────────
a_C = K_B(T)·Q で、気相の組成（Q = p_CO²/p_CO2）が両者で同じだとすると、比の対数は
**温度だけの関数**になる:

    ln r = ln(a_本実装 / a_特許) = ln K_本実装 − ln K_特許 = −ΔΔG(T) / RT

これを各点で計算し、次の候補を当てはめる（ΔΔX = X_本実装 − X_特許）:

    M0  ln r = B                 定数倍。|B| ≈ 0.013 なら atm と bar の違い
    M1  ln r = A/T               エンタルピーだけのずれ
    M2  ln r = A/T + B           エンタルピー＋エントロピー
    M3  ln r = A/T + B + C·ln T  熱容量の差（相手が log K = a + b/T 型の簡略相関）
    M4  M2 + E·P/T               圧力依存。K ではなくフガシティ（非理想性）の扱いの差

係数の読み方: ln K = −ΔH/RT + ΔS/R (+ ΔCp/R·ln T) なので
    ΔΔH = −R·A,  ΔΔS = +R·B,  ΔΔCp = +R·C
（⚠️ ΔΔS は **+R·B**。−R·B ではない。ln r の定数項は +ΔΔS/R のため）

────────────────────────────────────────────────────────────────────────
結論（main() がそのまま出す）
────────────────────────────────────────────────────────────────────────
**「K だけが違う」という仮定そのものが棄却される。M0〜M4 のどれを選ぶか以前の問題。**

1. **同じ (T, P) で ln r がばらつく。** K だけの差なら、M4（圧力依存）まで含めて
   ln r は (T, P) だけの関数なので、同じ (T, P) の点は同じ値になるはず。ところが
   850 ℃・2.0 MPaG の 8 点（群 A 7 点 + WO 実施例2）で ln r は −0.090〜+0.036 と
   **0.126 も開き**、各点の丸め由来の 1σ（0.001〜0.014）の 9 倍以上。
   この「純誤差」だけで、K(T, P) の形をしたどんなモデルにも越えられない残差の床
   （RMS 0.025・重み付き χ²/dof ≈ 53）ができる。M0〜M4 の RMS（0.081〜0.088）は
   床の 3 倍以上で、AIC/BIC の差も 3 以内＝**どのモデルも M0 から有意に良くならない**。
2. **WO 実施例2 と群 A は同じ 850 ℃・2.0 MPaG なのに比の符号が逆**（+0.036 対
   −0.03〜−0.09）。同じ会社・同じ定義・同じ出典（Stull 1969）。文献ごとに K を
   変えたとしても群 A の内部（7 点で 0.058 の幅、1σ は 0.001〜0.014）は説明できない。
   群 A だけに M2 を当てると RMS 0.017 まで下がるが、温度が 650/850 ℃ の 2 水準しか
   無いので 2 パラメータで 2 つのセル平均をなぞっただけ（検定になっていない）。
   しかも係数は ΔΔH = −7.5 kJ/mol・ΔΔS = −7.2 J/(mol·K) で、まともな出典どうしの差
   （≤ 0.8 kJ/mol、§6 検算 3）の 10 倍。875 ℃ に外挿すると ln r = −0.09 で、
   群 B の実測 +0.07〜+0.15 と符号まで逆になる。
3. **WO 実施例 1〜4 だけなら M0 で収まる**（ln r = +0.018〜+0.036、平均 +0.030
   ＝本実装が 3.0% 高い）。ただし (a) 符号が atm/bar（本実装が 1.3% **低く**なる向き）と
   逆、(b) 0.90 は設計目標値で、供給比の 2 桁丸め（S/C 0.85 など）だけで ±1〜2% 動く。
   K の違いを示す証拠としては弱い。
4. **入口の CH4 分解（有識者 2 点）の ΔΔG_M ≈ −3.2 kJ/mol は、出口 Boudouard の
   ΔΔG_B ≈ −0.3 kJ/mol（M0）と桁が違う。** 3 反応に共通する炭素の項が原因なら両者は
   一致するはずなので、共通原因ではない（ただし出典が別人なので比較自体が弱い）。

→ §6 の結論（残差の主因は特許側の組成の前提・丸め・社内シミュレータ）を補強する。
  K の補正は入れない（§7 の決定どおり）。

実行: PYTHONPATH=. python3 -u examples/example_chiyoda_K_model_fit.py
      （図: references/chiyoda_K_model_fit.png）
要 Cantera: pip install chemflow2[gibbs]
"""

import math
import os
import sys

import matplotlib
import numpy as np
from scipy.stats import spearmanr

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from example_chiyoda_patent_check import ALL_POINTS, _gas, equilibrium_state  # noqa: E402
from example_inlet_carbon_limit import CASES as INLET_CASES  # noqa: E402
from example_inlet_carbon_limit import inlet_partials  # noqa: E402

from chemflow2.core.carbon_activity import R_GAS, carbon_activity  # noqa: E402

#: WO 表5 実施例4。S/C が「−」（水蒸気なし）なので §6 の 21 点からは外していた。
#: CH4:CO2 = 1:2.4 と読むと出口 H2/CO = 0.53（表の 0.5）・CH4 2.3 dry%（表の 2.8）で
#: 表の他の欄とも辻褄が合うので、ここでは補助点として加える。
WO_EX4 = ("C", "WO 表5 実施例4", dict(CH4=1, CO2=2.4), 850, 1.30, 0.90)
POINTS = ALL_POINTS + [WO_EX4]

#: 供給比の丸め幅（最終桁の半分）。群 A はモル比が整数・0.5 刻みなので厳密とみなす。
#: 群 B は 2 桁表記・群 C は表5 の表記どおり。CH4 は基準（= 1）なので振らない。
FEED_HALF = {
    "JP2012 実施例1": dict(H2O=0.005, CO2=0.005),
    "JP2012 実施例2": dict(H2O=0.005, CO2=0.005),
    "JP2012 比較例3": dict(H2O=0.005, CO2=0.005),
    "JP2012 比較例4": dict(H2O=0.005),
    "JP5738743 実施例1": dict(H2O=0.005, CO2=0.005),
    "JP5738743 実施例2": dict(H2O=0.005, CO2=0.005),
    "WO 表5 実施例1": dict(H2O=0.05, CO2=0.05),
    "WO 表5 実施例2": dict(H2O=0.005, CO2=0.005),
    "WO 表5 実施例3": dict(H2O=0.05, CO2=0.05),
    "WO 表5 実施例4": dict(CO2=0.05),
}

MODELS = ("M0", "M1", "M2", "M3", "M4")
MODEL_TEXT = {"M0": "B", "M1": "A/T", "M2": "A/T + B", "M3": "A/T + B + C ln T",
              "M4": "A/T + B + E P/T"}

# --- 配色: example_carbon_limit.py / example_carbon_activity_vs_T.py と揃える ---
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
GROUP_COLOR = {"A": "#2a78d6", "B": "#c9a227", "C": "#eb6834"}
GROUP_LABEL = {"A": "A: JP2004 (pure CH$_4$, lab)", "B": "B: JP2012 / JP5738743 (natural gas)",
               "C": "C: WO2012140994 Table 5 (design Ac = 0.90)"}
MODEL_COLOR = {"M0": "#52514e", "M2": "#7fb069", "M3": "#8b2f8b"}
OUT_PNG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "references",
                       "chiyoda_K_model_fit.png")


def ac(feed, T, mpag, gas):
    p, X = equilibrium_state(feed, T, mpag, gas)
    return carbon_activity(p, T, per_route=True)["Boudouard"], p, X


def build_rows():
    """各点の ln r と、その 1σ（Ac の丸め＋供給比の丸め。いずれも一様分布→半幅/√3）。"""
    gas = _gas()
    rows = []
    for grp, name, feed, T, mp, patent in POINTS:
        a, p, X = ac(feed, T, mp, gas)
        var = (0.005 / patent / math.sqrt(3)) ** 2
        for sp, h in FEED_HALF.get(name, {}).items():
            d = [math.log(ac({**feed, sp: feed.get(sp, 0) + s * h}, T, mp, gas)[0])
                 for s in (+1, -1)]
            var += ((d[0] - d[1]) / 2 / math.sqrt(3)) ** 2
        rows.append(dict(group=grp, name=name, T=T + 273.15, Tc=T, mpag=mp,
                         P_bar=(mp * 1e6 + 101325) / 1e5, patent=patent, ours=a,
                         lnr=math.log(a / patent), sig=math.sqrt(var),
                         co_co2=p["CO"] / p["CO2"], y_h2o=X["H2O"], y_ch4=X["CH4"],
                         supplementary=(name == WO_EX4[1])))
    return rows


def design(model, T, P):
    one = np.ones_like(T)
    return {"M0": np.c_[one], "M1": np.c_[1 / T], "M2": np.c_[1 / T, one],
            "M3": np.c_[1 / T, one, np.log(T)], "M4": np.c_[1 / T, one, P / T]}[model]


def fit(rows, model, weighted=False):
    T = np.array([r["T"] for r in rows])
    P = np.array([r["P_bar"] for r in rows])
    y = np.array([r["lnr"] for r in rows])
    w = 1 / np.array([r["sig"] for r in rows]) if weighted else np.ones_like(y)
    X = design(model, T, P)
    c, *_ = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)
    res = y - X @ c
    n, k = len(y), X.shape[1]
    rss = float(np.sum(res**2))
    chi2 = float(np.sum((res * (1 / np.array([r["sig"] for r in rows]))) ** 2))
    aic = n * math.log(rss / n) + 2 * k
    s2 = rss / (n - k)
    se = np.sqrt(np.diag(s2 * np.linalg.pinv(X.T @ X))) if not weighted else None
    return dict(c=c, se=se, res=res, rss=rss, rms=math.sqrt(rss / n), chi2=chi2,
                chi2_dof=chi2 / (n - k), aic=aic,
                aicc=aic + 2 * k * (k + 1) / (n - k - 1) if n - k - 1 > 0 else math.nan,
                bic=n * math.log(rss / n) + k * math.log(n), k=k, n=n,
                rank=int(np.linalg.matrix_rank(X)))


def cells(rows):
    """(T, P) が同じ点の集まり。K(T, P) 型のモデルはこの中で同じ値しか出せない。"""
    out = {}
    for r in rows:
        out.setdefault((r["Tc"], r["mpag"]), []).append(r)
    return out


def pure_error(rows):
    """各セルの平均を当てた飽和モデル＝K(T, P) 型モデルが到達できる最良の残差。"""
    rss, chi2, dev = 0.0, 0.0, []
    for _key, sub in cells(rows).items():
        y = np.array([r["lnr"] for r in sub])
        w = 1 / np.array([r["sig"] for r in sub]) ** 2
        m, mw = y.mean(), np.sum(w * y) / np.sum(w)
        rss += float(np.sum((y - m) ** 2))
        chi2 += float(np.sum(w * (y - mw) ** 2))
        dev += [(r, r["lnr"] - m) for r in sub if len(sub) > 1]
    n, k = len(rows), len(cells(rows))
    return dict(rss=rss, rms=math.sqrt(rss / n), chi2=chi2, dof=n - k, dev=dev)


def print_fits(rows, label):
    print(f"\n--- {label}（n = {len(rows)}、温度 {sorted({r['Tc'] for r in rows})} ℃）---")
    print(f"  {'':4s}{'形':18s}{'RMS':>8s}{'AIC':>9s}{'AICc':>9s}{'BIC':>9s}"
          f"{'χ²/dof':>9s}   係数（重みなし）")
    for m in MODELS:
        f = fit(rows, m)
        fw = fit(rows, m, weighted=True)
        if f["rank"] < f["k"]:
            print(f"  {m:4s}{MODEL_TEXT[m]:18s}  識別不能（rank {f['rank']} < {f['k']}）")
            continue
        coef = "  ".join(f"{v:+.4g}±{s:.2g}" for v, s in zip(f["c"], f["se"]))
        print(f"  {m:4s}{MODEL_TEXT[m]:18s}{f['rms']:8.4f}{f['aic']:9.2f}{f['aicc']:9.2f}"
              f"{f['bic']:9.2f}{fw['chi2_dof']:9.1f}   {coef}")
    pe = pure_error(rows)
    print(f"  {'床':4s}{'セル平均（飽和）':16s}{pe['rms']:8.4f}"
          f"{'':27s}{pe['chi2'] / max(pe['dof'], 1):9.1f}   ← K(T,P) 型モデルの到達限界"
          f"（{len(cells(rows))} セル）")
    return pe


def main():
    rows = build_rows()
    main21 = [r for r in rows if not r["supplementary"]]

    print("=== 「K だけが違う」と仮定したときの ln r = ln(a_本実装 / a_特許) ===\n")
    print(f"{'群':3s}{'点':20s}{'T℃':>5s}{'MPaG':>6s}{'特許':>6s}{'本実装':>8s}"
          f"{'ln r':>8s}{'1σ':>7s}{'pCO/pCO2':>10s}")
    for r in rows:
        mark = "  ※補助点" if r["supplementary"] else ""
        print(f"{r['group']:3s}{r['name']:20s}{r['Tc']:5.0f}{r['mpag']:6.2f}"
              f"{r['patent']:6.2f}{r['ours']:8.4f}{r['lnr']:+8.4f}{r['sig']:7.4f}"
              f"{r['co_co2']:10.3f}{mark}")

    print("\n--- WO2012140994A1 表5 実施例 1〜4（全部 Ac = 0.90）---")
    wo = [r for r in rows if r["group"] == "C"]
    for r in wo:
        print(f"  {r['name']:16s} 本実装 {r['ours']:.4f} → 本実装が "
              f"{(r['ours'] / r['patent'] - 1) * 100:+.1f}%   ln r {r['lnr']:+.4f} ± {r['sig']:.4f}")
    f0 = fit(wo, "M0")
    print(f"  M0 当てはめ: B = {f0['c'][0]:+.4f} ± {f0['se'][0]:.4f}"
          f"（本実装が {math.expm1(f0['c'][0]) * 100:+.1f}%）")
    print(f"  ・atm/bar の違いなら |B| = ln 1.01325 = {math.log(1.01325):.4f}。しかも特許が bar なら")
    print("    a_特許 = 1.013·a_本実装 で**本実装が低く**出る向き。観測は逆符号なので単位では説明できない。")
    print("  ・0.90 は設計目標値。供給比の 2 桁丸めだけで上表の 1σ（0.3〜1.4%）動くので、")
    print("    4 点が M0 に乗ること自体は K の違いの証拠として弱い。")

    pe21 = print_fits(main21, "全 21 点（§6 と同じ点）")
    print_fits(rows, "22 点（WO 実施例4 を加える）")
    print_fits([r for r in main21 if r["group"] == "A"], "群 A のみ（最も信頼できる 12 点）")

    print("\n--- lack-of-fit: 同じ (T, P) の点が同じ ln r になっているか ---")
    print(f"  {'T℃':>5s}{'MPaG':>6s}{'n':>3s}{'ln r の範囲':>22s}{'幅':>8s}{'最小 1σ':>9s}  内訳")
    for (Tc, mp), sub in sorted(cells(main21).items()):
        y = [r["lnr"] for r in sub]
        grp = "".join(sorted({r["group"] for r in sub}))
        print(f"  {Tc:5.0f}{mp:6.2f}{len(sub):3d}{min(y):+11.4f} 〜{max(y):+8.4f}"
              f"{max(y) - min(y):8.4f}{min(r['sig'] for r in sub):9.4f}  群 {grp}")
    print(f"  → 純誤差の RMS {pe21['rms']:.4f}、重み付き χ²/dof = "
          f"{pe21['chi2'] / pe21['dof']:.0f}（丸めだけなら ≈ 1）。")
    print("     **K(T, P) の形をしたモデルは、どれを選んでもこの床より下には行けない。**")
    print("     M0〜M4 の RMS が床に届かず AIC/BIC もほぼ横並び＝温度で説明できる成分が無い。")

    print("\n--- 残差（同じセル内の偏差）が温度以外の量で説明できるか（Spearman）---")
    dev = pe21["dev"]
    d = np.array([x[1] for x in dev])
    for key, lab, tf in (("co_co2", "ln(pCO/pCO2)", math.log), ("y_h2o", "y_H2O", float),
                         ("y_ch4", "y_CH4", float), ("patent", "ln Ac(特許)", math.log)):
        v = np.array([tf(x[0][key]) for x in dev])
        rho, p = spearmanr(v, d)
        print(f"  {lab:14s} ρ = {rho:+.2f}  (p = {p:.2f}, n = {len(d)})")
    print("  → 単一の組成変数ではきれいに説明できない。特許側の組成が完全平衡でない（§6 検算 4）など、")
    print("     点ごとに違う前提が効いていると読むのが自然。")

    print("\n--- 入口の CH4 分解（有識者 2 点）と比べる: 炭素の項は共通か ---")
    fB = fit(main21, "M0")
    fB2 = fit(main21, "M2")
    for name, (feed, P, T_exp) in INLET_CASES.items():
        a_ours = carbon_activity(inlet_partials(feed, P), T_exp)   # 有識者は a = 1
        Tk = T_exp + 273.15
        lnr_m = math.log(a_ours)
        dGM = -R_GAS * Tk * lnr_m
        dGB0 = -R_GAS * Tk * fB["c"][0]
        dGB2 = -R_GAS * Tk * (fB2["c"][0] / Tk + fB2["c"][1])
        print(f"  {name:10s} {T_exp:.0f}℃: ln r_M = {lnr_m:+.3f} → ΔΔG_M = {dGM / 1000:+.2f} kJ/mol"
              f"   対 Boudouard ΔΔG_B（M0）{dGB0 / 1000:+.2f} /（M2 外挿）{dGB2 / 1000:+.2f}")
    print("  → 桁が違う。3 反応に共通する炭素の項（黒鉛の g）が原因なら一致するはずなので、")
    print("     共通原因ではない。ただし出口は特許・入口は有識者と出典が別で、比較自体が弱い。")

    fA = fit([r for r in main21 if r["group"] == "A"], "M2")
    A_, B_ = fA["c"]
    print("\n--- 群 A だけなら M2 が良く見える件 ---")
    print(f"  M2: A = {A_:+.1f} K, B = {B_:+.4f} → ΔΔH = −R·A = {-R_GAS * A_ / 1000:+.2f} kJ/mol、"
          f"ΔΔS = +R·B = {R_GAS * B_:+.2f} J/(mol·K)")
    print("  ・温度が 650/850 ℃ の 2 水準だけなので、2 パラメータで 2 つのセル平均をなぞっただけ。")
    print("  ・|ΔΔH| 7.5 kJ/mol は、まともな出典どうしの差（≤ 0.8 kJ/mol、§6 検算 3）の 10 倍。")
    lnr875 = A_ / 1148.15 + B_
    obs = [r["lnr"] for r in main21 if r["group"] == "B"]
    print(f"  ・875 ℃ に外挿すると ln r = {lnr875:+.3f}。群 B の実測 {min(obs):+.3f}〜{max(obs):+.3f}"
          " と符号まで逆。")

    print("\n=== 結論 ===")
    print("  ・「K だけが違う」なら同じ (T, P) で ln r は同じになるはずだが、850 ℃・2.0 MPaG の")
    print("    8 点で 0.126 開く（各点の 1σ の 9 倍以上）。**仮定そのものが棄却される**。")
    print("  ・したがって M0〜M4 の選択は意味を持たない（どれも M0 から有意に良くならない）。")
    print("  ・WO 実施例 1〜4 だけなら本実装が一律 +3.0%（M0）。ただし atm/bar とは逆符号で、")
    print("    設計目標値の供給比丸めでも 1〜2% 動くので、K の版差を示す証拠にはならない。")
    print("  ・§6/§7 の結論（K は gri30/JANAF のまま・補正しない）を変える材料は無い。")

    plot(rows, main21)


def plot(rows, main21):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 5.6), dpi=130,
                                   gridspec_kw=dict(width_ratios=[1.35, 1]))
    fig.patch.set_facecolor(SURFACE)
    for ax in (ax1, ax2):
        ax.set_facecolor(SURFACE)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelcolor=INK2)
        ax.yaxis.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
    ax1.axhline(0, color=MUTED, lw=1.0)

    # 左: ln r 対 1000/T と、21 点に当てはめたモデル
    rng = np.random.default_rng(0)
    for g in ("A", "B", "C"):
        sub = [r for r in rows if r["group"] == g]
        x = np.array([1000 / r["T"] for r in sub]) + rng.uniform(-0.004, 0.004, len(sub))
        ax1.errorbar(x, [r["lnr"] for r in sub], yerr=[r["sig"] for r in sub], fmt="none",
                     ecolor=GROUP_COLOR[g], elinewidth=1.0, alpha=0.7)
        face = [SURFACE if r["supplementary"] else GROUP_COLOR[g] for r in sub]
        ax1.scatter(x, [r["lnr"] for r in sub], s=42, c=face, edgecolors=GROUP_COLOR[g],
                    linewidths=1.6, zorder=3, label=GROUP_LABEL[g])
    Tg = np.linspace(800, 1170, 200)
    for m, ls in (("M0", "-"), ("M2", "--"), ("M3", ":")):
        f = fit(main21, m)
        y = design(m, Tg, np.full_like(Tg, 21.0)) @ f["c"]
        ax1.plot(1000 / Tg, y, color=MODEL_COLOR[m], lw=2.0, ls=ls,
                 label=f"{m}: ln r = {MODEL_TEXT[m]}  (RMS {f['rms']:.3f})")
    ax1.axhline(-math.log(1.01325), color=MUTED, lw=1.0, ls=(0, (2, 2)))
    ax1.text(1.235, -math.log(1.01325) - 0.004, "patent in bar (−0.013)", color=MUTED,
             fontsize=7.5, ha="right", va="top")
    ax1.set_xlabel("1000 / T  [1/K]", color=INK2)
    ax1.set_ylabel("ln r = ln(a$_{C,\\,this\\ work}$ / a$_{C,\\,patent}$)", color=INK2)
    ax1.set_title("ln r vs temperature — no T-only model fits", color=INK, fontsize=10.5,
                  loc="left")
    top = ax1.secondary_xaxis("top", functions=(lambda x: 1000 / x - 273.15,
                                                lambda t: 1000 / (t + 273.15)))
    top.set_xlabel("T [°C]", color=MUTED, fontsize=8.5)
    top.tick_params(colors=MUTED, labelcolor=MUTED, labelsize=8)
    ax1.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper left",
               bbox_to_anchor=(0.0, -0.13), ncol=2)
    ax1.set_ylim(-0.17, 0.2)

    # 右: 850 ℃・2.0 MPaG の 8 点。K(T, P) 型モデルはここで 1 つの値しか出せない
    sub = sorted([r for r in main21 if r["Tc"] == 850 and r["mpag"] == 2.0],
                 key=lambda r: r["lnr"])
    ypos = np.arange(len(sub))
    for i, r in enumerate(sub):
        c = GROUP_COLOR[r["group"]]
        ax2.errorbar(r["lnr"], i, xerr=r["sig"], fmt="o", color=c, ms=7, mec=SURFACE,
                     mew=1.5, elinewidth=1.4, capsize=0)
    m = np.mean([r["lnr"] for r in sub])
    ax2.axvline(m, color=INK2, lw=1.5, ls="--")
    ax2.text(m + 0.004, 0.5, "best any K(T,P) model\ncan do: one value", color=INK2,
             fontsize=7.5, va="center")
    names = [r["name"].replace("JP2004 表1 ", "JP2004 ").replace("WO 表5 ", "WO ")
             .replace("実施例", "Ex.").replace("比較例", "Comp.").replace("改質条件", "reforming")
             .replace("表3 Ni", "Tab.3 Ni") for r in sub]
    ax2.set_yticks(ypos)
    ax2.set_yticklabels(names, fontsize=8)
    ax2.yaxis.grid(False)
    ax2.xaxis.grid(True, color=GRID, lw=0.8)
    ax2.axvline(0, color=MUTED, lw=1.0)
    ax2.set_xlabel("ln r  (error bar = 1σ from rounding)", color=INK2)
    ax2.set_title("Same T and P (850 °C, 2.0 MPaG): spread 0.126", color=INK,
                  fontsize=10.5, loc="left")

    fig.tight_layout()
    fig.savefig(OUT_PNG, facecolor=SURFACE)
    print(f"\nsaved {os.path.relpath(OUT_PNG)}")


if __name__ == "__main__":
    main()
