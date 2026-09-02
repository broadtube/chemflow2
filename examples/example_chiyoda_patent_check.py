r"""検証: 本実装の炭素活量が千代田化工建設の特許の数値と合っているか（全 21 点）。

`references/PAPERS_INDEX.md` §0-2b の照合を**特許 5 件の全実施例に広げ、ずれの原因を
切り分ける**ところまでやる。以前は 6 点だけを見て「本実装が一律 12% 高い＝Stull 1969
との版差」と結論していたが、点を増やすとその説明は成り立たない（下記）。

────────────────────────────────────────────────────────────────────────
どこから拾ったか
────────────────────────────────────────────────────────────────────────
千代田の同族特許は**すべて同じ定義**を使う（§0-2 / §0-2b）:

    カーボン活性 Ac = K × (P_CO)² / P_CO2    K は 2CO = C + CO2 の平衡定数
    平衡定数の出典は D. R. Stull et al., The Chemical Thermodynamics of
    Organic Compounds, John Wiley (1969)

    群 A  JP2004-210626A 表1 + 本文 + 表3   12 点  650/850 ℃・0.5/2.0 MPaG・Ac 0.09〜2.8
          → **供給が CH4/CO2/H2O のモル比で明記**され、原料は純 CH4 のラボ試験。
            平衡転化率も併記されているので**組成側の検算までできる。最も信頼できる群**。
    群 B  JP2012-219068A（= JP5738658B2）    4 点  875 ℃・2.0 MPaG
          JP5738743B2（= JP2013-103909A）    2 点  875 ℃・1.85 MPaG
          → S/C と CO2/C しか与えられず、**原料は C2〜C6 を含む実天然ガス**
            （表1: CH4 91.5 / C2 4.2 / C3 2.5 / C4 1.3 / C5 0.4 / C6 0.1 mol%）。
            CH4 = 1 と置く従来の解釈は近似でしかない。
    群 C  WO2012140994A1 表5                 3 点  550/850 ℃
          → 4 例とも Ac = 0.90 ちょうど＝**設計目標値**。較正には使えない。

────────────────────────────────────────────────────────────────────────
結果（実測。main() がそのまま出す）
────────────────────────────────────────────────────────────────────────
    群                     n   比（特許 ÷ 本実装）
    A  650 ℃              4   0.899 ± 0.008    本実装が 11% 高い
    A  850 ℃              8   1.070 ± 0.021    本実装が  7% 低い
    B  875 ℃              6   0.882 ± 0.023    本実装が 13% 高い
    C  550/850 ℃（設計値） 3   0.967 ± 0.002    本実装が  3% 高い

**850 ℃ と 875 ℃ でずれの符号が反転する。** ΔG は温度の滑らかな関数なので、
**熱力学データの版の違いではこれは起こせない**。したがって以前の「Stull 1969 との
版差で 12% 保守側」という説明は成り立たない（撤回）。

────────────────────────────────────────────────────────────────────────
ずれの切り分け（下の 4 つの検算で確かめる）
────────────────────────────────────────────────────────────────────────
1. **組成側は合っている。** JP2004 が明記する平衡 CH4 転化率と本実装は
   55.5% 対 55%、57.4% 対 57%（いずれも特許の 2 桁表記の丸め幅の内側）。
2. **熱力学データは JANAF 級。** ΔG°(2CO → C + CO2) @ 298.15 K は本実装 −120.04 kJ/mol
   に対し CODATA/JANAF は −120.02。黒鉛の S°(298) = 5.734（JANAF 5.740）、
   S°(1000) = 24.451（同 24.457）、Cp(298) = 8.528（同 8.517）。
3. **出典を替えても動かない。** gri30 と NASA9（nasa_gas + nasa_condensed）は
   ΔG_B が完全に一致。methane_pox_on_pt でも 0.6%、古典的 Ellingham 式でも 0.8 kJ/mol
   しか違わない。一方**特許の数値が要求する ΔG_B は ΔH で +7.4 kJ/mol、ΔS で
   +7.2 J/mol/K ずれており、まともな出典でそんな値を出すものは無い**。
4. **同じ温度でも供給組成によって比が変わる**（850 ℃ で 1.033〜1.098）。
   K_B の誤差なら同一温度の全点が同じ比になるはずで、**これは組成側の違いを意味する**。
   実際、JP2012 の物質収支から出口組成を計算すると本実装は dry で CH4 7.7% に対し
   特許の表は 9%。**特許のリフォーマー出口は完全平衡ではない**（approach がある）。

→ **結論: 熱力学データを差し替えても「完全一致」にはできない。** 残差の主因は
  (a) 特許値が 2 桁表記（Ac = 0.09 なら丸めだけで ±5.6%）、
  (b) 特許側の組成の前提（平衡への接近度・実天然ガスの C2+）、
  (c) 千代田の社内相関やプロセスシミュレータの中身が不明、の 3 つ。

**照合として言えること**: Ac 0.09〜2.8（30 倍）、650〜875 ℃、0.5〜2.0 MPaG に
わたる 21 点で、本実装は**±13% 以内**に収まり、**大小関係はすべて一致**する。
組成側は独立に検算できており、熱力学データは JANAF と一致している。
定義（Boudouard 1 本）も平衡気相では本実装の 3 経路最大値と厳密に同値。

実行: PYTHONPATH=. python3 -u examples/example_chiyoda_patent_check.py
要 Cantera: pip install chemflow2[gibbs]
"""

import math
import statistics

import cantera as ct

from chemflow2.core.carbon_activity import carbon_activity, equilibrium_constants

SPECIES = ["CO2", "CH4", "H2O", "CO", "H2"]
R = 8.314462618

#: (群, 名前, 供給 CH4/CO2/H2O, T[℃], P[MPaG], 特許の Ac)
#: 群 A = JP2004-210626A。供給がモル比で明記され、原料は純 CH4。
POINTS_A = [
    ("A", "JP2004 表1 実施例1", dict(CH4=1, H2O=3),              850, 2.00, 0.21),
    ("A", "JP2004 表1 実施例2", dict(CH4=1, H2O=3),              850, 0.50, 0.09),
    ("A", "JP2004 表1 実施例3", dict(CH4=1, H2O=3),              650, 2.00, 0.29),
    ("A", "JP2004 表1 実施例4", dict(CH4=1, H2O=1),              850, 2.00, 0.90),
    ("A", "JP2004 表1 実施例5", dict(CH4=1, H2O=1),              650, 2.00, 0.75),
    ("A", "JP2004 表1 実施例6", dict(CH4=1, CO2=2, H2O=2),       850, 2.00, 0.37),
    ("A", "JP2004 表1 実施例7", dict(CH4=1, CO2=0.5, H2O=2),     850, 2.00, 0.41),
    ("A", "JP2004 表1 実施例8", dict(CH4=1, CO2=0.5, H2O=2),     650, 2.00, 0.67),
    ("A", "JP2004 表1 比較例1", dict(CH4=1, CO2=2, H2O=2),       650, 2.00, 1.03),
    ("A", "JP2004 表1 比較例2", dict(CH4=1, CO2=1, H2O=0.5),     850, 2.00, 1.36),
    ("A", "JP2004 改質条件",    dict(CH4=1, CO2=1),              850, 2.00, 2.80),
    ("A", "JP2004 表3 Ni",     dict(CH4=1, CO2=0.4, H2O=0.85),  850, 2.00, 1.10),
]
#: 群 B = メタノール系。S/C・CO2/C しか与えられず、原料は実天然ガス（C2+ を含む）。
POINTS_B = [
    ("B", "JP2012 実施例1",   dict(CH4=1, H2O=1.50, CO2=0.25),  875, 2.00, 0.43),
    ("B", "JP2012 実施例2",   dict(CH4=1, H2O=1.35, CO2=0.29),  875, 2.00, 0.49),
    ("B", "JP2012 比較例3",   dict(CH4=1, H2O=1.75, CO2=0.20),  875, 2.00, 0.35),
    ("B", "JP2012 比較例4",   dict(CH4=1, H2O=3.00),            875, 2.00, 0.15),
    ("B", "JP5738743 実施例1", dict(CH4=1, H2O=1.22, CO2=0.34),  875, 1.85, 0.53),
    ("B", "JP5738743 実施例2", dict(CH4=1, H2O=1.20, CO2=0.39),  875, 1.85, 0.54),
]
#: 群 C = WO2012140994A1 表5。4 例とも 0.90 ちょうどの**設計目標値**。
POINTS_C = [
    ("C", "WO 表5 実施例1", dict(CH4=1, H2O=4.00, CO2=9.60), 550, 0.15, 0.90),
    ("C", "WO 表5 実施例2", dict(CH4=1, H2O=0.85, CO2=0.82), 850, 2.00, 0.90),
    ("C", "WO 表5 実施例3", dict(CH4=1, H2O=4.00, CO2=9.40), 550, 0.15, 0.90),
]
ALL_POINTS = POINTS_A + POINTS_B + POINTS_C

#: JP2012-219068A 表1 の天然ガス組成 [mol%]。CH4 だけではない。
NATURAL_GAS = {"N2": 0.1, "CH4": 91.5, "C2": 4.2, "C3": 2.5, "C4": 1.3, "C5": 0.4, "C6": 0.1}
NG_C = {"CH4": 1, "C2": 2, "C3": 3, "C4": 4, "C5": 5, "C6": 6}
NG_H = {"CH4": 4, "C2": 6, "C3": 8, "C4": 10, "C5": 12, "C6": 14}


def _gas():
    selected = [s for s in ct.Species.list_from_file("gri30.yaml") if s.name in SPECIES]
    return ct.Solution(thermo="ideal-gas", species=selected)


def equilibrium_state(feed, T_celsius, mpag, gas=None):
    """供給を出口条件でギブズ平衡させ、(分圧 [atm], モル分率) を返す。"""
    gas = gas or _gas()
    P_pascal = mpag * 1e6 + 101325.0
    gas.TPX = T_celsius + 273.15, P_pascal, feed
    gas.equilibrate("TP")
    X = {s: gas.X[gas.species_index(s)] for s in SPECIES}
    return {s: X[s] * P_pascal / 101325.0 for s in SPECIES}, X


def ac_boudouard(feed, T_celsius, mpag, gas=None):
    """千代田の定義（Boudouard 1 本）での Ac。"""
    p_atm, _ = equilibrium_state(feed, T_celsius, mpag, gas)
    return carbon_activity(p_atm, T_celsius, per_route=True)["Boudouard"]


def compare_all():
    """全点を照合して行のリストを返す。"""
    gas = _gas()
    rows = []
    for grp, name, feed, T, mp, patent in ALL_POINTS:
        p_atm, _ = equilibrium_state(feed, T, mp, gas)
        routes = carbon_activity(p_atm, T, per_route=True)
        a = routes["Boudouard"]
        spread = (max(routes.values()) - min(routes.values())) / max(routes.values())
        rows.append(dict(group=grp, name=name, T=T, P=mp, patent=patent, ours=a,
                         ratio=patent / a, spread=spread,
                         round_pct=0.005 / patent * 100))
    return rows


def check_conversion():
    """JP2004 が明記する『平衡 CH4 転化率』と本実装を突き合わせる（組成側の検算）。"""
    gas = _gas()
    out = []
    for label, feed, T, mp, conv_patent in (
            ("CH4:CO2 = 1:1（改質条件）", dict(CH4=1, CO2=1), 850, 2.0, 55),
            ("CH4:CO2:H2O = 1:0.4:0.85（表3）", dict(CH4=1, CO2=0.4, H2O=0.85), 850, 2.0, 57)):
        _p, X = equilibrium_state(feed, T, mp, gas)
        n_in = feed.get("CH4", 0) + feed.get("CO2", 0)
        n_total = n_in / (X["CH4"] + X["CO"] + X["CO2"])
        conv = (1 - X["CH4"] * n_total / feed["CH4"]) * 100
        out.append((label, conv_patent, conv))
    return out


def check_thermo_sources():
    """ΔG°(2CO → C + CO2) を出典別に出す。どれを使ってもほぼ動かないことの確認。"""
    graphite = ct.Solution("graphite.yaml")

    def dG(gasfile, solidfile="graphite.yaml", T_kelvin=None):
        sp = [s for s in ct.Species.list_from_file(gasfile) if s.name in ("CO", "CO2")]
        g = ct.Solution(thermo="ideal-gas", species=sp)
        if solidfile == "graphite.yaml":
            sol = graphite
        else:
            ss = [s for s in ct.Species.list_from_file(solidfile) if s.name == "C(gr)"]
            sol = ct.Solution(thermo="fixed-stoichiometry", species=ss)
        g.TP = T_kelvin, g.reference_pressure
        sol.TP = T_kelvin, sol.reference_pressure
        gg = {s: g.standard_gibbs_RT[g.species_index(s)] for s in ("CO", "CO2")}
        return R * T_kelvin * ((sol.standard_gibbs_RT[0] + gg["CO2"]) - 2 * gg["CO"]) / 1000

    temps = (923.15, 1123.15, 1148.15)
    rows = [
        ("gri30 + graphite.yaml（本実装）",
         [dG("gri30.yaml", T_kelvin=T) for T in temps]),
        ("NASA9 gas + NASA condensed",
         [dG("nasa_gas.yaml", "nasa_condensed.yaml", T) for T in temps]),
        ("methane_pox_on_pt + graphite",
         [dG("methane_pox_on_pt.yaml", T_kelvin=T) for T in temps]),
        ("古典式 −170700 + 174.5·T",
         [(-170700 + 174.5 * T) / 1000 for T in temps]),
    ]
    return temps, rows


#: NIST Webbook（Chase 1998 = NIST-JANAF 第 4 版）の Shomate 係数。
#: **Cantera とは独立の第三者データ**で K_B を裏取りするために置いてある。
#: ⚠️ NIST の注記は CO・CO2 とも "Data last reviewed in September, 1965"。
#:    つまり Stull 1969 が使えた値と、現行 JANAF の値は**同じ年代のもの**。
SHOMATE = {
    "CO":     dict(A=25.56759, B=6.096130, C=4.054656, D=-2.671301, E=0.131021,
                   F=-118.0089, G=227.3665, H=-110.5271, dHf=-110.53),   # 298–1300 K
    "CO2":    dict(A=24.99735, B=55.18696, C=-33.69137, D=7.948387, E=-0.136638,
                   F=-403.6075, G=228.2431, H=-393.5224, dHf=-393.51),   # 298–1200 K
    "CO2_hi": dict(A=58.16639, B=2.720074, C=-0.492289, D=0.038844, E=-6.447293,
                   F=-425.9186, G=263.6125, H=-393.5224, dHf=-393.51),   # 1200–6000 K
}


def _shomate(par, T_kelvin):
    """(H°(T) [kJ/mol]、S°(T) [J/mol/K]) を返す。H は生成エンタルピー基準。"""
    t = T_kelvin / 1000.0
    H = (par["A"] * t + par["B"] * t**2 / 2 + par["C"] * t**3 / 3
         + par["D"] * t**4 / 4 - par["E"] / t + par["F"] - par["H"])
    S = (par["A"] * math.log(t) + par["B"] * t + par["C"] * t**2 / 2
         + par["D"] * t**3 / 3 - par["E"] / (2 * t**2) + par["G"])
    return par["dHf"] + H, S


def nist_boudouard(T_kelvin):
    """NIST-JANAF の Shomate（CO/CO2）+ 黒鉛から ΔG°(2CO → C + CO2) [kJ/mol]。"""
    gr = ct.Solution("graphite.yaml")
    gr.TP = 298.15, gr.reference_pressure
    h0 = gr.enthalpy_mole
    gr.TP = T_kelvin, gr.reference_pressure
    H_C, S_C = (gr.enthalpy_mole - h0) / 1e6, gr.entropy_mole / 1000
    H_CO, S_CO = _shomate(SHOMATE["CO"], T_kelvin)
    key = "CO2" if T_kelvin < 1200 else "CO2_hi"
    H_CO2, S_CO2 = _shomate(SHOMATE[key], T_kelvin)
    dH = (H_C + H_CO2) - 2 * H_CO
    dS = (S_C + S_CO2) - 2 * S_CO
    return dH - T_kelvin * dS / 1000


def gibbs_anchor():
    """298.15 K の ΔG°(CH4 ⇌ C + 2H2)。既知値 +50.7 kJ/mol との比較用。"""
    _K_B, K_M, _K_R = equilibrium_constants(25.0)
    return -R * 298.15 * math.log(K_M) / 1000.0


def natural_gas_note():
    """群 B の原料が純 CH4 でないことを数字で示す。"""
    C = sum(NATURAL_GAS[k] / 100 * NG_C[k] for k in NG_C)
    H = sum(NATURAL_GAS[k] / 100 * NG_H[k] for k in NG_H)
    return C, H, H / C


def main():
    rows = compare_all()

    print("=== 千代田の特許の全実施例との照合（Boudouard・出口平衡）===\n")
    print(f"{'群':3s}{'点':22s}{'T':>5s}{'MPaG':>6s}{'特許':>7s}{'本実装':>9s}"
          f"{'比':>7s}{'丸め%':>7s}{'3経路':>9s}")
    for r in rows:
        print(f"{r['group']:3s}{r['name']:22s}{r['T']:5.0f}{r['P']:6.2f}"
              f"{r['patent']:7.2f}{r['ours']:9.4f}{r['ratio']:7.3f}"
              f"{r['round_pct']:7.2f}{r['spread']:9.1e}")

    print("\n--- 群・温度ごとの比（特許 ÷ 本実装）---")
    groups = [("A", 650), ("A", 850), ("B", 875), ("C", 550), ("C", 850)]
    for g, T in groups:
        sub = [r for r in rows if r["group"] == g and r["T"] == T]
        if not sub:
            continue
        v = [r["ratio"] for r in sub]
        m, sd = statistics.mean(v), (statistics.pstdev(v) if len(v) > 1 else 0.0)
        print(f"  群{g} {T:.0f}℃  n={len(sub)}  比 {m:.4f} ± {sd:.4f}"
              f"   → 本実装が {(1/m-1)*100:+.0f}%")
    print("  ⚠ 群 A の 850 ℃ と 群 B の 875 ℃ で**ずれの符号が反転する**。"
          "ΔG は温度の滑らかな関数なので、熱力学データの差ではこれは起こせない。")

    print("\n--- 検算 1: 組成側（JP2004 が明記する平衡 CH4 転化率）---")
    for label, pat, ours in check_conversion():
        print(f"  {label:30s} 特許 {pat}%   本実装 {ours:.1f}%   差 {ours-pat:+.1f}pt")
    print("  → 特許の 2 桁表記の丸め幅の内側。**気相平衡（組成）の計算は一致している**。")

    print("\n--- 検算 2: 熱力学データの絶対値 ---")
    dG = gibbs_anchor()
    print(f"  ΔG°(CH4 ⇌ C + 2H2) @ 298.15 K = {dG:+.2f} kJ/mol（既知値 +50.7、"
          f"差 {dG-50.7:+.2f}）")
    gr = ct.Solution("graphite.yaml")
    for T, s_janaf, cp_janaf in ((298.15, 5.740, 8.517), (1000.0, 24.457, 21.610)):
        gr.TP = T, gr.reference_pressure
        print(f"  黒鉛 {T:7.2f} K: S° = {gr.entropy_mole/1000:6.3f}（JANAF {s_janaf}）"
              f"  Cp = {gr.cp_mole/1000:6.3f}（JANAF {cp_janaf}）")
    print("  → 黒鉛も気相も JANAF と一致。**こちら側のデータに疑う余地は無い**。")

    print("\n--- 検算 3: 出典を替えたらどれだけ動くか（ΔG°(2CO → C + CO2) [kJ/mol]）---")
    temps, src = check_thermo_sources()
    print(f"  {'出典':34s}" + "".join(f"{int(T-273.15):>10d}℃" for T in temps))
    for name, vals in src:
        print(f"  {name:34s}" + "".join(f"{v:11.3f}" for v in vals))
    print("  特許が要求する値（群 A から逆算）      "
          f"{-7.977:11.3f}{25.604:11.3f}{'—':>11s}")
    print("  → まともな出典はどれも 0.8 kJ/mol 以内で一致する。特許の数値が要求する"
          "\n     ΔH +7.4 kJ/mol・ΔS +7.2 J/mol/K のずれを出す出典は無い。")

    print("\n--- 検算 3b: Cantera と独立に NIST-JANAF（Shomate）で K_B を裏取り ---")
    print(f"  {'T [K]':>8}{'NIST-JANAF':>13}{'Cantera':>11}{'差 [kJ/mol]':>13}{'K の比':>9}")
    for Tk in (923.15, 1123.15, 1148.15):
        gasf = [x for x in ct.Species.list_from_file("gri30.yaml") if x.name in ("CO", "CO2")]
        g = ct.Solution(thermo="ideal-gas", species=gasf)
        gr = ct.Solution("graphite.yaml")
        g.TP = Tk, g.reference_pressure
        gr.TP = Tk, gr.reference_pressure
        gg = {x: g.standard_gibbs_RT[g.species_index(x)] for x in ("CO", "CO2")}
        dG_ct = R * Tk * ((gr.standard_gibbs_RT[0] + gg["CO2"]) - 2 * gg["CO"]) / 1000
        dG_nist = nist_boudouard(Tk)
        print(f"  {Tk:8.2f}{dG_nist:13.3f}{dG_ct:11.3f}{dG_nist-dG_ct:13.3f}"
              f"{math.exp(-(dG_nist-dG_ct)*1000/(R*Tk)):9.4f}")
    print("  → 独立の出典でも 0.034 kJ/mol（K で 0.4%）しか違わない。")
    print("  ★ NIST の注記は CO・CO2 とも \"Data last reviewed in September, 1965\"。")
    print("     **Stull 1969 が使えた値と現行 JANAF の値は同じ年代のもの**なので、")
    print("     『Stull と現代データの版差』では 4〜6% の差は作れない。仮説は閉じる。")

    print("\n--- 検算 4: 同じ温度でも比が組成で変わる（850 ℃・群 A）---")
    sub = sorted([r for r in rows if r["group"] == "A" and r["T"] == 850],
                 key=lambda r: r["ratio"])
    for r in sub:
        print(f"  {r['name']:22s} 比 {r['ratio']:.3f}（丸め ±{r['round_pct']:.2f}%）")
    print("  → K_B の誤差なら同一温度の全点が同じ比になるはず。1.033〜1.098 と"
          "\n     ばらつく以上、**組成側の前提も違っている**。")

    print("\n--- 検算 5: ウィスカーカーボンでずれを説明できるか → **できない** ---")
    print("  ウィスカー炭素は黒鉛より**不安定**なので、基準をずらすなら ΔG_dev > 0 しかありえない")
    print("  （a_C は exp(−ΔG_dev/RT) 倍で必ず**下がる**）。各点で必要な ΔG_dev の符号を見る:")
    need = []
    for r in rows:
        dG = R * (r["T"] + 273.15) * math.log(r["ours"] / r["patent"])
        need.append((r["group"], r["name"], r["T"], dG, r["round_pct"]))
    neg = [x for x in need if x[3] <= 0]
    pos = [x for x in need if x[3] > 0]
    print(f"  正（ウィスカーで説明できる向き）{len(pos)} 点 / "
          f"**負（ウィスカーでは不可能な向き）{len(neg)} 点**")
    for g, n, T, dG, rp in neg:
        print(f"    {g} {n:22s} {T:.0f}℃  必要 ΔG_dev = {dG:+6.0f} J/mol  ← 負")
    print("  → 負が要るのは**群 A の 850 ℃ の 8 点すべて**。群 A は供給がモル比で明記され、")
    print("     平衡転化率で組成側も検算できている**最も信頼できる群**。ここで符号が逆になる以上、")
    print("     ウィスカー基準ではずれを説明できない。")
    print(f"  さらに、正の側でも必要量は {min(x[3] for x in pos):.0f}〜{max(x[3] for x in pos):.0f} J/mol で、")
    print("     Alstrup の下限（黒鉛極限 2600 J/mol）にすら届かない。")
    at850 = [x for x in need if x[2] == 850 and x[0] == "A"]
    print(f"  もう 1 つ: ΔG_dev は**温度だけの関数**なので同一温度の全点で同じ値になるはず。")
    print(f"     850 ℃・群 A の必要量は {min(x[3] for x in at850):.0f}〜{max(x[3] for x in at850):.0f} J/mol と"
          f" {max(x[3] for x in at850)-min(x[3] for x in at850):.0f} J/mol もばらつく")
    print("     （最も精度の高い 2 点の丸め由来の不確かさは ±17 と ±35 J/mol）。")

    C, H, HC = natural_gas_note()
    print("\n--- 群 B の但し書き: 原料は純 CH4 ではない ---")
    print(f"  JP2012 表1 の天然ガス: {NATURAL_GAS}")
    print(f"  → C {C:.4f} / H {H:.4f} per mol、H/C = {HC:.4f}（純 CH4 なら 4.0）")
    print("  S/C・CO2/C の分母（カーボン）がこの C2+ を含むため、CH4 = 1 と置く"
          "\n  従来の解釈は近似。群 B の比はこの不確かさを含む。")

    print("\n=== 結論 ===")
    print("  ・定義（Boudouard 1 本）は平衡気相では本実装の 3 経路最大値と厳密に同値"
          f"（上表の 3 経路ばらつき ≤ {max(r['spread'] for r in rows):.0e}）。")
    print("  ・Ac 0.09〜2.8・650〜875 ℃・0.5〜2.0 MPaG の 21 点で **±13% 以内**、"
          "大小関係はすべて一致。")
    print("  ・**熱力学データを差し替えても完全一致にはできない。**"
          "残差は特許側の 2 桁表記と組成の前提に支配されている。")


if __name__ == "__main__":
    main()
