"""68 ケースの結合 xlsx を、用途の違う 3 ファイルに組み直す。

merge_xlsx.py が 1 ケース = 1 ファイル（co2removal_{label}_merged.xlsx）まで作る。
このスクリプトはその 68 個を横断して、次の 3 ファイルを書き出す:

    all_summary.xlsx   1 シート  68 行 × 主要指標          … ケース間の比較用
    all_balance.xlsx   68 シート 1 ケース = 1 物質収支表    … 個別ケースを読む用
    all_streams.xlsx   1 シート  縦持ちの全ストリームデータ … ピボット集計用

3 つに分けてあるのは用途が違うため。1 ファイルに混ぜると、比較したいだけのときに
68 シートを抱え込むことになるし、ピボットしたいだけのときに横持ちの表が邪魔になる。

────────────────────────────────────────────────────────────────────────
ケース番号 (No.1〜68)
────────────────────────────────────────────────────────────────────────
樹形図の順に通し番号を振る。この順序が 3 ファイル共通のキーになる:

    1. 反応管の数   plant3 (28 ケース) → plant4 (40 ケース)
    2. 改質条件     A → B → C → D
    3. 吸収塔       baseline → a_upstream → b_recycle → c_interstage
    4. 吸収率 η     0.5 → 0.9 → 0.99

    plant3  A: No. 1〜 7      plant4  A: No.29〜38
            B: No. 8〜14              B: No.39〜48
            C: No.15〜21              C: No.49〜58
            D: No.22〜28              D: No.59〜68

plant3 が 7 ケース刻み・plant4 が 10 ケース刻みなのは、

  - baseline は η に依らないので 1 ケースだけ（3 ケースではない）
  - c_interstage（段間除去）は plant4 にしか無い

の 2 点による。3 + 3 = 6 に baseline 1 を足して 7、c_interstage の 3 を足して 10。

────────────────────────────────────────────────────────────────────────
サマリ指標はどこから来るか
────────────────────────────────────────────────────────────────────────
run_case() が返す指標を、各スクリプトが 2 つの形で出す。

    co2removal_{label}_metrics.json   ← 主。xlsx の隣に必ず残る
    標準出力の「除去位置の比較」表      ← 従。人が読む用（列で並ぶほうが比べやすい）

xlsx から再計算することもできるが、そうすると**定義がスクリプトと二重管理**になる。
「新鮮CO基準 炭素収率」の分母や「CO2 正味生成」の符号のような判断が 2 か所に散ると、
片方だけ直したときに黙ってズレる。どちらも同じ dict から出るので食い違わない。

**当初はログの表だけを読んでいた。** だが表は画面に出るだけなので、`| tee` を
付け忘れて回すと **xlsx は全部揃っているのに指標だけ再現できない**。別の環境で
回した結果を持ってきたときにも同じことが起きる。JSON なら結果と一緒に運べる。

古い結果（JSON が無い）は --export-metrics でログから移植できる。

表の列順は results dict の挿入順 = plan の順 = 上の No. の順なので、列見出しを
解釈する必要がない（見出しは 18 桁固定幅に収まらず隣とくっつくため解釈できない）。
値は "%18.4f" なので、ラベル部を落として float を拾えばそのまま No. 順に並ぶ。

────────────────────────────────────────────────────────────────────────
68 ケースを揃えるまでの手順
────────────────────────────────────────────────────────────────────────
**--reform は 1 条件しか取らない**ので 8 回に分けて実行する（本体スクリプト側の
docstring と同じ理由: 1 回の比較の中で改質条件を固定しないと、改質条件の効果と
除去位置の効果が混ざる）。

サマリ指標は各実行が xlsx の隣に JSON で書くので、**ログを取らなくても集計は
できる**。下の tee は「除去位置の比較」表を後から読み返すためのもので、必須では
ない（要らなければ `2>&1 | tee -a "$LOG"` と直前の echo 行を落とす）。

    cd chemflow2
    LOG=examples/output/reformer/sweep.log
    : > "$LOG"                       # 追記なので最初に空にする
    export PYTHONPATH=.:../reaction_rate/src

    for R in A B C D; do             # plant3: 4 改質条件 × 7 ケース = 28
      echo "######## plant3 REFORM=$R ########" >> "$LOG"
      python3 -u examples/example_reformer_co2removal.py --reform "$R" \
        --cases baseline,a_upstream,b_recycle --eta 0.5,0.9,0.99 \
        --solve-tol 1e-4 --stop-at-tol 2>&1 | tee -a "$LOG"
    done

    for R in A B C D; do             # plant4: 4 改質条件 × 10 ケース = 40
      echo "######## plant4 REFORM=$R ########" >> "$LOG"
      python3 -u examples/example_reformer_co2removal_plant4.py --reform "$R" \
        --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99 \
        --solve-tol 1e-4 --stop-at-tol 2>&1 | tee -a "$LOG"
    done

    python3 examples/merge_all_cases.py          # ← 最後にこれ 1 本

**plant4 は触媒 sizing の規約を --sv-mode で選べる。** 上の手順は既定 (plant3) で
回すので、plant3 と同じ規則になり STY を横比較できる。examples/output/reformer/ に
入っている既存 40 ケースは旧規約 (per-bed) で回したものなので、そちらを再現する
なら --sv-mode per-bed を付ける。規約を変えて回すときは --outdir を分けること
（ラベルが同じなので同じ場所だと上書きになり、新旧の混ざった集計ができる）:

    OUT=examples/output/reformer_sv3
    python3 -u examples/example_reformer_co2removal_plant4.py --reform A \
      --cases baseline,a_upstream,b_recycle,c_interstage --eta 0.5,0.9,0.99 \
      --sv-mode plant3 --outdir "$OUT" --solve-tol 1e-4 --stop-at-tol
    python3 examples/merge_all_cases.py --outdir "$OUT"

plant3 側は規約が変わらないので再計算不要。同じ集計に載せるなら
co2removal_R*_{reformer,plant3,caustic}.xlsx と co2removal_R*_metrics.json を
$OUT にコピーすればよい。

`######## plant3 REFORM=A ########` の見出し行が要るのは、ログのどの区間がどの
(反応管, 改質条件) のものかを parse_log() が判別できないため。スクリプト自身は
自分が plant3 か plant4 かを出力しない。ログを取らないなら不要。

**Windows (cmd) で回すなら**本体スクリプト側の docstring にある set 版を使う
（`set PYTHONPATH=.;../reaction_rate/src` の区切りを円記号にしたもの）。
上の for ループは bash 用。

やり直すと JSON は上書きされる。失敗したケースは JSON が消えるだけで、古い値が
残ることはない（本体スクリプトがケース開始時に消してから解く）。ただし
**結合 xlsx と all_*.xlsx は残る**ので、ケースを減らして回したときは前回の
結果が混ざる。まっさらにするなら output/reformer/ を空にしてから回すこと。

所要時間は実測で plant3 が約 3.3 時間、plant4 が約 2.8 時間（計 6 時間強）。

────────────────────────────────────────────────────────────────────────
実行
────────────────────────────────────────────────────────────────────────
    python3 examples/merge_all_cases.py

引数なしで完結する。既定で

  - ログは examples/output/reformer/sweep.log を読む（上の手順の出力先）
  - merge_xlsx.py の結合を**内部で先に走らせる**（別途叩かなくてよい）

ようにしてあるのは、この 2 つを忘れると**黙って歯抜けのファイルが出る**ため。
結合していなければ物質収支シートが空になり、ログが無ければサマリの指標列が
空になる。どちらもエラーにはならないので、既定で正しい順に走らせている。

    python3 examples/merge_all_cases.py --log other.log   # 別のログを読む
    python3 examples/merge_all_cases.py --no-merge        # 結合済みなら省略できる
    python3 examples/merge_all_cases.py --prefix try1     # try1_summary.xlsx …

要 openpyxl。
"""

import argparse
import os
import re
import sys

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "reformer")
LOG = os.path.join(OUT, "sweep.log")      # 既定のログ（docstring の手順の出力先）

CASES3 = ["baseline", "a_upstream", "b_recycle"]
CASES4 = ["baseline", "a_upstream", "b_recycle", "c_interstage"]
ETAS = [0.5, 0.9, 0.99]

# サマリ表の指標名。**plant3 と plant4 で違う**ので別々に持つ。
# plant4 は 2 段反応器なので R2 入口の組成と触媒 2 床の内訳が増え、代わりに
# plant3 の「反応器入口 CH4」「CO2 パージ」が無い。1 つのリストで両方を読もうと
# すると行数が合わずブロックごと捨てることになる。
METRICS3 = [
    "供給 CO2 [mol/h]",
    "除去 CO2 [mol/h]",
    "反応器入口 [mol/h]",
    "反応器入口 CO2 [mol%]",
    "反応器入口 CO  [mol%]",
    "反応器入口 CH4 [mol%]",
    "全触媒体積 [mL]",
    "酢酸メチル [mol/h]",
    "CO パージ損失 [mol/h]",
    "CO2 パージ [mol/h]",
    "新鮮CO基準 炭素収率 [%]",
    "CO2 正味生成 [mol/h]",
]

METRICS4 = [
    "供給 CO2 [mol/h]",
    "除去 CO2 [mol/h]",
    "反応器入口 [mol/h]",
    "反応器入口 CO2 [mol%]",
    "反応器入口 CO  [mol%]",
    "R2入口 CO  [mol%]",
    "R2入口 CO2 [mol%]",
    "R2入口 H2O [mol%]",
    "触媒 V1 hybrid [mL]",
    "触媒 V2 MA床 [mL]",
    "触媒 合計 [mL]",
    "酢酸メチル [mol/h]",
    "CO パージ損失 [mol/h]",
    "DME パージ [mol/h]",
    "CO2 正味生成 [mol/h]",
    "新鮮CO基準 炭素収率 [%]",
    "投入炭素基準 収率 [%]",
    "CO 収支計 [mol/h]",
    "供給 CO [mol/h]",
]

# 同じ量に別名が付いているものは 1 列に寄せる。触媒総量を両プラントで比べたいのに
# 名前が違うだけで列が割れると、サマリの意味が薄れる。
ALIAS = {"全触媒体積 [mL]": "触媒 合計 [mL]"}

# サマリの指標列 = 両者の和集合（初出順）。片方にしか無い指標は他方で空欄になる。
SUMMARY_METRICS: list[str] = []
for _m in METRICS3 + METRICS4:
    _m = ALIAS.get(_m, _m)
    if _m not in SUMMARY_METRICS:
        SUMMARY_METRICS.append(_m)


# 成分あたりの炭素数。炭素収支を率で見るために要る。
CARBON = {"CO": 1, "CO2": 1, "CH4": 1, "CH3OH": 1, "CH3OCH3": 2,
          "CH3COOCH3": 3, "CH3CHO": 2, "CH3COOH": 2, "Na2CO3": 1}

# 結合 xlsx から算出する指標。ログの 12 指標には無いものだけを足す。
#
# ── 炭素の行き先を「捨」と「戻」で分ける ──────────────────────────
# 吸収塔が抜いた CO2 は**回収**であって損失ではない。改質器へ戻せるので、系から
# 永久に失われるのはパージだけ。この前提の違いで順位が入れ替わるため、両方を出す:
#
#     C効率(捨)  = MA×3 / 投入炭素                 … 除去 CO2 を捨てる前提
#     C効率(戻)  = MA×3 / (投入炭素 − 回収CO2)     … 改質器へ戻す前提
#
# 用語は example_reformer_co2removal_plant4.py の docstring の比較表に合わせてある
# （そこでも「plant3 では戻す前提で評価が反転する」と同じ現象を扱っている）。
# 分母から回収分を引くのは、戻した CO2 は**新鮮炭素として買わなくて済む**ため。
# 改質器の解は H2/CO 制約で決まるので回収 CO2 を戻しても変わらず、新鮮 CO2 の
# 勘定だけが減る、というのが plant4 側 docstring の説明。
#
# ── 率で見る理由 ──────────────────────────────────────────────
# 絶対量だけだと規模の違うケースを比べられない。「CO2 パージ 1.5355 mol/h」は
# baseline では投入炭素の 56% だが、b_recycle η=0.99 では 0.0009 mol/h（0.03%）。
# ただしこれは損失が消えたのではなく**吸収塔へ移った**だけで、回収率の列と
# 並べて初めて読める。
#
# ── 流量列 ────────────────────────────────────────────────────
# 装置サイズが流量で決まるため。触媒体積は既にあるが、循環ガスは圧縮機とループ
# 配管、吸収塔ガス処理量は塔径、パージは放出処理の規模にそれぞれ効く。
DERIVED = [
    "C効率(戻) [%]",
    "C効率(捨) [%]",
    "パージ炭素ロス率 [%]",
    "CO2 回収率 [%]",
    "CO パージ損失率 [%]",
    "循環ガス [mol/h]",
    "循環比 [-]",
    "パージガス [mol/h]",
    "反応器入口 [NL/h]",
    "吸収塔ガス処理量 [mol/h]",
    "NaOH 消費 [mol/h]",
]


# ────────────────────────────────────────────────────────────────
# ケース索引
# ────────────────────────────────────────────────────────────────
def build_index() -> list[dict]:
    """No.1〜68 の索引。label は出力 xlsx のファイル名接頭辞と一致する。"""
    rows: list[dict] = []
    for plant, cases, pref in (("plant3", CASES3, "R"), ("plant4", CASES4, "p4R")):
        for reformer in "ABCD":
            for case in cases:
                etas = [None] if case == "baseline" else ETAS
                for eta in etas:
                    label = (f"{pref}{reformer}_baseline" if eta is None
                             else f"{pref}{reformer}_{case}_eta{eta:g}")
                    rows.append({
                        "no": len(rows) + 1,
                        "plant": plant,
                        "reformer": reformer,
                        "case": case,
                        "eta": eta,
                        "label": label,
                    })
    return rows


def sheet_name(rec: dict) -> str:
    """Excel のシート名 = No. の数字だけ。

    条件はシート先頭行に書いてあるので、タブに詰め込む必要がない。数字だけなら
    all_summary / all_streams の "No." 列からそのまま引ける。
    """
    return str(rec["no"])


# ────────────────────────────────────────────────────────────────
# ログからサマリ指標を拾う
# ────────────────────────────────────────────────────────────────
def parse_log(path: str) -> dict[tuple[str, str], dict]:
    """(plant, 改質条件) → {"values": 指標×ケース, "ok": 成功ラベル順, "status": ...}。

    ログは `######## plant3 REFORM=A ########` で区切られ、各区間の末尾に
    「除去位置の比較」表がある。表が無い区間（実行中・失敗）は結果に載せない。

    表の列は results dict の順 = **成功したケースだけ**が実行順に並ぶ。収束しな
    かったケースがあると列が欠けるので、`[i/N] ラベル 完了/失敗` 行を拾って
    成功ラベルの順序を作り、それと列を突き合わせる。No. 順に決め打ちすると、
    1 件失敗しただけでブロック全体が 1 つずつズレる（あるいは列数不一致で丸ごと
    捨てられる）。
    """
    text = open(path, encoding="utf-8", errors="replace").read()
    blocks: dict[tuple[str, str], dict] = {}

    # 「1回目失敗」は再試行の途中経過なので、最終状態だけを残す
    res = re.compile(r"^\[\d+/\d+\] (\S+) (完了（緩和判定）|完了|1回目失敗|失敗)", re.M)
    head = re.compile(r"^#+ (plant[34]) REFORM=([A-D]) #+$", re.M)
    marks = list(head.finditer(text))
    for i, m in enumerate(marks):
        seg = text[m.end():marks[i + 1].start() if i + 1 < len(marks) else len(text)]

        status, ok = {}, []
        for lb, st in res.findall(seg):
            if st == "1回目失敗":
                status.setdefault(lb, "緩和判定")     # 再試行に入った印
                continue
            if st.startswith("完了"):
                if lb not in ok:
                    ok.append(lb)
                status[lb] = "緩和判定" if "緩和判定" in st or status.get(lb) == "緩和判定" else ""
            else:
                status[lb] = "収束せず"
                if lb in ok:
                    ok.remove(lb)

        if "除去位置の比較" not in seg:
            if status:                                # 実行中でも状態だけは持ち帰る
                blocks[(m.group(1), m.group(2))] = {"values": None, "ok": ok,
                                                    "status": status}
            continue
        tail = seg.split("除去位置の比較", 1)[1]
        metrics = METRICS3 if m.group(1) == "plant3" else METRICS4
        values = []
        for line in tail.splitlines():
            for metric in metrics:
                if line.startswith(" ") and line.strip().startswith(metric):
                    off = len(line) - len(line.lstrip()) + len(metric)
                    values.append([float(x) for x in
                                   re.findall(r"-?\d+\.\d{4}", line[off:])])
                    break
            if len(values) == len(metrics):
                break
        key = (m.group(1), m.group(2))
        new = {"values": values if len(values) == len(metrics) else None,
               "names": metrics, "ok": ok, "status": status}
        # 同じ (反応管, 改質条件) が 2 回出るのは、収束しなかったケースを直して
        # そのブロックを解き直したとき。**表を持つ区間を、持たない区間で上書き
        # しない**（後から 1 ケースだけ回した区間には比較表が出ない）。
        prev = blocks.get(key)
        if prev and prev["values"] is not None and new["values"] is None:
            prev["status"].update(status)
            continue
        if prev:
            new["status"] = {**prev["status"], **status}
        blocks[key] = new
    return blocks


def attach_metrics(index: list[dict], blocks: dict) -> tuple[int, list[str]]:
    """索引に指標と収束状態を書き込む。(埋まった件数, 警告) を返す。"""
    warn, filled = [], 0
    by_label = {r["label"]: r for r in index}

    for (plant, reformer), blk in sorted(blocks.items()):
        for lb, st in blk["status"].items():
            if lb in by_label:
                by_label[lb]["備考"] = st
            else:
                warn.append(f"{plant} {reformer}: ログのラベル {lb} が索引に無い")

        values, ok = blk["values"], blk["ok"]
        if values is None:
            continue
        if len(values[0]) != len(ok):
            warn.append(f"{plant} {reformer}: 表 {len(values[0])} 列 / 成功 {len(ok)} 件"
                        f" — 突き合わせずスキップ")
            continue
        for col, lb in enumerate(ok):
            rec = by_label.get(lb)
            if rec is None:
                continue
            for name, row in zip(blk["names"], values):
                rec[ALIAS.get(name, name)] = row[col]
            filled += 1

    missing = [f"{p} {r}" for p in ("plant3", "plant4") for r in "ABCD"
               if blocks.get((p, r), {}).get("values") is None]
    if missing:
        warn.append(f"サマリ表が未了のブロック: {', '.join(missing)}")
    bad = [r["label"] for r in index if r.get("備考") == "収束せず"]
    if bad:
        warn.append(f"収束しなかったケース {len(bad)} 件: {', '.join(bad)}")
    soft = [r["label"] for r in index if r.get("備考") == "緩和判定"]
    if soft:
        warn.append(f"緩和判定で解いたケース {len(soft)} 件: {', '.join(soft)}")
    return filled, warn


def attach_sidecars(index: list[dict], outdir: str) -> int:
    """各ケースが書いた co2removal_{label}_metrics.json を読む。

    指標名は plant3 / plant4 で違うが、JSON は dict なので名前で引ける
    （ログの表と違って列位置に依存しない）。ALIAS で同義の列に寄せる。
    """
    import json
    n = 0
    for rec in index:
        p = os.path.join(outdir, f"co2removal_{rec['label']}_metrics.json")
        if not os.path.exists(p):
            continue
        with open(p, encoding="utf-8") as f:
            d = json.load(f)
        for k, v in d.get("metrics", {}).items():
            rec[ALIAS.get(k, k)] = v
        if d.get("note"):
            rec["備考"] = d["note"]
        n += 1
    return n


# ────────────────────────────────────────────────────────────────
# 結合済み xlsx の読み取り
# ────────────────────────────────────────────────────────────────
def read_merged(path: str) -> dict:
    """merge_xlsx.py が書いた merged シートを構造化して返す。

    レイアウト（merge_xlsx.merge_label と対応）:
        行1 ブロック帯 / 行2 グローバル番号 / 行3 ストリーム名
        行4 以降が本体。A 列 = ラベル、B 列 = MW、C 列以降が各ストリーム。
    """
    ws = load_workbook(path, data_only=True).active
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    band, numbers, names = rows[0], rows[1], rows[2]
    body = rows[3:]

    # 帯はセル結合で先頭列にしか値が無いので、右へ引き伸ばす
    bands, cur = [], ""
    for v in band:
        cur = v if v else cur
        bands.append(cur)

    streams = []
    for c in range(2, len(names)):
        if names[c] is None:
            continue
        streams.append({"col": c, "block": bands[c] or "",
                        "gno": numbers[c], "name": names[c]})

    # A 列を走査して T/P/Phase と各単位セクションの成分行を拾う
    meta, sections, totals, unit, comps = {}, {}, {}, None, []
    for r, row in enumerate(body):
        a = row[0]
        if a in ("T [°C]", "P", "Phase"):
            meta[a] = r
        elif isinstance(a, str) and a.startswith("[") and a.endswith("]"):
            if unit:
                sections[unit] = comps
            unit, comps = a.strip("[]"), []
        elif isinstance(a, str) and a.startswith("total"):
            if unit:
                totals[unit] = r
        elif a and unit:
            comps.append((a, r))
    if unit:
        sections[unit] = comps

    return {"rows": rows, "body": body, "streams": streams,
            "meta": meta, "sections": sections, "totals": totals}


# ────────────────────────────────────────────────────────────────
# 結合 xlsx から算出する指標
# ────────────────────────────────────────────────────────────────
def derive(data: dict) -> dict:
    """DERIVED の指標を 1 ケース分計算する。

    ストリームは**名前の番号を落として**引く。plant3 と plant4 で採番が違う
    （Purge は plant3 が ⑦、plant4 が ⑫）が、名前部分は共通だから。
    ブロックの区別はグローバル通番の百の位で付く（100=改質器 / 200,400=吸収塔 /
    300=プラント）。帯の文字列で判定するとケース種別ごとに表記が変わって壊れる。
    """
    body, sections, totals = data["body"], data["sections"], data["totals"]

    def find(block: int, name: str):
        for st in data["streams"]:
            gno = st["gno"]
            if not isinstance(gno, int) or gno // 100 != block:
                continue
            if str(st["name"]).split(". ", 1)[-1] == name:
                return st["col"]
        return None

    def cell(col, r):
        if col is None or r is None:
            return None
        row = body[r]
        v = row[col] if col < len(row) else None
        return v if isinstance(v, (int, float)) else None

    def total(col, unit="mol/h"):
        return cell(col, totals.get(unit))

    def comp(col, name, unit="mol/h"):
        r = dict(sections.get(unit, [])).get(name)
        return cell(col, r)

    def carbon(col):
        """そのストリームが運ぶ炭素 [mol-C/h]。"""
        if col is None:
            return None
        tot = 0.0
        for name, r in sections.get("mol/h", []):
            n = CARBON.get(name)
            v = cell(col, r)
            if n and v:
                tot += n * v
        return tot

    feed = find(1, "Feed")                  # 改質器の原料 = 系に入る炭素の全量
    steam1 = find(3, "Steam1")              # プラント入口
    rin = find(3, "ReactorIn")
    purge = find(3, "Purge")
    recyc = find(3, "RecycleGas")
    gasin = find(2, "GasIn") or find(4, "GasIn")        # 吸収塔（前段/後段）
    caustic = find(2, "Caustic") or find(4, "Caustic")

    c_in = carbon(feed)
    out: dict = {}

    def pct(x, base):
        return None if not base or x is None else x / base * 100

    ma_c = (comp(find(3, "MethylAcetate"), "CH3COOCH3") or 0) * 3
    # 吸収塔が抜いた CO2 = 回収量。CO2Vent は結合シートに列が無いケースがあるので、
    # 吸収塔の入口−出口で数える（baseline は吸収塔が無いので 0）。
    scrub = find(2, "ScrubbedGas") or find(4, "ScrubbedGas")
    removed = 0.0
    if gasin is not None and scrub is not None:
        gi, so = comp(gasin, "CO2"), comp(scrub, "CO2")
        if gi is not None and so is not None:
            removed = gi - so

    out["C効率(戻) [%]"] = pct(ma_c, None if c_in is None else c_in - removed)
    out["C効率(捨) [%]"] = pct(ma_c, c_in)
    out["パージ炭素ロス率 [%]"] = pct(carbon(purge), c_in)
    out["CO2 回収率 [%]"] = pct(removed, c_in)
    out["CO パージ損失率 [%]"] = pct(comp(purge, "CO"), comp(steam1, "CO"))

    out["循環ガス [mol/h]"] = total(recyc)
    fresh = total(steam1)
    rec_v = total(recyc)
    out["循環比 [-]"] = None if not fresh or rec_v is None else rec_v / fresh
    out["パージガス [mol/h]"] = total(purge)
    out["反応器入口 [NL/h]"] = total(rin, "NL/h")
    out["吸収塔ガス処理量 [mol/h]"] = total(gasin)
    out["NaOH 消費 [mol/h]"] = comp(caustic, "NaOH")
    return out


# ────────────────────────────────────────────────────────────────
# 出力 1: サマリ
# ────────────────────────────────────────────────────────────────
def write_summary(index: list[dict], path: str) -> None:
    wb = Workbook()
    ws = wb.active
    ws.title = "summary"
    head = (["No.", "反応管", "改質条件", "吸収塔", "吸収率 η", "ラベル", "備考"]
            + SUMMARY_METRICS + DERIVED)
    ws.append(head)
    for c in ws[1]:
        c.font = Font(bold=True)
        c.alignment = Alignment(horizontal="center", wrap_text=True)
    for rec in index:
        ws.append([rec["no"], rec["plant"], rec["reformer"], rec["case"],
                   rec["eta"], rec["label"], rec.get("備考", "")]
                  + [rec.get(m) for m in SUMMARY_METRICS + DERIVED])
    ws.freeze_panes = "H2"
    for i, w in enumerate([6, 9, 10, 14, 10, 24, 11], 1):
        ws.column_dimensions[get_column_letter(i)].width = w
    for i in range(8, 8 + len(SUMMARY_METRICS) + len(DERIVED)):
        ws.column_dimensions[get_column_letter(i)].width = 15
    ws.auto_filter.ref = ws.dimensions
    wb.save(path)


# ────────────────────────────────────────────────────────────────
# 出力 2: 68 シートの物質収支
# ────────────────────────────────────────────────────────────────
def write_balance(index: list[dict], merged: dict[int, dict], path: str) -> None:
    wb = Workbook()
    wb.remove(wb.active)
    bold, center = Font(bold=True), Alignment(horizontal="center")
    for rec in index:
        data = merged.get(rec["no"])
        ws = wb.create_sheet(sheet_name(rec))
        eta = "—" if rec["eta"] is None else f"{rec['eta']:g}"
        ws.append([f"No.{rec['no']}", rec["plant"], f"改質条件 {rec['reformer']}",
                   rec["case"], f"η={eta}", rec["label"]])
        ws["A1"].font = bold
        if data is None:
            ws.append(["(結合 xlsx が見つかりません)"])
            continue
        ws.append([])
        for r, row in enumerate(data["rows"]):
            ws.append(list(row))
            if r < 3:
                for c in ws[ws.max_row]:
                    c.font, c.alignment = bold, center
        ws.freeze_panes = "C6"
        ws.column_dimensions["A"].width = 18
        ws.column_dimensions["B"].width = 8
        for i in range(3, len(data["rows"][2]) + 1):
            ws.column_dimensions[get_column_letter(i)].width = 13
    wb.save(path)


# ────────────────────────────────────────────────────────────────
# 出力 3: 縦持ちの全ストリームデータ
# ────────────────────────────────────────────────────────────────
def write_long(index: list[dict], merged: dict[int, dict], path: str) -> list[str]:
    """1 行 = (ケース, ストリーム, 成分)。単位は列に展開する。

    単位を行に落とす（No/stream/component/unit/value の完全縦持ち）と行数が 4 倍に
    なるうえ、ピボットで「mol/h と mol% を並べて見る」のに unit をフィルタする手間が
    増える。単位は 4 つで固定なので列に置くほうが素直。
    """
    wb = Workbook()
    ws = wb.active
    ws.title = "streams"
    units = ["mol/h", "mol%", "g/h", "NL/h"]
    head = ["No.", "反応管", "改質条件", "吸収塔", "吸収率 η", "ラベル",
            "ブロック", "通番", "ストリーム", "T [°C]", "P", "Phase", "成分"] + units
    ws.append(head)
    for c in ws[1]:
        c.font = Font(bold=True)
        c.alignment = Alignment(horizontal="center", wrap_text=True)

    warn, n = [], 0
    for rec in index:
        data = merged.get(rec["no"])
        if data is None:
            warn.append(rec["label"])
            continue
        body, meta, sections = data["body"], data["meta"], data["sections"]
        base = [rec["no"], rec["plant"], rec["reformer"], rec["case"],
                rec["eta"], rec["label"]]
        for st in data["streams"]:
            c = st["col"]

            def at(r):
                row = body[r]
                return row[c] if c < len(row) else None

            head3 = [st["block"], st["gno"], st["name"],
                     at(meta["T [°C]"]) if "T [°C]" in meta else None,
                     at(meta["P"]) if "P" in meta else None,
                     at(meta["Phase"]) if "Phase" in meta else None]
            # 成分行は単位セクションで同じ並びなので mol/h を基準にする
            for comp, r in sections.get("mol/h", []):
                vals = []
                for u in units:
                    rr = dict(sections.get(u, [])).get(comp)
                    vals.append(at(rr) if rr is not None else None)
                ws.append(base + head3 + [comp] + vals)
                n += 1

    ws.freeze_panes = "N2"
    for i, w in enumerate([6, 9, 10, 14, 10, 24, 30, 8, 18, 10, 11, 8, 14], 1):
        ws.column_dimensions[get_column_letter(i)].width = w
    for i in range(14, 14 + len(units)):
        ws.column_dimensions[get_column_letter(i)].width = 14
    ws.auto_filter.ref = ws.dimensions
    wb.save(path)
    print(f"  縦持ち {n} 行")
    return warn


# ────────────────────────────────────────────────────────────────
def main() -> int:
    ap = argparse.ArgumentParser(
        description="68 ケースを 3 ファイル（サマリ / 物質収支 / 縦持ち）に組み直す",
        epilog="引数なしで実行すれば結合から出力まで通しで走る。")
    ap.add_argument("--outdir", default=OUT, help=f"入出力ディレクトリ（既定 {OUT}）")
    ap.add_argument("--log", default=LOG,
                    help=f"スイープのログ（既定 {LOG}）。サマリ指標の取得元")
    ap.add_argument("--no-merge", action="store_true",
                    help="merge_xlsx.py の結合を省く（結合済みのときだけ）")
    ap.add_argument("--prefix", default="all", help="出力ファイル名の接頭辞（既定 all）")
    ap.add_argument("--export-metrics", action="store_true",
                    help="ログから読んだ指標を JSON サイドカーに書き出す"
                         "（tee を付けずに回した古い結果を移植するとき用）")
    args = ap.parse_args()

    index = build_index()
    print(f"ケース索引 {len(index)} 件"
          f"（plant3 {sum(r['plant'] == 'plant3' for r in index)} / "
          f"plant4 {sum(r['plant'] == 'plant4' for r in index)}）")

    # ── 1. ケースごとの結合（merge_xlsx.py と同じ処理を先に済ませる）──────
    # 別コマンドにすると忘れたときに物質収支シートが黙って空になる。
    if not args.no_merge:
        import merge_xlsx                       # 同じ examples/ にある
        labels = merge_xlsx.discover_labels(args.outdir)
        done = [lb for lb in labels
                if merge_xlsx.merge_label(lb, args.outdir) is not None]
        print(f"ケース結合: {len(done)} / {len(labels)} ケース")

    # ── 2. サマリ指標: JSON サイドカーを主、ログを従とする ──────────
    # サイドカーは各ケースが xlsx の隣に書くので、結果と一緒に必ず残る。ログは
    # `| tee` を付けて回したときにしか無いので、古い結果を読むときの保険。
    n_json = attach_sidecars(index, args.outdir)
    print(f"サマリ指標: {n_json} / {len(index)} ケース分を JSON から取得")
    if n_json < len(index):
        if os.path.exists(args.log):
            filled, warn = attach_metrics(index, parse_log(args.log))
            print(f"  ログから追加: {filled} ケース分（{args.log}）")
            for w in warn:
                print(f"  ⚠ {w}")
        elif n_json == 0:
            print(f"  ⚠ JSON もログも無い — サマリの指標列は空になります\n"
                  f"     計算をやり直すか、計算した環境から "
                  f"co2removal_*_metrics.json をコピーしてください")

    if args.export_metrics:
        import json
        n = 0
        for rec in index:
            # 各スクリプトが書く JSON と同じ名前で出す。サマリ側は ALIAS で
            # 寄せた名前（触媒 合計）を持っているので、plant3 の元の名前
            # （全触媒体積）に戻してから書く。混在すると diff が読みにくい。
            names = METRICS3 if rec["plant"] == "plant3" else METRICS4
            got = {n: rec[ALIAS.get(n, n)] for n in names
                   if ALIAS.get(n, n) in rec}
            if not got:
                continue
            with open(os.path.join(args.outdir,
                                   f"co2removal_{rec['label']}_metrics.json"),
                      "w", encoding="utf-8") as f:
                json.dump({"label": rec["label"], "plant": rec["plant"],
                           "reform": rec["reformer"], "case": rec["case"],
                           "eta": rec["eta"], "note": rec.get("備考", ""),
                           "metrics": got}, f, ensure_ascii=False, indent=1)
            n += 1
        print(f"✓ metrics.json を {n} 件書き出しました")

    # ── 3. 結合済み xlsx を読む ───────────────────────────────
    merged, missing = {}, []
    for rec in index:
        p = os.path.join(args.outdir, f"co2removal_{rec['label']}_merged.xlsx")
        if os.path.exists(p):
            merged[rec["no"]] = read_merged(p)
            rec.update(derive(merged[rec["no"]]))
        else:
            missing.append(f"No.{rec['no']} {rec['label']}")
    print(f"結合 xlsx: {len(merged)} / {len(index)} ケース")
    if missing:
        print(f"  ⚠ 未検出 {len(missing)} 件: {', '.join(missing[:6])}"
              + (" …" if len(missing) > 6 else ""))

    # ── 4. 書き出し ───────────────────────────────────────
    os.makedirs(args.outdir, exist_ok=True)
    p1 = os.path.join(args.outdir, f"{args.prefix}_summary.xlsx")
    p2 = os.path.join(args.outdir, f"{args.prefix}_balance.xlsx")
    p3 = os.path.join(args.outdir, f"{args.prefix}_streams.xlsx")

    write_summary(index, p1)
    print(f"✓ {os.path.basename(p1)}  (1 シート / {len(index)} 行)")
    write_balance(index, merged, p2)
    print(f"✓ {os.path.basename(p2)}  ({len(index)} シート)")
    write_long(index, merged, p3)
    print(f"✓ {os.path.basename(p3)}")
    return 0 if len(merged) == len(index) else 1


if __name__ == "__main__":
    sys.exit(main())
