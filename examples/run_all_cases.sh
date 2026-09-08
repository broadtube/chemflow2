#!/usr/bin/env bash
# 68 ケース（plant3 28 + plant4 40）を一括で解き直し、集計まで通す。
#
# これまで手順は merge_all_cases.py の docstring にコピペ用として書いてあるだけで、
# 実行可能な形が無かった。--reform が 1 条件しか取らないので 8 回に分けて回す必要があり、
# 手で打つと取りこぼす（実際、K3 修正後の再計算で必要になった）。
#
# ────────────────────────────────────────────────────────────────────────
# 使い方
# ────────────────────────────────────────────────────────────────────────
#   bash examples/run_all_cases.sh                     # 既定（下の環境変数を参照）
#   JOBS=3 bash examples/run_all_cases.sh              # 3 並列
#   SV_MODE=plant3 OUT=examples/output/reformer_sv3 \
#     bash examples/run_all_cases.sh                   # plant4 を新規約で回す
#
#   OUT       出力先ディレクトリ（既定 examples/output/reformer）
#   SV_MODE   plant4 の触媒 sizing 規約: per-bed（既定・この検討の元の設計）/ plant3
#             ⚠️ 規約を変えるときは OUT も変えること。ラベルが同じなので上書きになる
#   JOBS      同時実行数（既定 1）。各ソルバは単スレッドなので JOBS ≤ コア数-2 が目安
#   PYTHONPATH  既定 .:../reaction_rate/src（Cantera と reaction_rate の両方が要る）
#
# ⚠️ 既存の JSON / xlsx は上書きされる。前の結果を残すなら OUT を分けるかコピーしておく。
#
# ────────────────────────────────────────────────────────────────────────
# ログ
# ────────────────────────────────────────────────────────────────────────
# merge_all_cases.py の parse_log() は「######## plant3 REFORM=A ########」の見出しで
# 区間を判別する。並列実行すると 1 本のログでは混ざるので、**実行ごとに別ログに書き、
# 最後に正しい順で結合する**。集計自体は JSON から作れるのでログは必須ではない。
#
# 所要時間の実測: 逐次で plant3 約 3.3 h + plant4 約 2.8 h = 6 h 強。
set -euo pipefail

OUT=${OUT:-examples/output/reformer}
SV_MODE=${SV_MODE:-per-bed}
JOBS=${JOBS:-1}
export PYTHONPATH=${PYTHONPATH:-.:../reaction_rate/src}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

COMMON=(--eta 0.5,0.9,0.99 --solve-tol 1e-4 --stop-at-tol)
mkdir -p "$OUT"

run_one() {                       # $1=plant3|plant4  $2=改質条件
    local plant=$1 R=$2 script cases extra=()
    if [ "$plant" = plant3 ]; then
        script=examples/example_reformer_co2removal.py
        cases=baseline,a_upstream,b_recycle
    else
        script=examples/example_reformer_co2removal_plant4.py
        cases=baseline,a_upstream,b_recycle,c_interstage
        extra=(--sv-mode "$SV_MODE" --outdir "$OUT")
    fi
    echo "[$(date +%H:%M:%S)] start $plant REFORM=$R"
    python3 -u "$script" --reform "$R" --cases "$cases" \
        "${COMMON[@]}" "${extra[@]}" > "$OUT/sweep_${plant}_${R}.log" 2>&1
    echo "[$(date +%H:%M:%S)] done  $plant REFORM=$R"
}

for plant in plant3 plant4; do
    for R in A B C D; do
        run_one "$plant" "$R" &
        while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
    done
done
wait

# --- ログを正しい順に結合してから集計 ---
LOG="$OUT/sweep.log"
: > "$LOG"
for plant in plant3 plant4; do
    for R in A B C D; do
        echo "######## $plant REFORM=$R ########" >> "$LOG"
        cat "$OUT/sweep_${plant}_${R}.log" >> "$LOG"
    done
done

n=$(ls "$OUT"/*_metrics.json 2>/dev/null | wc -l)
echo "metrics JSON: $n 件（揃っていれば 68）"
python3 examples/merge_all_cases.py --outdir "$OUT" --log "$LOG"
