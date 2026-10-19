#!/bin/bash
# =============================================================================
# 全量消融训练：数据/h5data 下所有数据集（可排除） × 6 种消融（不含 BiLSTM）
# 消融：baseline, no_se, no_multiscale, no_gru, no_combined_loss, no_conservative
# 输出：消融实验/全量消融运行/模型权重、各个具体实验的指标、运行时间等
# 默认槽位：GPU0=15，其它每卡=20（15 20 20 20，共 75 槽），常驻并行；完成一个补一个
# =============================================================================
set +e
# 任意退出时在日志中打印退出码，便于排查「自动退出」（nohup 会把这些打到 run_ablation.log）
trap 'echo "[全量消融] 主进程退出, exit_code=$?"' EXIT
BASE="/data/CARBON_PROJECT"
DATA_ROOT="${BASE}/数据/h5data"
OUT_ROOT="${BASE}/消融实验/全量消融运行"
MODEL_12="${BASE}/MARCP"

# 消融列表（不含 BiLSTM）：简称 -> 传给 marcp_train 的额外参数
# baseline: 无额外参数
# 构建任务列表时用 ABLATION_LIST_SCHEDULE：非 baseline 优先，避免本批多为「仅缺 baseline」时只跑少量 GPU、槽位闲置
ABLATION_BASELINE="baseline"
ABLATION_LIST="baseline no_se no_multiscale no_gru no_combined_loss no_conservative"
ABLATION_LIST_SCHEDULE="no_se no_multiscale no_gru no_combined_loss no_conservative baseline"
ABLATION_FLAG_no_se="--ablation_no_se"
ABLATION_FLAG_no_multiscale="--ablation_no_multiscale"
ABLATION_FLAG_no_gru="--ablation_no_gru"
ABLATION_FLAG_no_combined_loss="--ablation_no_combined_loss"
ABLATION_FLAG_no_conservative="--ablation_no_conservative"

# GPU 槽位：GPU0 最多 SLOTS_GPU0，其它每卡 SLOTS_OTHER（4 卡共 15+20+20+20=75 槽）
SLOTS_GPU0="${SLOTS_GPU0:-15}"
SLOTS_OTHER="${SLOTS_OTHER:-20}"
# 排除的任务（数据已删除或不做消融）；空格分隔，如 "TASK1 TASK2"
TASK_EXCLUDE="${TASK_EXCLUDE:-HNRNPC_HepG2_ENCSR550DVK HNRNPU_K562_ENCSR520BZQ}"
SORT_TASKS_BY_SIZE="${SORT_TASKS_BY_SIZE:-1}"
CLEAR_BEFORE_RETRAIN="${CLEAR_BEFORE_RETRAIN:-0}"
# 每进程 DataLoader workers；槽位≥24 时脚本会自动降为 2 以防 OOM 杀进程
DATA_WORKERS_PER_JOB="${DATA_WORKERS_PER_JOB:-0}"
RESOURCE_FRAC_PER_JOB="${RESOURCE_FRAC_PER_JOB:-0.9}"
# 批大小：默认 256，保证与 baseline 公平对比；显存充足时可改为 512、1024 等
BATCH_SIZE="${BATCH_SIZE:-256}"
# 固定任务列表：指向 batch11_tasks.txt 等文件时，本批仅跑文件内各行，且不跳过「模型+指标已存在」项（用于指定 N 条同批并发重跑）
# 由 run_ablation_background.sh 在 USE_FIXED_BATCH11=1 时自动设为 batch11_tasks.txt

VAL_CHRS="2,9,16"
TEST_CHRS="1,8,15"
TRAIN_CHRS="3,4,5,6,7,10,11,12,13,14,17,18,19,20,21,22"

if [ -z "${N_GPU}" ]; then
    N_GPU=$(nvidia-smi -L 2>/dev/null | wc -l)
    [ "$N_GPU" -lt 1 ] && N_GPU=1
fi
# 可选：仅使用指定 GPU（与其它任务隔离，避免本程序一直抢不到卡）。例: ABLATION_GPU_IDS=2,3
ABLATION_GPU_IDS="${ABLATION_GPU_IDS:-}"
if [ -n "$ABLATION_GPU_IDS" ]; then
    N_GPU=0
    for g in ${ABLATION_GPU_IDS//,/ }; do [ -n "$g" ] && N_GPU=$((N_GPU + 1)); done
    [ "$N_GPU" -lt 1 ] && N_GPU=1
    echo "ABLATION_GPU_IDS=$ABLATION_GPU_IDS → 仅使用以上 $N_GPU 张 GPU，其它卡留给其它任务"
fi
export ABLATION_GPU_IDS

# 确保输出根目录为绝对路径（后台子进程可能继承不到相对路径）
OUT_ROOT="$(cd "${OUT_ROOT}" 2>/dev/null && pwd)" || true
[ -z "$OUT_ROOT" ] && OUT_ROOT="${BASE}/消融实验/全量消融运行"
mkdir -p "${OUT_ROOT}/模型权重/MARCP"
mkdir -p "${OUT_ROOT}/各个具体实验的指标/MARCP"
mkdir -p "${OUT_ROOT}/logs"
WORKSPACE_ROOT="${OUT_ROOT}/workspace/MARCP"
mkdir -p "${WORKSPACE_ROOT}"
# 若存在 venv 则使用其 Python（需先执行 bash setup_venv.sh）
VENV_PYTHON=""
if [ -x "${OUT_ROOT}/venv/bin/python3" ]; then
    VENV_PYTHON="${OUT_ROOT}/venv/bin/python3"
    echo "使用 venv: ${OUT_ROOT}/venv"
fi
export VENV_PYTHON

PID_FILE="${OUT_ROOT}/.run_ablation.pid"
STATE_FILE="${OUT_ROOT}/.run_ablation.state"
COMPLETED_LOG="${OUT_ROOT}/.run_ablation.completed"
COMPLETED_LOCK="${OUT_ROOT}/.run_ablation.completed.lock"
SUMMARY="${OUT_ROOT}/汇总说明.md"
SUMMARY_LOCK="${OUT_ROOT}/.summary.lock"

if [ "$CLEAR_BEFORE_RETRAIN" = "1" ]; then
    echo "CLEAR_BEFORE_RETRAIN=1: 清空旧权重、指标与 workspace..."
    rm -rf "${OUT_ROOT}/模型权重/MARCP"/*
    rm -rf "${OUT_ROOT}/各个具体实验的指标/MARCP"/*
    rm -rf "${WORKSPACE_ROOT:?}"/*
    rm -f "$COMPLETED_LOG" "$COMPLETED_LOCK"
    echo "已清空。"
fi

# 尽早写入 PID，避免建任务列表失败或异常退出时 query 显示「未运行」；仅正常结束或 stop 时再删
echo $$ > "$PID_FILE"
rm -f "${OUT_ROOT}/.run_ablation.done"

EPOCHS=100

get_tasks() {
    for d in "${DATA_ROOT}"/*/; do
        [ -d "$d" ] || continue
        task=$(basename "$d")
        [ -d "${d}/_5ch" ] || continue
        echo "$task"
    done | sort -u
}

file_list() {
    local task=$1
    local chrs=$2
    local sep=""
    local out=""
    for c in ${chrs//,/ }; do
        f="${DATA_ROOT}/${task}/_5ch/${task}_chr${c}_5ch.h5"
        [ -f "$f" ] || continue
        out="${out}${sep}${f}"
        sep=","
    done
    echo "$out"
}

discover_files() {
    local task=$1
    local dir="${DATA_ROOT}/${task}/_5ch"
    local pattern="${task}_chr*_5ch.h5"
    if [ ! -d "$dir" ]; then
        echo ""; echo ""; echo ""
        return
    fi
    local train="" val="" test="" st="" sv="" ss=""
    for f in "${dir}"/${pattern}; do
        [ -f "$f" ] || continue
        base=$(basename "$f")
        c=$(echo "$base" | sed -n "s/.*_chr\([0-9]*\)_5ch.h5/\1/p")
        [ -n "$c" ] || continue
        case "$c" in
            2|9|16) val="${val}${sv}${f}"; sv="," ;;
            1|8|15) test="${test}${ss}${f}"; ss="," ;;
            *)      train="${train}${st}${f}"; st="," ;;
        esac
    done
    echo "$train"
    echo "$val"
    echo "$test"
}

# 返回某消融的 MARCP 额外参数（baseline 为空）
get_ablation_flag() {
    local ab="$1"
    case "$ab" in
        no_se) echo "--ablation_no_se" ;;
        no_multiscale) echo "--ablation_no_multiscale" ;;
        no_gru) echo "--ablation_no_gru" ;;
        no_combined_loss) echo "--ablation_no_combined_loss" ;;
        no_conservative) echo "--ablation_no_conservative" ;;
        *) echo "" ;;
    esac
}

TASKS=$(get_tasks)
if [ -z "$TASKS" ]; then
    echo "未找到任何任务（需在 ${DATA_ROOT} 下存在 <TASK>/_5ch）。退出."
    exit 0
fi

# 按 GPU 槽位限制分配：GPU0 最多 SLOTS_GPU0 个，其它每卡最多 SLOTS_OTHER 个，避免全挤到一张卡
# 使用 .gpu_slot_counts 记录每卡当前任务数，.gpu_done 记录完成时归还的 gpu_id
GPU_COUNTS_FILE="${OUT_ROOT}/.gpu_slot_counts"
GPU_DONE_FILE="${OUT_ROOT}/.gpu_done"
GPU_COUNTS_LOCK="${OUT_ROOT}/.gpu_slot_counts.lock"
GPU_DONE_LOCK="${OUT_ROOT}/.gpu_done.lock"
# 总并发槽位数
TOTAL_SLOTS=$((SLOTS_GPU0 + (N_GPU - 1) * SLOTS_OTHER))
[ "$TOTAL_SLOTS" -lt 1 ] && TOTAL_SLOTS=1
# 槽位多时总进程数 = 槽位×(1+DATA_WORKERS_PER_JOB)，易触发系统 OOM 杀进程；自动降低每任务 worker 数
if [ "$TOTAL_SLOTS" -ge 40 ] && [ "$DATA_WORKERS_PER_JOB" -gt 1 ]; then
    echo "[OOM 防护] 槽位=$TOTAL_SLOTS，将每任务 DataLoader 进程数由 $DATA_WORKERS_PER_JOB 降为 1，避免内存耗尽导致 DataLoader worker 被系统杀死；需更多 worker 可显式设 DATA_WORKERS_PER_JOB 后重跑"
    DATA_WORKERS_PER_JOB=1
elif [ "$TOTAL_SLOTS" -ge 24 ] && [ "$DATA_WORKERS_PER_JOB" -gt 2 ]; then
    echo "[OOM 防护] 槽位=$TOTAL_SLOTS，将每任务 DataLoader 进程数由 $DATA_WORKERS_PER_JOB 降为 2，避免内存耗尽导致进程被系统杀死；需更多 worker 可显式设 DATA_WORKERS_PER_JOB 后重跑"
    DATA_WORKERS_PER_JOB=2
fi
echo "常驻并行: 最多 $TOTAL_SLOTS 个任务，GPU0 最多 $SLOTS_GPU0 个、其它每卡最多 $SLOTS_OTHER 个，按槽位分配"

# 汇总表头
echo -e "# 全量消融训练汇总\n\n| 任务 | 消融 | Mean_PCC | Spearman | Peak_AUC | Peak_AP | 训练时长(min) | 备注 |" > "$SUMMARY"
echo -e "|------|------|----------|----------|----------|---------|----------------|------|" >> "$SUMMARY"

# 构建任务列表：每行 task|train_5ch|val_5ch|test_5ch|ablation_short|ablation_flag|run_index
TASK_LIST_FILE=$(mktemp)
trap "rm -f '$TASK_LIST_FILE'" EXIT
run_index=0
_FIXED_LIST_ACTIVE=0

_fixed_list="${ABLATION_FIXED_LIST_FILE:-}"
if [ -n "$_fixed_list" ] && [ ! -f "$_fixed_list" ] && [ -f "${OUT_ROOT}/$(basename "$_fixed_list")" ]; then
    _fixed_list="${OUT_ROOT}/$(basename "$_fixed_list")"
fi

if [ -n "$_fixed_list" ] && [ -f "$_fixed_list" ]; then
    grep -v '^[[:space:]]*#' "$_fixed_list" | grep -v '^[[:space:]]*$' > "$TASK_LIST_FILE"
    n_tasks_after_exclude=$(wc -l < "$TASK_LIST_FILE")
    _FIXED_LIST_ACTIVE=1
    echo "任务来源: 固定列表 $_fixed_list（共 $n_tasks_after_exclude 条；不跳过「已完成」，全部进入本批调度以利同批并发）"
else
    if [ -n "${ABLATION_FIXED_LIST_FILE:-}" ] && [ ! -f "${_fixed_list}" ]; then
        echo "ERROR: ABLATION_FIXED_LIST_FILE=${ABLATION_FIXED_LIST_FILE} 未找到文件（已尝试: ${_fixed_list}）。" >&2
        exit 1
    fi
    # 先按数据量排序任务，再对每个任务生成各消融
    TASK_SIZES=$(mktemp)
    for t in $TASKS; do
        case " $TASK_EXCLUDE " in *" $t "*) continue ;; esac
        sz=$(du -sb "${DATA_ROOT}/${t}/_5ch" 2>/dev/null | awk '{print $1}')
        [ -z "$sz" ] && sz=0
        echo "${sz}|${t}"
    done > "$TASK_SIZES"
    n_tasks_after_exclude=$(wc -l < "$TASK_SIZES")
    echo "任务来源: ${DATA_ROOT}，排除: ${TASK_EXCLUDE:-（无）}，参与消融任务数: $n_tasks_after_exclude"
    if [ "$SORT_TASKS_BY_SIZE" = "1" ]; then
        sort -t'|' -k1 -rn "$TASK_SIZES" > "${TASK_SIZES}.sorted"
        mv "${TASK_SIZES}.sorted" "$TASK_SIZES"
        echo "任务已按数据量从大到小排序（大任务优先）"
    fi

    while IFS='|' read -r sz task; do
        [ -n "$task" ] || continue
        train_5ch=$(file_list "$task" "$TRAIN_CHRS")
        val_5ch=$(file_list "$task" "$VAL_CHRS")
        test_5ch=$(file_list "$task" "$TEST_CHRS")
        if [ -z "$train_5ch" ] || [ -z "$val_5ch" ] || [ -z "$test_5ch" ]; then
            { read -r train_5ch; read -r val_5ch; read -r test_5ch; } < <(discover_files "$task")
        fi
        if [ -z "$train_5ch" ] || [ -z "$val_5ch" ] || [ -z "$test_5ch" ]; then
            echo "[SKIP] $task: 缺少 5ch 文件"
            echo "| $task | - | - | - | - | - | - | 缺少数据 |" >> "$SUMMARY"
            continue
        fi
        for ab in $ABLATION_LIST_SCHEDULE; do
            weight_dir="${OUT_ROOT}/模型权重/MARCP/${task}/${ab}"
            metric_dir="${OUT_ROOT}/各个具体实验的指标/MARCP/${task}/${ab}"
            # 以指标文件为准判定「已完成」；仅权重存在而指标缺失时仍会重跑（与 marcp_train 成功写出 summary 对齐）
            if [ "${FORCE_RERUN}" != "1" ] && [ -f "${metric_dir}/test_metrics_summary.txt" ]; then
                echo "[SKIP] $task $ab: 已完成"
                continue
            fi
            run_index=$((run_index + 1))
            flag=$(get_ablation_flag "$ab")
            echo "${task}|${train_5ch}|${val_5ch}|${test_5ch}|${ab}|${flag}|${run_index}" >> "$TASK_LIST_FILE"
        done
    done < "$TASK_SIZES"
    rm -f "$TASK_SIZES"
fi

total_to_run=$(wc -l < "$TASK_LIST_FILE")
if [ "$total_to_run" -eq 0 ]; then
    echo "所有消融任务均已完成或已跳过，无需运行。"
    if [ "$n_tasks_after_exclude" -eq 0 ]; then
        echo "原因: 参与任务数=0（数据目录无 _5ch 或全部被 TASK_EXCLUDE 排除）。数据目录: ${DATA_ROOT}，排除: ${TASK_EXCLUDE:-（无）}"
    else
        echo "原因: 本批 $n_tasks_after_exclude 个任务×6 种消融中，均已存在 模型权重 与 指标 文件。若需强制重跑全部，请设置 FORCE_RERUN=1 后重新执行: FORCE_RERUN=1 bash run_ablation_background.sh"
    fi
    # 仍写入 STATE 便于 query 显示「本批 0」
    echo "start $(date +%s)" > "$STATE_FILE"
    echo "total 0" >> "$STATE_FILE"
    echo "total_slots $TOTAL_SLOTS" >> "$STATE_FILE"
    echo "gpu_slots_desc GPU0最多${SLOTS_GPU0}、其它每卡${SLOTS_OTHER}(共${TOTAL_SLOTS}槽)" >> "$STATE_FILE"
    echo "jobs" >> "$STATE_FILE"
    exit 0
fi

# 立即写入 STATE，便于 query 在启动阶段就能看到本批任务数（避免读到上一轮的 total 0）
echo "start $(date +%s)" > "$STATE_FILE"
echo "total $total_to_run" >> "$STATE_FILE"
echo "n_gpu $N_GPU" >> "$STATE_FILE"
echo "total_slots $TOTAL_SLOTS" >> "$STATE_FILE"
echo "gpu_slots_desc GPU0最多${SLOTS_GPU0}、其它每卡${SLOTS_OTHER}(共${TOTAL_SLOTS}槽)" >> "$STATE_FILE"
echo "jobs" >> "$STATE_FILE"
cat "$TASK_LIST_FILE" >> "$STATE_FILE"

export OUT_ROOT MODEL_12 EPOCHS WORKSPACE_ROOT COMPLETED_LOG COMPLETED_LOCK SUMMARY SUMMARY_LOCK VENV_PYTHON GPU_DONE_FILE GPU_DONE_LOCK

_run_one() {
    local task="$1" train_5ch="$2" val_5ch="$3" test_5ch="$4" gpu_id="$5" ablation_short="$6" ablation_flag="$7" run_index="$8" total_runs="$9" base_out="${10}"
    trap '( flock 201; echo "$gpu_id" >> "$GPU_DONE_FILE"; ) 201>"$GPU_DONE_LOCK"' EXIT
    [ -z "$base_out" ] && base_out="${OUT_ROOT}/workspace/MARCP/${task}_${ablation_short}"
    local weight_dir="${OUT_ROOT}/模型权重/MARCP/${task}/${ablation_short}"
    local metric_dir="${OUT_ROOT}/各个具体实验的指标/MARCP/${task}/${ablation_short}"
    mkdir -p "$base_out" "$weight_dir" "$metric_dir"
    (
        cd "${MODEL_12}"
        # CUDA_VISIBLE_DEVICES 必须在 Python 进程内 import torch 之前生效，故在 shell 中先 export，否则脚本内设置无效、进程可能无法正确上卡
        export CUDA_VISIBLE_DEVICES="$gpu_id"
        ${VENV_PYTHON:-python3} -u marcp_train.py \
            --lr 0.000785 \
            --CUDA "$gpu_id" \
            --signalMax 0 \
            --seqLength 101 \
            --bsz "${BATCH_SIZE}" \
            --encode OneHot \
            --modelFramework Treat \
            --epochs ${EPOCHS} \
            --max_stopping_step 9 \
            --input_channels 5 \
            --num_workers "${DATA_WORKERS_PER_JOB}" \
            --prefetch_factor "${PREFETCH_FACTOR:-4}" \
            --resource_fraction "${RESOURCE_FRAC_PER_JOB}" \
            -train "$train_5ch" \
            -validation "$val_5ch" \
            -test "$test_5ch" \
            --task "$task" \
            --lossM CombinedProfileLoss \
            --alpha 1.42 --beta 0.54 --gamma 1.80 \
            --weight_decay 0.000114 \
            --grad_clip 2.0 \
            --output "$base_out" \
            $ablation_flag \
            2>&1 | tee "${OUT_ROOT}/logs/MARCP_${task}_${ablation_short}.log" | grep -E --line-buffered 'Traceback|Error|Exception|WARNING|NameError|ModuleNotFoundError|CUDA out of memory' | while IFS= read -r line; do echo "$(date '+%Y-%m-%dT%H:%M:%S') [${task}_${ablation_short}] $line"; done >> "${OUT_ROOT}/run_ablation_errors.log" 2>/dev/null || true
    ); local train_ret=$?
    local actual_out="$base_out"
    [ "$ablation_short" != "baseline" ] && actual_out="$base_out/$ablation_short"
    # marcp_train 在仅有内存权重时也会写出 test_metrics_summary.txt；成功以指标文件为准，权重存在则一并拷贝
    if [ -f "${actual_out}/test_metrics_summary.txt" ]; then
        [ -f "${actual_out}/model_params.pth" ] && cp "${actual_out}/model_params.pth" "${weight_dir}/"
        for f in test_metrics_summary.txt test_metric.txt run_info.txt val_metric.txt train_metric.txt ablation_config.txt; do
            [ -f "${actual_out}/$f" ] && cp "${actual_out}/$f" "${metric_dir}/"
        done
        rm -rf "${base_out}"
        ( flock 10; echo "$task $ablation_short $(date +%s)" >> "$COMPLETED_LOG"; ) 10>"$COMPLETED_LOCK"
    else
        # 失败/异常退出时也清理 workspace 目录，避免 query 用 workspace 目录数误判“仍在运行”
        rm -rf "${base_out}" 2>/dev/null || true
        if [ "$train_ret" -ne 0 ]; then
            echo "[FAIL] $task $ablation_short: train_exit_code=$train_ret（详见 logs/MARCP_${task}_${ablation_short}.log 或 run_ablation_errors.log）"
        fi
    fi
    local pcc="N/A" train_min="N/A"
    [ -f "${metric_dir}/test_metrics_summary.txt" ] && pcc=$(grep -E "^Mean_PCC" "${metric_dir}/test_metrics_summary.txt" 2>/dev/null | awk -F'\t' '{print $2}')
    [ -f "${metric_dir}/run_info.txt" ] && train_min=$(grep -E "^training_time_min" "${metric_dir}/run_info.txt" 2>/dev/null | awk -F'\t' '{print $2}')
    ( flock 9; echo "| $task | $ablation_short | ${pcc:-N/A} | - | - | - | ${train_min:-N/A} | 完成 |" >> "$SUMMARY"; ) 9>"$SUMMARY_LOCK"
}

: > "$COMPLETED_LOG"
echo "start $(date +%s)" > "$STATE_FILE"
echo "total $total_to_run" >> "$STATE_FILE"
echo "n_gpu $N_GPU" >> "$STATE_FILE"
echo "total_slots $TOTAL_SLOTS" >> "$STATE_FILE"
echo "gpu_slots_desc GPU0最多${SLOTS_GPU0}、其它每卡${SLOTS_OTHER}(共${TOTAL_SLOTS}槽)" >> "$STATE_FILE"
echo "jobs" >> "$STATE_FILE"
cat "$TASK_LIST_FILE" >> "$STATE_FILE"

# 保留 EXIT 时打印退出码（不在此处删 PID/STATE，由上方正常流程处理）
trap 'echo "[全量消融] 主进程退出, exit_code=$?"' EXIT

# 按 GPU 槽位分配：初始化每卡当前任务数，任务结束时通过 .gpu_done 归还槽位
# 仅根据「本程序当前在各卡上的任务数」做最少任务数优先的轮询分配，不再依赖 nvidia-smi，
# 确保任务均匀落在 GPU0/1/2/3 上（例如 8/4/4/4），避免全部挤到 GPU0。
printf "0 0 0 0\n" > "$GPU_COUNTS_FILE"
: > "$GPU_DONE_FILE"
acquire_gpu_slot() {
    local c0 c1 c2 c3
    read -r c0 c1 c2 c3 < "$GPU_COUNTS_FILE" 2>/dev/null || { c0=0; c1=0; c2=0; c3=0; }
    # 在仍有空位的 GPU 中，选「当前任务数最少」的一张；并列时取编号最小者（避免旧逻辑在 c0=c1=… 时只能选 GPU0）
    local min_c=999999 g c cap
    for g in 0 1 2 3; do
        case "$g" in
            0) c=$c0; cap=$SLOTS_GPU0 ;;
            1) [ "$N_GPU" -lt 2 ] && continue; c=$c1; cap=$SLOTS_OTHER ;;
            2) [ "$N_GPU" -lt 3 ] && continue; c=$c2; cap=$SLOTS_OTHER ;;
            3) [ "$N_GPU" -lt 4 ] && continue; c=$c3; cap=$SLOTS_OTHER ;;
        esac
        [ "$c" -lt "$cap" ] || continue
        [ "$c" -lt "$min_c" ] && min_c=$c
    done
    local best_gpu=-1
    for g in 0 1 2 3; do
        case "$g" in
            0) c=$c0; cap=$SLOTS_GPU0 ;;
            1) [ "$N_GPU" -lt 2 ] && continue; c=$c1; cap=$SLOTS_OTHER ;;
            2) [ "$N_GPU" -lt 3 ] && continue; c=$c2; cap=$SLOTS_OTHER ;;
            3) [ "$N_GPU" -lt 4 ] && continue; c=$c3; cap=$SLOTS_OTHER ;;
        esac
        [ "$c" -lt "$cap" ] || continue
        [ "$c" -eq "$min_c" ] || continue
        best_gpu=$g
        break
    done
    # 若所有卡都达上限（极端情况），退化为 GPU0
    [ "$best_gpu" -lt 0 ] && best_gpu=0
    case "$best_gpu" in
        0) c0=$((c0+1)) ;;
        1) c1=$((c1+1)) ;;
        2) c2=$((c2+1)) ;;
        3) c3=$((c3+1)) ;;
    esac
    printf "%s %s %s %s\n" "$c0" "$c1" "$c2" "$c3" > "$GPU_COUNTS_FILE"
    echo "$best_gpu"
}
release_gpu_slot() {
    local done_gpu=$1 c0 c1 c2 c3
    read -r c0 c1 c2 c3 < "$GPU_COUNTS_FILE" 2>/dev/null || true
    case "$done_gpu" in
        0) c0=$((c0-1)); [ "$c0" -lt 0 ] && c0=0 ;;
        1) c1=$((c1-1)); [ "$c1" -lt 0 ] && c1=0 ;;
        2) c2=$((c2-1)); [ "$c2" -lt 0 ] && c2=0 ;;
        3) c3=$((c3-1)); [ "$c3" -lt 0 ] && c3=0 ;;
    esac
    printf "%s %s %s %s\n" "$c0" "$c1" "$c2" "$c3" > "$GPU_COUNTS_FILE"
}

# 非交互式 bash 下 jobs -p 常为 0，不可用；从 $$ 起 DFS 子进程树，统计 cmdline 含 marcp_train 的进程（实际训练数）
_walk_descendants_count_marcp() {
    local p="$1" child
    [ -r "/proc/$p/cmdline" ] || return 0
    if tr '\0' ' ' < "/proc/$p/cmdline" 2>/dev/null | grep -q "marcp_train"; then
        count_training_children__n=$((count_training_children__n + 1))
    fi
    for child in $(pgrep -P "$p" 2>/dev/null); do
        _walk_descendants_count_marcp "$child"
    done
}
count_training_children() {
    count_training_children__n=0
    local root
    for root in $(pgrep -P $$ 2>/dev/null); do
        _walk_descendants_count_marcp "$root"
    done
    echo "$count_training_children__n"
}

TASK_LINES=()
while IFS= read -r line; do [ -n "$line" ] && TASK_LINES+=("$line"); done < "$TASK_LIST_FILE"
total_tasks=${#TASK_LINES[@]}
next=0
# 已调度且尚未被 wait 回收的后台任务数（勿用「进程树里是否已有 marcp_train」代替：调度后 Python 尚未 exec 时 count=0，若据此跳过 wait 主脚本会立刻退出并 SIGHUP 杀光子任务）
inflight=0
_sched_first_logged=0
echo "常驻并行: 最多 $TOTAL_SLOTS 个任务，完成即补下一个 (共 $total_tasks 个消融任务，大任务优先)"

# 固定列表（如 batch11）：一次调度满本批（受 TOTAL_SLOTS 限制）后使用「wait」等待全部子进程。
# 避免仅依赖 wait -n + inflight 时在部分 bash/环境下出现竞态，主循环提前退出并 SIGHUP 杀死全部 marcp_train。
if [ "$_FIXED_LIST_ACTIVE" = "1" ]; then
    while [ "$next" -lt "$total_tasks" ] && [ "$inflight" -lt "$TOTAL_SLOTS" ]; do
        line="${TASK_LINES[$next]}"
        IFS='|' read -r task train_5ch val_5ch test_5ch ablation_short ablation_flag idx <<< "$line"
        gpu_id=$( ( flock 200; acquire_gpu_slot ) 200>"$GPU_COUNTS_LOCK" )
        [ -z "$gpu_id" ] && gpu_id=0
        base_out_abs="${OUT_ROOT}/workspace/MARCP/${task}_${ablation_short}"
        echo "========== [$idx/$total_tasks] $task | 消融: $ablation_short (GPU $gpu_id) =========="
        _run_one "$task" "$train_5ch" "$val_5ch" "$test_5ch" "$gpu_id" "$ablation_short" "$ablation_flag" "$idx" "$total_tasks" "$base_out_abs" &
        next=$((next + 1))
        inflight=$((inflight + 1))
    done
    echo "[调度] 固定列表模式: 已启动 $inflight 个后台任务，主进程将 wait 直至全部结束（防止提前退出导致训练被终止）"
    [ "$_sched_first_logged" = "0" ] && echo "[调度] 首波内层结束: next=$next inflight=$inflight total_tasks=$total_tasks TOTAL_SLOTS=$TOTAL_SLOTS"
    wait
    printf "0 0 0 0\n" > "$GPU_COUNTS_FILE"
    : > "$GPU_DONE_FILE"
else
while [ "$next" -lt "$total_tasks" ] || [ "$inflight" -gt 0 ]; do
    # 内层并发上界：用 inflight（已调度且尚未被 wait 回收）与 TOTAL_SLOTS 比较。
    # 勿用 count_training_children：子进程 cmdline 常含 marcp_train（含 DataLoader worker 等），会严重高估，误判已达槽位上限，表现为约 5 路即停。
    while [ "$next" -lt "$total_tasks" ] && [ "$inflight" -lt "$TOTAL_SLOTS" ]; do
        line="${TASK_LINES[$next]}"
        IFS='|' read -r task train_5ch val_5ch test_5ch ablation_short ablation_flag idx <<< "$line"
        gpu_id=$( ( flock 200; acquire_gpu_slot ) 200>"$GPU_COUNTS_LOCK" )
        [ -z "$gpu_id" ] && gpu_id=0
        base_out_abs="${OUT_ROOT}/workspace/MARCP/${task}_${ablation_short}"
        echo "========== [$idx/$total_tasks] $task | 消融: $ablation_short (GPU $gpu_id) =========="
        _run_one "$task" "$train_5ch" "$val_5ch" "$test_5ch" "$gpu_id" "$ablation_short" "$ablation_flag" "$idx" "$total_tasks" "$base_out_abs" &
        next=$((next + 1))
        inflight=$((inflight + 1))
    done
    if [ "$_sched_first_logged" = "0" ]; then
        echo "[调度] 首波内层结束: next=$next inflight=$inflight total_tasks=$total_tasks TOTAL_SLOTS=$TOTAL_SLOTS"
        if [ "$next" -lt "$total_tasks" ] || [ "$inflight" -lt "$total_tasks" ]; then
            echo "[调度] 警告: 首波未一次调度满本批（若 total_tasks 应为 11 且此处 next<11，多为服务器上脚本未更新为 inflight 并发，或任务列表行数不足）"
        fi
        _sched_first_logged=1
    fi
    if [ "$inflight" -gt 0 ]; then
        # 必须直接取 wait -n 的退出码：127=无可等待子进程。仅非 127 表示成功回收一个后台任务，inflight-- 并回收 GPU 槽位
        wait -n 2>/dev/null
        wait_ret=$?
        ct=$(count_training_children)
        if [ "$wait_ret" -eq 127 ]; then
            if [ "$inflight" -gt 0 ] && [ "$ct" -gt 0 ]; then
                echo "[全量消融] wait -n 返回 127 但仍有 marcp_train 进程，跳过本次槽位回收并重试"
                sleep 1
                continue
            fi
            if [ "$inflight" -gt 0 ] && [ "$ct" -eq 0 ]; then
                echo "[全量消融] wait -n 返回 127、inflight=${inflight} 但尚未检测到 marcp_train（启动中），跳过一次槽位回收"
                sleep 1
                continue
            fi
            continue
        fi
        if [ "$wait_ret" -lt 0 ] || [ "$wait_ret" -gt 255 ]; then
            echo "[全量消融] wait -n 异常码 $wait_ret，等待当前全部子进程后继续"
            wait 2>/dev/null || true
            continue
        fi
        inflight=$((inflight - 1))
        done_gpu=$( ( flock 201
          done_gpu=$(head -1 "$GPU_DONE_FILE" 2>/dev/null)
          [ -n "$done_gpu" ] && sed -i '1d' "$GPU_DONE_FILE" 2>/dev/null
          echo "$done_gpu"
        ) 201>"$GPU_DONE_LOCK" )
        if [ -n "$done_gpu" ]; then
            ( flock 200; release_gpu_slot "$done_gpu" ) 200>"$GPU_COUNTS_LOCK"
        fi
    fi
done
fi

# 仅当「本批全部成功」时才删 PID/STATE；否则保留 STATE 以便 query 显示本批总数与未完成数，并提示重新执行补跑
completed_count=0
[ -f "$COMPLETED_LOG" ] && completed_count=$(wc -l < "$COMPLETED_LOG")
if [ "$completed_count" -eq 0 ]; then
    echo "WARNING: 本轮 0 个任务成功完成，请检查 run_ablation.log 与 logs/ 下各任务 log。不删除 PID/STATE。"
    echo "completed_zero" > "${OUT_ROOT}/.run_ablation.done"
elif [ "$completed_count" -lt "$total_tasks" ]; then
    echo "本批 $total_tasks 个已全部调度完毕，其中 $completed_count 个成功；未成功 $((total_tasks - completed_count)) 个。保留 STATE，请重新执行: bash run_ablation_background.sh（将自动补跑未完成的）。"
    echo "incomplete $total_tasks $completed_count" >> "$STATE_FILE"
    rm -f "$PID_FILE"
    rm -f "${OUT_ROOT}/.run_ablation.done"
else
    rm -f "$PID_FILE" "$STATE_FILE"
    rm -f "${OUT_ROOT}/.run_ablation.done"
fi

echo "" >> "$SUMMARY"
echo "## 详细指标（从各实验 test_metrics_summary.txt / run_info.txt 汇总）" >> "$SUMMARY"
echo "" >> "$SUMMARY"
echo "| 任务 | 消融 | Mean_PCC | Spearman | Signal_AUC | Signal_AP | Peak_AUC | Peak_AP | 训练时长(min) |" >> "$SUMMARY"
echo "|------|------|----------|----------|------------|-----------|----------|---------|----------------|" >> "$SUMMARY"
for mdir in "${OUT_ROOT}/各个具体实验的指标/MARCP"/*/*/; do
    [ -d "$mdir" ] || continue
    task=$(basename "$(dirname "$mdir")")
    ab=$(basename "$mdir")
    f="${mdir}/test_metrics_summary.txt"
    r="${mdir}/run_info.txt"
    [ -f "$f" ] || continue
    pcc=$(grep -E "^Mean_PCC" "$f" 2>/dev/null | awk -F'\t' '{print $2}')
    spearman=$(grep -E "^Spearman" "$f" 2>/dev/null | awk -F'\t' '{print $2}')
    sig_auc=$(grep -E "^Signal_AUC" "$f" 2>/dev/null | awk -F'\t' '{print $2}')
    sig_ap=$(grep -E "^Signal_AP" "$f" 2>/dev/null | awk -F'\t' '{print $2}')
    peak_auc=$(grep -E "^Peak_AUC" "$f" 2>/dev/null | awk -F'\t' '{print $2}')
    peak_ap=$(grep -E "^Peak_AP" "$f" 2>/dev/null | awk -F'\t' '{print $2}')
    train_min="N/A"
    [ -f "$r" ] && train_min=$(grep -E "^training_time_min" "$r" 2>/dev/null | awk -F'\t' '{print $2}')
    echo "| $task | $ab | ${pcc:-N/A} | ${spearman:-N/A} | ${sig_auc:-N/A} | ${sig_ap:-N/A} | ${peak_auc:-N/A} | ${peak_ap:-N/A} | ${train_min:-N/A} |" >> "$SUMMARY"
done
echo "" >> "$SUMMARY"
echo "全量消融训练结束。总运行: $total_tasks 个。指标见 各个具体实验的指标/MARCP/，权重见 模型权重/MARCP/。"
echo "Done. 汇总: $SUMMARY"
