#!/bin/bash
# =============================================================================
# MARCP 全量消融：查询训练进度（状态、本轮完成、当前在跑任务、ETA、GPU/CPU 资源）
# 用法：bash query_ablation_progress.sh
# 输出目录：模型权重/MARCP、各个具体实验的指标/MARCP、workspace/MARCP
# =============================================================================
BASE="/data/CARBON_PROJECT"
OUT_ROOT="${BASE}/消融实验/全量消融运行"
DATA_ROOT="${BASE}/数据/h5data"
PID_FILE="${OUT_ROOT}/.run_ablation.pid"
STATE_FILE="${OUT_ROOT}/.run_ablation.state"
COMPLETED_LOG="${OUT_ROOT}/.run_ablation.completed"
WEIGHT_BASE="${OUT_ROOT}/模型权重/MARCP"
METRIC_BASE="${OUT_ROOT}/各个具体实验的指标/MARCP"
WORKSPACE="${OUT_ROOT}/workspace/MARCP"

echo "========== 全量消融训练进度 =========="
# 判断某 PID 是否属于「本程序」：cmdline 含「全量消融运行」或（主调度在跑且该进程为主调度的后代，用 /proc 父链追溯）
is_our_pid() {
    local c=$1
    local mp=$2
    [ -z "$c" ] || [ -z "$mp" ] && return 1
    if [ -r "/proc/$c/cmdline" ]; then
        cmd=$(tr '\0' ' ' < "/proc/$c/cmdline" 2>/dev/null)
        echo "$cmd" | grep -q "全量消融运行" && return 0
    fi
    # 沿父进程链向上，若到达主调度 PID 则视为本程序（/proc/pid/stat 格式: pid (comm) state ppid ...）
    local p="$c"
    while [ -n "$p" ] && [ "$p" -gt 0 ] 2>/dev/null; do
        [ "$p" = "$mp" ] && return 0
        [ -r "/proc/$p/stat" ] || break
        p=$(sed 's/^[0-9]* ([^)]*) //' "/proc/$p/stat" 2>/dev/null | awk '{print $2}')  # 去掉 "pid (comm) " 后取 ppid
        [ -z "$p" ] && break
    done
    return 1
}
our_proc_count=0
main_pid=""
[ -f "$PID_FILE" ] && main_pid=$(cat "$PID_FILE")
for c in $(pgrep -f "marcp_train" 2>/dev/null); do
    is_our_pid "$c" "$main_pid" && our_proc_count=$((our_proc_count + 1))
done

if [ -f "$PID_FILE" ]; then
    pid=$(cat "$PID_FILE")
    if kill -0 "$pid" 2>/dev/null; then
        echo "状态: 运行中 (主调度 PID=$pid，后台由 run_ablation_background.sh 启动)"
    else
        echo "状态: 主调度已退出 (PID=$pid 已结束)"
        [ "$our_proc_count" -gt 0 ] && echo "        当前仍有 $our_proc_count 个本程序训练进程在跑；跑完后不会自动接新任务，需重新运行 run_ablation_background.sh 继续排队"
        if [ -f "${OUT_ROOT}/.run_ablation.done" ] && grep -q "completed_zero" "${OUT_ROOT}/.run_ablation.done" 2>/dev/null; then
            echo "        上轮 0 个任务成功完成，请检查 run_ablation.log 与 logs/ 下各任务 log 排查原因"
        fi
    fi
else
    echo "状态: 未运行（未检测到本次运行的 PID 文件）"
    [ "$our_proc_count" -gt 0 ] && echo "        但检测到 $our_proc_count 个本程序训练进程在跑，可能为主调度退出后残留；若需继续排队请重新运行 run_ablation_background.sh"
    # 无 PID 但存在大量 workspace 目录时，说明曾跑过消融且主调度已退出
    if [ -d "$WORKSPACE" ]; then
        n_ws=$(find "$WORKSPACE" -maxdepth 1 -type d ! -path "$WORKSPACE" 2>/dev/null | wc -l)
        [ "$n_ws" -gt 0 ] && echo "        检测到 $n_ws 个 workspace 任务目录；若需继续排队请执行: bash run_ablation_background.sh（请始终用该方式启动以免关终端后状态丢失）"
    fi
fi
# 下面会读 STATE 得到 batch_total；若未运行且本批有未完成，在输出完「本批未完成」后提示补跑（见下方）

# 本批任务数、本轮已完成、槽位信息
batch_total=""
start_ts=""
n_gpu=""
total_slots=""
gpu_slots_desc=""
if [ -f "$STATE_FILE" ]; then
    in_jobs=0
    while read -r line; do
        case "$line" in
            start*) start_ts="${line#start }" ;;
            total_slots*) total_slots="${line#total_slots }" ;;
            total*) batch_total="${line#total }" ;;
            n_gpu*) n_gpu="${line#n_gpu }" ;;
            gpu_slots_desc*) gpu_slots_desc="${line#gpu_slots_desc }" ;;
            jobs) in_jobs=1 ;;
            *)
                if [ "$in_jobs" = "1" ] && [ -n "$line" ]; then
                    : # 可用来统计本批任务
                fi
                ;;
        esac
    done < "$STATE_FILE"
fi
batch_total=$((batch_total + 0))
this_run_completed=0
[ -f "$COMPLETED_LOG" ] && this_run_completed=$(wc -l < "$COMPLETED_LOG")
pending_this_run=0
[ "$batch_total" -gt 0 ] && pending_this_run=$((batch_total - this_run_completed))
[ "$pending_this_run" -lt 0 ] && pending_this_run=0

echo "本批消融任务数: $batch_total"
# 运行中但本批为 0：多为刚启动尚未写入 STATE，提示稍后重试
if [ "$batch_total" -eq 0 ] && [ -f "$PID_FILE" ]; then
    p=$(cat "$PID_FILE"); if kill -0 "$p" 2>/dev/null; then
        echo "        （若刚启动，本批任务数可能尚未写入，请几秒后重试: bash query_ablation_progress.sh）"
    fi
fi
echo "本轮已完成: $this_run_completed"
echo "本批未完成（尚在跑/队列）: $pending_this_run"
[ "$batch_total" -gt 0 ] && [ "$pending_this_run" -gt 0 ] && [ ! -f "$PID_FILE" ] && echo "        提示: 本批已调度完毕，尚有 $pending_this_run 个未完成，请执行 bash run_ablation_background.sh 继续补跑。"
[ -n "$n_gpu" ] && echo "当前 GPU 数: $n_gpu"
[ -n "$total_slots" ] && echo "总并发槽位: $total_slots"
[ -n "$gpu_slots_desc" ] && echo "槽位说明: $gpu_slots_desc"

# 当前在跑哪些消融：以 marcp_train 进程为准（从 --output 解析 任务_消融）；workspace 目录在失败/成功 rm 后会变少，不能单独代表并发数
echo ""
echo "--- 当前正在运行的消融实验（任务_消融）---"
echo "  本程序 marcp_train 进程数（真实并发）: $our_proc_count"
echo "  （下列名单从各进程命令行 --output 目录名解析；仅数 workspace 子目录会误判，因任务结束会 rm 临时目录）"
SHOW_RUNNING_HEAD=8
RUNNING_TMP=$(mktemp)
if [ "$our_proc_count" -gt 0 ] && [ -n "$main_pid" ]; then
    for c in $(pgrep -f "marcp_train" 2>/dev/null); do
        is_our_pid "$c" "$main_pid" || continue
        [ -r "/proc/$c/cmdline" ] || continue
        cmd=$(tr '\0' ' ' < "/proc/$c/cmdline" 2>/dev/null)
        out=$(echo "$cmd" | sed 's/.*--output //;s/ .*//')
        if [ -n "$out" ]; then
            bn=$(basename "$out" 2>/dev/null)
            [ -n "$bn" ] && echo "$bn" >> "$RUNNING_TMP"
        fi
    done
fi
if [ -s "$RUNNING_TMP" ]; then
    nproc_run=$(sort -u "$RUNNING_TMP" | wc -l)
    sort -u "$RUNNING_TMP" | head -n "$SHOW_RUNNING_HEAD" | while IFS= read -r line; do
        [ -n "$line" ] && echo "  $line"
    done
    if [ "$nproc_run" -gt "$SHOW_RUNNING_HEAD" ]; then
        echo "  ... 等共 $nproc_run 个（按进程解析）"
    fi
else
    if [ "$our_proc_count" -gt 0 ]; then
        echo "  (进程在跑但未能从 --output 解析目录名，可能正在启动)"
    else
        echo "  (无)"
    fi
fi
rm -f "$RUNNING_TMP"
if [ -d "$WORKSPACE" ]; then
    ws_n=$(find "$WORKSPACE" -maxdepth 1 -type d ! -path "$WORKSPACE" 2>/dev/null | wc -l)
    echo "  workspace 下尚未删除的临时目录数（仅作参考，常小于真实并发）: $ws_n"
fi
if [ -f "${OUT_ROOT}/run_ablation.log" ]; then
    _sched=$(grep '^\[调度\] 首波内层结束:' "${OUT_ROOT}/run_ablation.log" 2>/dev/null | tail -n 1)
    if [ -n "$_sched" ]; then
        echo ""
        echo "--- 调度首波诊断（是否一次调度满本批）---"
        echo "  $_sched"
        echo "  解读: 若 next 与 total_tasks 相等 → 首波已同时调度满本批，并非「只开 5 路」。"
        echo "  若本批 marcp 进程数 < 本批任务数: 常见原因是部分任务已训练结束（成功后脚本会删除 workspace，进程退出）。"
        echo "        小数据集每 epoch 步数极少、早停快时，单任务可能数秒～数分钟就跑完（如 SERBP1 与大体量任务混跑时，会先只剩大任务在跑）。"
        echo "        只有当下方出现 [FAIL] 或 logs/MARCP_*.log 有 Traceback 时，才按「失败」排查。"
    fi
    _fl=$(grep '\[FAIL\]' "${OUT_ROOT}/run_ablation.log" 2>/dev/null | tail -n 8)
    if [ -n "$_fl" ]; then
        echo ""
        echo "--- 最近失败任务（run_ablation.log 中 [FAIL]）---"
        echo "$_fl" | sed 's/^/  /'
    fi
fi
if [ "$batch_total" -gt 0 ] && [ "$our_proc_count" -gt 0 ] && [ "$our_proc_count" -lt "$batch_total" ] && [ -z "$_sched" ]; then
    echo "  说明: 本批共 $batch_total 个任务，当前仅 $our_proc_count 个 marcp_train；未找到 [调度] 首波行时，可 tail run_ablation.log 确认是否一次调度满批。"
fi

# 正式目录已完成数（按 任务/消融 计）
completed_count=0
[ -d "$METRIC_BASE" ] && completed_count=$(find "$METRIC_BASE" -mindepth 2 -maxdepth 3 -name "test_metrics_summary.txt" 2>/dev/null | wc -l)
echo ""
echo "正式目录累计已完成: $completed_count 个 (任务×消融)"

# 时间估算
echo ""
echo "--- 时间估算 ---"
if [ -n "$start_ts" ] && [ "$start_ts" -gt 0 ]; then
    now=$(date +%s)
    total_elapsed_sec=$((now - start_ts))
    total_elapsed_min=$((total_elapsed_sec / 60))
    if [ "$total_elapsed_sec" -lt 0 ]; then
        echo "总已运行: 未记录或时间异常"
    else
        echo "总已运行: $total_elapsed_min 分钟 ($((total_elapsed_sec / 3600)) 小时)"
    fi
    if [ "$this_run_completed" -gt 0 ] && [ "$total_elapsed_sec" -ge 60 ]; then
        avg=$((total_elapsed_sec / this_run_completed))
        echo "平均每任务约: $((avg / 60)) 分钟"
        if [ "$pending_this_run" -gt 0 ] && [ "$avg" -gt 0 ]; then
            remaining_sec=$((pending_this_run * avg))
            echo "预计剩余时间: 约 $((remaining_sec / 60)) 分钟 ($((remaining_sec / 3600)) 小时)"
            echo "预计总完成时间: $(date -d "@$((now + remaining_sec))" '+%Y-%m-%d %H:%M' 2>/dev/null || date -r $((now + remaining_sec)) '+%Y-%m-%d %H:%M' 2>/dev/null)"
        fi
    elif [ "$this_run_completed" -eq 0 ] && [ -f "$PID_FILE" ] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
        echo "平均每任务: 需至少完成 1 个后显示"
    fi
else
    echo "总已运行: 未记录（请通过 run_ablation_background.sh 启动）"
fi

# 资源
echo ""
echo "--- 资源使用 ---"
# CPU 使用情况
if command -v top >/dev/null 2>&1; then
    cpu_used=$(top -bn1 2>/dev/null | grep "Cpu(s)" | sed 's/.*, *\([0-9.]*\)%* id.*/\1/' | awk '{printf "%.1f", 100-$1}')
    [ -n "$cpu_used" ] && echo "CPU 使用率: ${cpu_used}%"
else
    echo "CPU 使用率: (无 top，可用 mpstat 或查看 /proc/stat)"
fi
if [ -r /proc/loadavg ]; then
    read -r load1 load5 load15 _ _ < /proc/loadavg 2>/dev/null
    [ -n "$load1" ] && echo "CPU 负载 (1/5/15 min): $load1 $load5 $load15"
fi
if [ -r /proc/cpuinfo ]; then
    cpu_cores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null)
    [ -n "$cpu_cores" ] && echo "CPU 核数: $cpu_cores"
fi
# 内存使用情况
if command -v free >/dev/null 2>&1; then
    mem_info=$(free -m 2>/dev/null | awk '/^Mem:/ {printf "%.1f%% (%d / %d MB)", ($3/$2)*100, $3, $2}')
    [ -n "$mem_info" ] && echo "内存使用率: $mem_info"
    swap_info=$(free -m 2>/dev/null | awk '/^Swap:/ {if($2>0) printf "%.1f%% (%d / %d MB)", ($3/$2)*100, $3, $2; else print "未启用"}')
    [ -n "$swap_info" ] && echo "Swap 使用率: $swap_info"
elif [ -r /proc/meminfo ]; then
    mem_total=$(awk '/^MemTotal:/ {print int($2/1024)}' /proc/meminfo 2>/dev/null)
    mem_avail=$(awk '/^MemAvailable:/ {print int($2/1024)}' /proc/meminfo 2>/dev/null)
    if [ -n "$mem_total" ] && [ -n "$mem_avail" ] && [ "$mem_total" -gt 0 ]; then
        mem_used=$((mem_total - mem_avail))
        mem_pct=$(awk "BEGIN {printf \"%.1f\", ($mem_used/$mem_total)*100}")
        echo "内存使用率: ${mem_pct}% ($mem_used / $mem_total MB)"
    fi
fi
if command -v nvidia-smi >/dev/null 2>&1; then
    echo "GPU:"
    # 本程序 = 全量消融运行的训练进程（cmdline 含「全量消融运行」或父链到主调度，或主调度在跑时其后代中 cmdline 含 marcp_train）
    OUR_PIDS_FILE=$(mktemp)
    main_pid_gpu=""
    [ -f "$PID_FILE" ] && main_pid_gpu=$(cat "$PID_FILE")
    for c in $(pgrep -f "marcp_train" 2>/dev/null); do
        is_our_pid "$c" "$main_pid_gpu" && echo "$c" >> "$OUR_PIDS_FILE"
    done
    # 主调度在跑时，用 ps 建树收集主调度所有后代中 cmdline 含 marcp_train 的进程（比递归 pgrep 更稳）
    if [ -n "$main_pid_gpu" ] && kill -0 "$main_pid_gpu" 2>/dev/null; then
        while read -r c; do
            [ -r "/proc/$c/cmdline" ] || continue
            tr '\0' ' ' < "/proc/$c/cmdline" 2>/dev/null | grep -q "marcp_train" && echo "$c" >> "$OUR_PIDS_FILE"
        done < <(ps -eo pid,ppid --no-headers 2>/dev/null | awk -v root="$main_pid_gpu" '
            { parent[$1]=$2 }
            END {
                for (p in parent) {
                    q=p
                    while (q != "") {
                        if (q+0 == root+0) { desc[p]=1; break }
                        if (!(q in parent)) break
                        q=parent[q]
                    }
                }
                for (d in desc) print d
            }')
    fi
    # 补漏：nvidia-smi 报告的 GPU 进程可能包含上轮遗留等；cmdline 含 marcp_train 且含「全量消融运行」或「全量消融」（路径可能被截断）即算本程序
    for gpu_idx in $(nvidia-smi --query-gpu=index --format=csv,noheader,nounits 2>/dev/null); do
        while IFS=',' read -r c _; do
            c=$(echo "$c" | tr -d ' ')
            [ -z "$c" ] || [ "$c" = "pid" ] && continue
            grep -q "^${c}$" "$OUR_PIDS_FILE" 2>/dev/null && continue
            [ -r "/proc/$c/cmdline" ] || continue
            cmd=$(tr '\0' ' ' < "/proc/$c/cmdline" 2>/dev/null)
            echo "$cmd" | grep -q "marcp_train" || continue
            # 全量消融运行 / 全量消融 / 消融实验（输出路径含任一词即视为本程序，排除 论文data/全量运行）
            (echo "$cmd" | grep -q "全量消融运行") || (echo "$cmd" | grep -q "全量消融") || (echo "$cmd" | grep -q "消融实验") || continue
            echo "$cmd" | grep -q "论文data/全量运行" && continue
            echo "$c" >> "$OUR_PIDS_FILE"
        done < <(nvidia-smi -i "$gpu_idx" --query-compute-apps=pid,used_memory --format=csv,noheader,nounits 2>/dev/null)
    done
    sort -u "$OUR_PIDS_FILE" 2>/dev/null > "${OUR_PIDS_FILE}.tmp" && mv "${OUR_PIDS_FILE}.tmp" "$OUR_PIDS_FILE"
    n_our_pids=0
    [ -s "$OUR_PIDS_FILE" ] && n_our_pids=$(wc -l < "$OUR_PIDS_FILE")
    # 占用 GPU 的才是训练主进程（其余多为 DataLoader 子进程：9 任务 × num_workers ≈ 63）
    n_our_gpu_pids=0
    for pid in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits 2>/dev/null | tr -d ' '); do
        grep -q "^${pid}$" "$OUR_PIDS_FILE" 2>/dev/null && n_our_gpu_pids=$((n_our_gpu_pids + 1))
    done
    [ "$n_our_pids" -gt 0 ] && echo "（本程序进程数: $n_our_pids，其中 $n_our_gpu_pids 个占用 GPU；其余为 DataLoader 子进程）"
    # 本程序 0% 时提示：可能无本程序进程在跑、或本程序尚未占用 GPU（如主调度已退出或其它任务占满 GPU）
    if [ "$n_our_pids" -eq 0 ]; then
        main_running=0
        [ -f "$PID_FILE" ] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null && main_running=1
        if [ "$main_running" -eq 0 ]; then
            echo "（说明: 未检测到本程序训练进程；若主调度已退出则本程序 0% 为正常，当前 GPU 占用可能来自其它任务。详见 说明文档和数据记录/为何会显示未运行.md 第六节）"
        fi
    fi
    # 将本程序 PID 列表放入变量，避免管道子 shell 内读文件异常；空格包围便于 [[ *" $pid "* ]] 匹配
    OUR_PIDS_SPACE=""
    [ -s "$OUR_PIDS_FILE" ] && OUR_PIDS_SPACE=" $(tr '\n' ' ' < "$OUR_PIDS_FILE" | tr -s ' ') "
    # 各 GPU 上本程序计算占比（pmon 的 sm 列）与显存：用进程替换让 while 在主 shell 运行并直接用 OUR_PIDS_SPACE
    OUR_SM_FILE=$(mktemp)
    nvidia-smi pmon -c 1 -s u 2>/dev/null | awk -v pf="$OUR_PIDS_FILE" '
        BEGIN { while((getline < pf)>0) { gsub(/\r$/,""); gsub(/\n$/,""); pids[$0+0]=1 }; close(pf) }
        /^[[:space:]]*[0-9]/ { gpu=$1; pid=$2+0; sm=$4; if (pid in pids && sm != "-") sum[gpu]+=sm }
        END { for (g in sum) print g, sum[g] }
    ' > "$OUR_SM_FILE"
    total_our_sm=0
    total_our_mem=0
    while IFS=',' read -r idx name util_gpu util_mem mem_used mem_total; do
        idx=$(echo "$idx" | sed 's/^ *//;s/ *$//')
        name=$(echo "$name" | sed 's/^ *//;s/ *$//')
        util_gpu=$(echo "$util_gpu" | sed 's/^ *//;s/ *$//')
        util_mem=$(echo "$util_mem" | sed 's/^ *//;s/ *$//')
        mem_used=$(echo "$mem_used" | sed 's/^ *//;s/ *$//')
        mem_total=$(echo "$mem_total" | sed 's/^ *//;s/ *$//')
        our_sm=$(awk -v i="$idx" '$1==i {print $2}' "$OUR_SM_FILE" 2>/dev/null)
        [ -z "$our_sm" ] && our_sm="0"
        our_mem=0
        while IFS=',' read -r pid used_mb; do
            pid=$(echo "$pid" | tr -d ' ')
            used_mb=$(echo "$used_mb" | tr -d ' ')
            [ -z "$used_mb" ] || [ "$used_mb" = "memory.used [MiB]" ] && continue
            [[ " $OUR_PIDS_SPACE " == *" $pid "* ]] && our_mem=$((our_mem + used_mb))
        done < <(nvidia-smi -i "$idx" --query-compute-apps=pid,used_memory --format=csv,noheader,nounits 2>/dev/null)
        total_our_sm=$((total_our_sm + our_sm))
        total_our_mem=$((total_our_mem + our_mem))
        our_pct="0"
        [ "$mem_total" -gt 0 ] 2>/dev/null && our_pct=$(awk "BEGIN {printf \"%.1f\", ($our_mem/$mem_total)*100}")
        note=""
        [ "$our_mem" -gt "$mem_used" ] 2>/dev/null && note=" (进程合计>整卡used，统计口径不同)"
        echo "  GPU $idx ($name): 计算 ${util_gpu}%, 显存 ${util_mem}% (${mem_used} / ${mem_total} MiB)${note}"
    done < <(nvidia-smi --query-gpu=index,name,utilization.gpu,utilization.memory,memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null)
    # 已识别本程序进程但各 GPU 上本程序占用均为 0：说明这些 PID 未出现在 nvidia-smi 的 GPU 进程列表中
    if [ "$n_our_pids" -gt 0 ] && [ "$total_our_sm" -eq 0 ] && [ "$total_our_mem" -eq 0 ]; then
        gpu_pids=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits 2>/dev/null | tr -d ' ' | sort -u | tr '\n' ',' | sed 's/,$//')
        [ -n "$gpu_pids" ] && echo "（说明: 已识别 $n_our_pids 个本程序进程，但当前无一出现在各 GPU 占用列表中 → 本程序计算/显存为 0。当前占用 GPU 的 PID: $gpu_pids；可用 ps -p <PID> -o args= 核对是否为本程序。若均为其它任务，说明本程序未上卡，可改用预留 GPU 启动，例如: ABLATION_GPU_IDS=2,3 CUDA_VISIBLE_DEVICES=2,3 bash run_ablation_background.sh）"
        [ -z "$gpu_pids" ] && echo "（说明: 已识别 $n_our_pids 个本程序进程，但当前无一出现在各 GPU 占用列表中；可稍后再次查询或运行 nvidia-smi 核对）"
    fi
    rm -f "$OUR_PIDS_FILE" "$OUR_SM_FILE"
else
    echo "GPU: (未检测到 nvidia-smi)"
fi
echo "===================================="
