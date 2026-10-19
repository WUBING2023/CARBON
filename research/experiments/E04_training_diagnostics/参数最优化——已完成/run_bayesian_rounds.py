"""
分阶段贝叶斯优化脚本

优化流程设计：
- 第1阶段（已完成）：25个trial随机探索（热身阶段，让TPE采样器积累先验知识）
- 第2-4阶段：每轮10个trial，基于历史数据进行贝叶斯优化
- 总共：25 + 10*3 = 55个trial

为什么分阶段：
1. 第1阶段的25个并行trial主要用于随机探索超参数空间
2. 后续阶段TPE采样器会基于已完成的trial建模，智能选择更有希望的超参数组合
3. 每轮10个trial同时运行，轮次间TPE有机会充分利用新数据更新模型
"""
import os
import sys
import time
import argparse
import subprocess

def get_completed_trials():
    """获取已完成的trial数量"""
    from optuna.storages import JournalStorage
    try:
        from optuna.storages.journal import JournalFileBackend
    except ImportError:
        from optuna.storages import JournalFileStorage as JournalFileBackend
    import optuna
    
    journal_path = "/data/CARBON_PROJECT/参数最优化/results/optuna_journal.log"
    if not os.path.exists(journal_path):
        return 0
    
    storage = JournalStorage(JournalFileBackend(journal_path))
    try:
        study = optuna.load_study(
            study_name="marcp_hyperparam_optimization",
            storage=storage,
        )
        completed = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
        return len(completed)
    except Exception as e:
        print(f"读取study失败: {e}")
        return 0


def run_optimization_round(target_trials: int, cuda_config: str, round_num: int):
    """
    运行一轮优化
    
    Args:
        target_trials: 目标完成的总trial数（不是本轮新增数）
        cuda_config: GPU配置字符串
        round_num: 当前轮次编号
    """
    current = get_completed_trials()
    new_trials = target_trials - current
    
    print(f"\n{'='*60}")
    print(f"第 {round_num} 轮贝叶斯优化")
    print(f"  当前已完成: {current} trials")
    print(f"  目标完成: {target_trials} trials")
    print(f"  本轮新增: {new_trials} trials")
    print(f"  GPU配置: {cuda_config}")
    print(f"{'='*60}\n")
    
    if new_trials <= 0:
        print(f"已达到目标，跳过本轮")
        return 0
    
    cmd = [
        "python3", "run_optimization_parallel.py",
        "--n_trials", str(target_trials),
        "--n_datasets", "5",
        "--cuda", cuda_config,
    ]
    
    process = subprocess.Popen(
        cmd,
        cwd="/data/CARBON_PROJECT/参数最优化",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    
    for line in process.stdout:
        print(line, end="", flush=True)
    
    process.wait()
    return process.returncode


def print_best_result():
    """打印当前最佳结果"""
    from optuna.storages import JournalStorage
    try:
        from optuna.storages.journal import JournalFileBackend
    except ImportError:
        from optuna.storages import JournalFileStorage as JournalFileBackend
    import optuna
    
    journal_path = "/data/CARBON_PROJECT/参数最优化/results/optuna_journal.log"
    storage = JournalStorage(JournalFileBackend(journal_path))
    study = optuna.load_study(
        study_name="marcp_hyperparam_optimization",
        storage=storage,
    )
    
    print(f"\n{'='*60}")
    print("当前最佳结果")
    print(f"{'='*60}")
    print(f"最佳Trial: {study.best_trial.number}")
    print(f"最佳值 (mean_val_pcc): {study.best_value:.4f}")
    print("\n最佳超参数:")
    for k, v in study.best_params.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.6g}")
        else:
            print(f"  {k}: {v}")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(description="分阶段贝叶斯优化")
    parser.add_argument("--initial_trials", type=int, default=25,
                       help="第1阶段已完成的trial数（默认25）")
    parser.add_argument("--rounds", type=int, default=3,
                       help="额外优化轮数（默认3）")
    parser.add_argument("--trials_per_round", type=int, default=10,
                       help="每轮trial数量（默认10）")
    parser.add_argument("--cuda", type=str, default="0:4,1:2,2:2,3:2",
                       help="GPU配置（默认0:4,1:2,2:2,3:2）")
    args = parser.parse_args()
    
    total_trials = args.initial_trials + args.rounds * args.trials_per_round
    
    print("="*60)
    print("分阶段贝叶斯超参数优化")
    print("="*60)
    print(f"优化配置:")
    print(f"  第1阶段（随机探索）: {args.initial_trials} trials（已完成）")
    print(f"  后续轮数: {args.rounds}")
    print(f"  每轮trials: {args.trials_per_round}")
    print(f"  总trials: {total_trials}")
    print(f"  GPU配置: {args.cuda}")
    print("="*60)
    
    current_trials = get_completed_trials()
    print(f"\n当前已完成 {current_trials} 个trial")
    
    if current_trials > 0:
        print_best_result()
    
    if current_trials < args.initial_trials:
        print(f"\n警告: 第1阶段尚未完成 ({current_trials}/{args.initial_trials})")
        print("请先完成第1阶段的随机探索")
        return
    
    for round_idx in range(args.rounds):
        round_num = round_idx + 2
        target = args.initial_trials + (round_idx + 1) * args.trials_per_round
        
        print(f"\n{'#'*60}")
        print(f"开始第 {round_num} 轮 (总目标: {target} trials)")
        print(f"{'#'*60}")
        
        run_optimization_round(
            target_trials=target,
            cuda_config=args.cuda,
            round_num=round_num,
        )
        
        completed = get_completed_trials()
        print(f"\n第 {round_num} 轮完成! 已完成 {completed} trials")
        print_best_result()
        
        if round_idx < args.rounds - 1:
            print(f"等待3秒后开始下一轮...")
            time.sleep(3)
    
    print("\n" + "="*60)
    print("所有优化轮次完成!")
    print("="*60)
    print_best_result()
    
    final_trials = get_completed_trials()
    print(f"\n总共完成 {final_trials} 个trial")
    print(f"预期: {total_trials} trials")


if __name__ == "__main__":
    main()
