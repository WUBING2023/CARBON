"""多GPU并行超参数优化脚本"""
import os
import sys
import json
import argparse
import subprocess
import multiprocessing
from typing import Dict, List, Any, Tuple

import optuna
from optuna.samplers import TPESampler
from optuna.pruners import MedianPruner
from optuna.storages import JournalStorage

try:
    from optuna.storages.journal import JournalFileBackend
except ImportError:
    from optuna.storages import JournalFileStorage as JournalFileBackend

from config import (
    MARCP_DIR, OUTPUT_ROOT, SEARCH_SPACE, FIXED_PARAMS, build_train_command
)
from dataset_selector import get_selected_datasets


class MaxTrialsCallback:
    """在达到最大trials数时停止优化"""
    def __init__(self, max_trials: int):
        self.max_trials = max_trials
    
    def __call__(self, study: optuna.Study, trial: optuna.trial.FrozenTrial):
        n_complete = len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])
        if n_complete >= self.max_trials:
            study.stop()


def sample_hyperparams(trial: optuna.Trial) -> Dict[str, Any]:
    """从搜索空间采样超参数"""
    params = {}
    for name, spec in SEARCH_SPACE.items():
        if spec["type"] == "float":
            params[name] = trial.suggest_float(
                name, spec["low"], spec["high"], log=spec.get("log", False)
            )
        elif spec["type"] == "int":
            params[name] = trial.suggest_int(name, spec["low"], spec["high"])
        elif spec["type"] == "categorical":
            params[name] = trial.suggest_categorical(name, spec["choices"])
    return params


def parse_validation_metrics(output_dir: str) -> Dict[str, float]:
    """解析验证集指标"""
    metrics = {
        "val_pcc": float("nan"),
        "val_prauc": float("nan"),
        "val_spearman": float("nan"),
        "val_loss": float("nan"),
    }
    val_metric_file = os.path.join(output_dir, "val_metric.txt")
    if os.path.isfile(val_metric_file):
        try:
            import pandas as pd
            df = pd.read_csv(val_metric_file, sep="\t")
            if len(df) > 0:
                last_row = df.iloc[-1]
                metrics["val_pcc"] = float(last_row.get("pearson", float("nan")))
                metrics["val_prauc"] = float(last_row.get("PR_AUC", float("nan")))
                metrics["val_spearman"] = float(last_row.get("spearman", float("nan")))
                metrics["val_loss"] = float(last_row.get("loss", float("nan")))
        except Exception:
            pass
    return metrics


def run_single_training(
    task_name: str,
    params: Dict[str, Any],
    output_dir: str,
    cuda_device: str = "0",
    timeout: int = 7200,
) -> Tuple[Dict[str, float], bool]:
    """运行单次训练"""
    os.makedirs(output_dir, exist_ok=True)
    cmd = build_train_command(task_name, params, output_dir, cuda_device)
    
    cmd_file = os.path.join(output_dir, "train_command.sh")
    with open(cmd_file, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(cmd + "\n")
    
    log_file = os.path.join(output_dir, "train.log")
    
    try:
        with open(log_file, "w") as log_f:
            process = subprocess.Popen(
                cmd.replace("\\\n  ", " "),
                shell=True,
                cwd=MARCP_DIR,
                stdout=log_f,
                stderr=subprocess.STDOUT,
            )
            process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        process.kill()
        return {}, False
    except Exception:
        return {}, False
    
    metrics = parse_validation_metrics(output_dir)
    val_metric_file = os.path.join(output_dir, "val_metric.txt")
    if os.path.isfile(val_metric_file):
        return metrics, True
    return {}, False


def run_multi_dataset_evaluation(
    params: Dict[str, Any],
    datasets: List[str],
    trial_dir: str,
    cuda_device: str = "0",
) -> Dict[str, float]:
    """在多个数据集上评估超参数"""
    all_metrics = []
    successful = 0
    
    for task_name in datasets:
        dataset_dir = os.path.join(trial_dir, task_name)
        metrics, success = run_single_training(
            task_name, params, dataset_dir, cuda_device
        )
        if success and not any(
            str(v) == "nan" for k, v in metrics.items() 
            if k in ["val_pcc", "val_prauc"]
        ):
            all_metrics.append(metrics)
            successful += 1
    
    import numpy as np
    if not all_metrics:
        return {
            "mean_val_pcc": float("nan"),
            "mean_val_prauc": float("nan"),
            "mean_val_spearman": float("nan"),
            "std_val_pcc": float("nan"),
            "successful_datasets": 0,
            "total_datasets": len(datasets),
        }

    aggregated = {
        "mean_val_pcc": np.nanmean([m["val_pcc"] for m in all_metrics]),
        "mean_val_prauc": np.nanmean([m["val_prauc"] for m in all_metrics]),
        "mean_val_spearman": np.nanmean([m["val_spearman"] for m in all_metrics]),
        "std_val_pcc": np.nanstd([m["val_pcc"] for m in all_metrics]),
        "successful_datasets": successful,
        "total_datasets": len(datasets),
    }
    return aggregated


def objective(trial: optuna.Trial, datasets: List[str], cuda_device: str) -> float:
    """优化目标函数"""
    params = sample_hyperparams(trial)
    
    trial_dir = os.path.join(OUTPUT_ROOT, "trial_outputs", f"trial_{trial.number:04d}")
    os.makedirs(trial_dir, exist_ok=True)
    
    with open(os.path.join(trial_dir, "hyperparameters.json"), "w") as f:
        json.dump(params, f, indent=2)
    
    print(f"[Worker | GPU {cuda_device}] Trial {trial.number}", flush=True)
    print(f"  超参数: lr={params['lr']:.5f}, bsz={params['bsz']}, alpha={params['alpha']:.2f}", flush=True)
    
    results = run_multi_dataset_evaluation(params, datasets, trial_dir, cuda_device)
    
    with open(os.path.join(trial_dir, "aggregated_metrics.json"), "w") as f:
        json.dump({k: (None if str(v) == "nan" else v) for k, v in results.items()}, f, indent=2)
    
    mean_pcc = results["mean_val_pcc"]
    successful = results["successful_datasets"]
    total = results["total_datasets"]
    
    pcc_str = f"{mean_pcc:.4f}" if not str(mean_pcc) == "nan" else "nan"
    print(f"[Worker | GPU {cuda_device}] Trial {trial.number} 完成: mean_val_pcc={pcc_str}, 成功={successful}/{total}", flush=True)
    
    if str(mean_pcc) == "nan":
        return float("-inf")
    
    return mean_pcc


def run_worker(
    worker_id: int,
    gpu_id: str,
    study_name: str,
    storage_path: str,
    datasets: List[str],
    n_trials_per_worker: int,
    daemon: bool = False,
    max_trials: int = None,
):
    """单个worker进程"""
    print(f"[Worker {worker_id}] 启动于 GPU {gpu_id}", flush=True)
    
    storage = JournalStorage(JournalFileBackend(storage_path))
    study = optuna.load_study(study_name=study_name, storage=storage)
    
    callbacks = []
    if max_trials:
        callbacks.append(MaxTrialsCallback(max_trials))
    
    try:
        study.optimize(
            lambda trial: objective(trial, datasets, gpu_id),
            n_trials=n_trials_per_worker,
            callbacks=callbacks if callbacks else None,
        )
    except Exception as e:
        print(f"[Worker {worker_id}] 错误: {e}", flush=True)
    
    print(f"[Worker {worker_id}] 完成", flush=True)


def run_parallel_optimization(
    n_trials: int = 50,
    n_datasets: int = 5,
    cuda_devices: str = "0:2,1:2,2:2,3:2",
    daemon: bool = False,
):
    """运行并行优化"""
    print(f"启动并行超参数优化")
    print(f"  总trials: {n_trials}")
    print(f"  每trial数据集: {n_datasets}")
    print(f"  GPU配置: {cuda_devices}")
    print(f"  Daemon模式: {daemon}")
    
    os.makedirs(OUTPUT_ROOT, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_ROOT, "trial_outputs"), exist_ok=True)
    
    all_datasets = get_selected_datasets(30)
    import random
    random.seed(42)
    eval_datasets = random.sample(all_datasets, min(n_datasets, len(all_datasets)))
    print(f"  评估数据集: {eval_datasets}")
    
    worker_gpus = []
    for part in cuda_devices.split(","):
        gpu_id, count = part.split(":")
        worker_gpus.extend([gpu_id] * int(count))
    
    n_workers = len(worker_gpus)
    
    journal_path = os.path.join(OUTPUT_ROOT, "optuna_journal.log")
    storage = JournalStorage(JournalFileBackend(journal_path))
    
    study = optuna.create_study(
        study_name="marcp_hyperparam_optimization",
        storage=storage,
        direction="maximize",
        sampler=TPESampler(seed=42),
        pruner=MedianPruner(),
        load_if_exists=True,
    )
    
    # 计算还需要多少个新trial（基于已完成的trial数）
    current_complete = len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])
    needed_trials = max(0, n_trials - current_complete)
    
    if daemon:
        n_trials_per_worker = None
    elif needed_trials == 0:
        print(f"  已达到目标 {n_trials} 个COMPLETE trial，无需运行")
        return
    else:
        # 每个worker分配的trial数 = ceil(需要的trial数 / worker数)
        n_trials_per_worker = max(1, (needed_trials + n_workers - 1) // n_workers)
    
    print(f"  总workers: {n_workers}")
    print(f"  当前已完成: {current_complete}")
    print(f"  目标完成: {n_trials}")
    print(f"  还需要: {needed_trials}")
    print(f"  每worker trials: {n_trials_per_worker if n_trials_per_worker else '无限'}")
    
    processes = []
    for i, gpu_id in enumerate(worker_gpus):
        p = multiprocessing.Process(
            target=run_worker,
            args=(
                i, gpu_id, "marcp_hyperparam_optimization",
                journal_path, eval_datasets, n_trials_per_worker,
                daemon, n_trials if not daemon else None,
            ),
        )
        p.start()
        processes.append(p)
        print(f"  启动 Worker {i} 于 GPU {gpu_id}")
    
    try:
        for p in processes:
            p.join()
    except KeyboardInterrupt:
        print("\n收到中断信号，停止所有workers...")
        for p in processes:
            p.terminate()
    
    print("\n优化完成!")
    
    try:
        print(f"最佳trial: {study.best_trial.number}")
        print(f"最佳值: {study.best_value:.4f}")
        print(f"最佳参数: {study.best_params}")
    except Exception:
        print("没有完成的trial")


def main():
    parser = argparse.ArgumentParser(description="并行超参数优化")
    parser.add_argument("--n_trials", type=int, default=50, help="总trial数")
    parser.add_argument("--n_datasets", type=int, default=5, help="每trial评估的数据集数")
    parser.add_argument("--cuda", type=str, default="0:2,1:2,2:2,3:2", help="GPU配置 (格式: gpu_id:count,...)")
    parser.add_argument("--daemon", action="store_true", help="Daemon模式(持续运行)")
    
    args = parser.parse_args()
    
    run_parallel_optimization(
        n_trials=args.n_trials,
        n_datasets=args.n_datasets,
        cuda_devices=args.cuda,
        daemon=args.daemon,
    )


if __name__ == "__main__":
    main()
