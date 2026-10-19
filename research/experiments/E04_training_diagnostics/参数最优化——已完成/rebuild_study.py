"""
重建Optuna Study

将已完成的trial（0-34）导入新的Optuna study
"""
import os
import json
import optuna
from optuna.storages import JournalStorage
try:
    from optuna.storages.journal import JournalFileBackend
except ImportError:
    from optuna.storages import JournalFileStorage as JournalFileBackend

from config import OUTPUT_ROOT, OPTIMIZATION_CONFIG

def rebuild_study():
    journal_path = os.path.join(OUTPUT_ROOT, OPTIMIZATION_CONFIG["storage_name"])
    
    # 删除旧的journal（如果存在）
    if os.path.exists(journal_path):
        os.remove(journal_path)
        print(f"已删除旧的journal: {journal_path}")
    
    # 创建新的storage和study
    storage = JournalStorage(JournalFileBackend(journal_path))
    study = optuna.create_study(
        study_name=OPTIMIZATION_CONFIG["study_name"],
        storage=storage,
        direction="maximize",
    )
    print(f"已创建新的study: {OPTIMIZATION_CONFIG['study_name']}")
    
    # 遍历trial_outputs目录，找到所有已完成的trial
    trial_outputs_dir = os.path.join(OUTPUT_ROOT, "trial_outputs")
    imported_count = 0
    
    for trial_dir in sorted(os.listdir(trial_outputs_dir)):
        if not trial_dir.startswith("trial_"):
            continue
        
        trial_num = int(trial_dir.split("_")[1])
        trial_path = os.path.join(trial_outputs_dir, trial_dir)
        
        # 读取超参数
        hp_file = os.path.join(trial_path, "hyperparameters.json")
        if not os.path.exists(hp_file):
            print(f"跳过 {trial_dir}: 缺少 hyperparameters.json")
            continue
        
        with open(hp_file, "r") as f:
            params = json.load(f)
        
        # 读取聚合指标
        metrics_file = os.path.join(trial_path, "aggregated_metrics.json")
        if not os.path.exists(metrics_file):
            print(f"跳过 {trial_dir}: 缺少 aggregated_metrics.json")
            continue
        
        with open(metrics_file, "r") as f:
            metrics = json.load(f)
        
        value = metrics.get("mean_val_pcc")
        if value is None:
            print(f"跳过 {trial_dir}: mean_val_pcc 为 None")
            continue
        
        # 添加trial到study
        # 【关键】必须与config.py中的SEARCH_SPACE完全一致！
        trial = optuna.trial.create_trial(
            params=params,
            distributions={
                "lr": optuna.distributions.FloatDistribution(1e-4, 1e-3, log=True),
                "bsz": optuna.distributions.CategoricalDistribution([256, 512]),
                "weight_decay": optuna.distributions.FloatDistribution(1e-4, 1e-2, log=True),
                "alpha": optuna.distributions.FloatDistribution(0.5, 2.0),
                "beta": optuna.distributions.FloatDistribution(0.1, 1.0),
                "gamma": optuna.distributions.FloatDistribution(0.5, 2.0),
                "grad_clip": optuna.distributions.CategoricalDistribution([0.5, 1.0, 2.0]),
                "max_stopping_step": optuna.distributions.IntDistribution(5, 15),
            },
            values=[value],
            state=optuna.trial.TrialState.COMPLETE,
        )
        study.add_trial(trial)
        imported_count += 1
        print(f"导入 {trial_dir}: value={value:.4f}")
    
    print(f"\n导入完成: {imported_count} 个trial")
    print(f"最佳trial: {study.best_trial.number}, value={study.best_value:.4f}")
    print(f"最佳参数: {study.best_params}")

if __name__ == "__main__":
    rebuild_study()
