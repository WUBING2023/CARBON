"""超参数优化配置"""
import os

MARCP_DIR = "/data/CARBON_PROJECT/MARCP"
DATA_ROOT = "/data/CARBON_PROJECT/数据/h5data"
OUTPUT_ROOT = "/data/CARBON_PROJECT/参数最优化/results"
TRAIN_SCRIPT = os.path.join(MARCP_DIR, "marcp_train.py")

SEARCH_SPACE = {
    "lr": {"type": "float", "low": 1e-4, "high": 1e-3, "log": True},
    "bsz": {"type": "categorical", "choices": [256, 512]},
    "weight_decay": {"type": "float", "low": 1e-4, "high": 1e-2, "log": True},
    "alpha": {"type": "float", "low": 0.5, "high": 2.0},
    "beta": {"type": "float", "low": 0.1, "high": 1.0},
    "gamma": {"type": "float", "low": 0.5, "high": 2.0},
    "grad_clip": {"type": "categorical", "choices": [0.5, 1.0, 2.0]},
    "max_stopping_step": {"type": "int", "low": 5, "high": 15},
}

FIXED_PARAMS = {
    "epochs": 100,
}

OPTIMIZATION_CONFIG = {
    "study_name": "marcp_hyperparam_optimization",
    "storage_name": "optuna_journal.log",
}


def get_data_paths(task_name: str):
    base = os.path.join(DATA_ROOT, task_name, "_5ch")
    train_path = os.path.join(base, f"{task_name}_chr1_5ch.h5")
    val_path = os.path.join(base, f"{task_name}_chr2_5ch.h5")
    test_path = os.path.join(base, f"{task_name}_chr3_5ch.h5")
    return train_path, val_path, test_path


def build_train_command(task_name: str, params: dict, output_dir: str, cuda_device: str = "0") -> str:
    train_path, val_path, test_path = get_data_paths(task_name)
    
    cmd_parts = [
        f"python3 -u {TRAIN_SCRIPT}",
        f"--CUDA {cuda_device}",
        f"--task {task_name}",
        f"-train {train_path}",
        f"-validation {val_path}",
        f"-test {test_path}",
        f"--output {output_dir}",
        f"--epochs {FIXED_PARAMS['epochs']}",
        "--input_channels 5",
        "--encode OneHot",
        "--modelFramework Treat",
        "--lossM CombinedProfileLoss",
        "--signalMax 0",
        "--seqLength 101",
        "--loss_scale 1.0",
        f"--lr {params['lr']:.6g}",
        f"--bsz {params['bsz']}",
        f"--weight_decay {params['weight_decay']:.6g}",
        f"--alpha {params['alpha']:.6g}",
        f"--beta {params['beta']:.6g}",
        f"--gamma {params['gamma']:.6g}",
        f"--grad_clip {params['grad_clip']}",
        f"--max_stopping_step {params['max_stopping_step']}",
    ]
    
    return " \\\n  ".join(cmd_parts)
