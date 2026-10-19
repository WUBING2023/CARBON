# 超参数优化代码结构详解

## 目录结构
```
参数最优化/
├── config.py                    # 配置文件（路径、搜索空间、命令构建）
├── dataset_selector.py          # 数据集选择器
├── run_optimization_parallel.py # 核心并行优化脚本
├── run_bayesian_rounds.py       # 分阶段优化控制脚本
├── rebuild_study.py             # Optuna Study重建脚本
├── analyze_results.py           # 结果分析脚本
├── README.md                    # 使用说明
├── CODE_STRUCTURE.md            # 本文档
└── results/                     # 输出目录
    ├── optuna_journal.log       # Optuna存储文件
    ├── selected_datasets.json   # 已选数据集列表
    └── trial_outputs/           # 每个trial的输出
        └── trial_XXXX/
            ├── hyperparameters.json
            ├── aggregated_metrics.json
            └── {dataset_name}/  # 各数据集训练结果
```

---

## 1. config.py - 配置文件

### 1.1 路径配置（第4-7行）
```python
MARCP_DIR = "/data/CARBON_PROJECT/MARCP"          # MARCP模型目录
DATA_ROOT = "/data/CARBON_PROJECT/数据/h5data"     # 数据集根目录
OUTPUT_ROOT = "/data/CARBON_PROJECT/参数最优化/results"  # 输出目录
TRAIN_SCRIPT = os.path.join(MARCP_DIR, "marcp_train.py")  # 训练脚本路径
```

### 1.2 搜索空间定义（第9-18行）
```python
SEARCH_SPACE = {
    "lr": {"type": "float", "low": 1e-4, "high": 1e-3, "log": True},
    # 学习率：对数均匀分布，范围[0.0001, 0.001]
    
    "bsz": {"type": "categorical", "choices": [256, 512]},
    # 批大小：从256或512中选择
    
    "weight_decay": {"type": "float", "low": 1e-4, "high": 1e-2, "log": True},
    # 权重衰减：对数均匀分布，范围[0.0001, 0.01]
    
    "alpha": {"type": "float", "low": 0.5, "high": 2.0},
    # Profile Shape Loss权重：均匀分布[0.5, 2.0]
    
    "beta": {"type": "float", "low": 0.1, "high": 1.0},
    # Jitter Loss权重：均匀分布[0.1, 1.0]
    
    "gamma": {"type": "float", "low": 0.5, "high": 2.0},
    # Count Loss权重：均匀分布[0.5, 2.0]
    
    "grad_clip": {"type": "categorical", "choices": [0.5, 1.0, 2.0]},
    # 梯度裁剪：从0.5, 1.0, 2.0中选择
    
    "max_stopping_step": {"type": "int", "low": 5, "high": 15},
    # 早停耐心值：整数范围[5, 15]
}
```

### 1.3 固定参数（第20-22行）
```python
FIXED_PARAMS = {
    "epochs": 100,  # 每次训练最大epoch数（可能被早停截断）
}
```

### 1.4 Optuna配置（第24-27行）
```python
OPTIMIZATION_CONFIG = {
    "study_name": "marcp_hyperparam_optimization",  # Study名称（唯一标识）
    "storage_name": "optuna_journal.log",           # 存储文件名
}
```

### 1.5 get_data_paths函数（第30-35行）
```python
def get_data_paths(task_name: str):
    """根据任务名获取训练/验证/测试数据路径"""
    base = os.path.join(DATA_ROOT, task_name, "_5ch")
    # 数据目录结构：h5data/{task_name}/_5ch/
    
    train_path = os.path.join(base, f"{task_name}_chr1_5ch.h5")  # 训练：chr1
    val_path = os.path.join(base, f"{task_name}_chr2_5ch.h5")    # 验证：chr2
    test_path = os.path.join(base, f"{task_name}_chr3_5ch.h5")   # 测试：chr3
    return train_path, val_path, test_path
```

### 1.6 build_train_command函数（第38-67行）
```python
def build_train_command(task_name: str, params: dict, output_dir: str, cuda_device: str = "0") -> str:
    """构建marcp_train.py的命令行"""
    train_path, val_path, test_path = get_data_paths(task_name)
    
    cmd_parts = [
        f"python3 -u {TRAIN_SCRIPT}",       # -u: unbuffered output
        f"--CUDA {cuda_device}",             # GPU设备号
        f"--task {task_name}",               # 任务名称
        f"-train {train_path}",              # 训练数据
        f"-validation {val_path}",           # 验证数据
        f"-test {test_path}",                # 测试数据
        f"--output {output_dir}",            # 输出目录
        f"--epochs {FIXED_PARAMS['epochs']}", # 最大epoch
        "--input_channels 5",                # 输入通道数（5ch数据）
        "--encode OneHot",                   # 编码方式
        "--modelFramework Treat",            # 模型框架
        "--lossM CombinedProfileLoss",       # 损失函数
        "--signalMax 0",                     # 信号最大值
        "--seqLength 101",                   # 序列长度
        "--loss_scale 1.0",                  # 损失缩放
        f"--lr {params['lr']:.6g}",          # 【优化参数】学习率
        f"--bsz {params['bsz']}",            # 【优化参数】批大小
        f"--weight_decay {params['weight_decay']:.6g}",  # 【优化参数】权重衰减
        f"--alpha {params['alpha']:.6g}",    # 【优化参数】alpha
        f"--beta {params['beta']:.6g}",      # 【优化参数】beta
        f"--gamma {params['gamma']:.6g}",    # 【优化参数】gamma
        f"--grad_clip {params['grad_clip']}", # 【优化参数】梯度裁剪
        f"--max_stopping_step {params['max_stopping_step']}",  # 【优化参数】早停
    ]
    
    return " \\\n  ".join(cmd_parts)  # 用反斜杠连接，便于阅读
```

---

## 2. dataset_selector.py - 数据集选择器

### 2.1 核心功能
- 从223个数据集中选择30个
- 按细胞系平衡分布（HepG2和K562各约15个）
- 结果保存到 `selected_datasets.json`

### 2.2 get_selected_datasets函数（第62-73行）
```python
def get_selected_datasets(n_select: int = 30):
    """获取已选择的数据集列表"""
    if os.path.exists(SELECTED_FILE):
        # 如果已存在选择文件，直接读取
        with open(SELECTED_FILE, "r") as f:
            data = json.load(f)
            return data.get("datasets", [])
    
    # 否则进行新的选择
    selected = select_balanced_datasets(n_select)
    with open(SELECTED_FILE, "w") as f:
        json.dump({"datasets": selected, "count": len(selected)}, f, indent=2)
    
    return selected
```

---

## 3. run_optimization_parallel.py - 核心并行优化脚本

### 3.1 MaxTrialsCallback类（第26-34行）
```python
class MaxTrialsCallback:
    """回调函数：当COMPLETE的trial数达到max_trials时停止优化"""
    def __init__(self, max_trials: int):
        self.max_trials = max_trials
    
    def __call__(self, study: optuna.Study, trial: optuna.trial.FrozenTrial):
        # 每个trial完成后调用
        n_complete = len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])
        if n_complete >= self.max_trials:
            study.stop()  # 停止优化
```
**作用**：防止超过目标trial数

### 3.2 sample_hyperparams函数（第37-49行）
```python
def sample_hyperparams(trial: optuna.Trial) -> Dict[str, Any]:
    """根据SEARCH_SPACE定义，让Optuna采样超参数"""
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
```
**作用**：TPE采样器根据历史trial智能选择参数

### 3.3 parse_validation_metrics函数（第52-73行）
```python
def parse_validation_metrics(output_dir: str) -> Dict[str, float]:
    """从val_metric.txt解析验证集指标"""
    # 读取val_metric.txt的最后一行
    # 提取: pearson (PCC), PR_AUC, spearman, loss
```
**作用**：获取训练结果用于优化

### 3.4 run_single_training函数（第76-114行）
```python
def run_single_training(task_name, params, output_dir, cuda_device, timeout=7200):
    """运行单个数据集的训练"""
    # 1. 创建输出目录
    # 2. 构建并保存训练命令
    # 3. 使用subprocess执行训练
    # 4. 解析结果
    # 5. 判断成功：检查val_metric.txt是否存在
```
**注意**：成功判断基于`val_metric.txt`存在，而非进程返回码

### 3.5 run_multi_dataset_evaluation函数（第117-158行）
```python
def run_multi_dataset_evaluation(params, datasets, trial_dir, cuda_device):
    """在多个数据集上评估同一组超参数"""
    # 对每个数据集调用run_single_training
    # 聚合结果：计算mean和std
    # 返回: mean_val_pcc, mean_val_prauc, etc.
```
**作用**：一个trial = 在5个数据集上的平均性能

### 3.6 objective函数（第161-189行）
```python
def objective(trial: optuna.Trial, datasets: List[str], cuda_device: str) -> float:
    """Optuna优化目标函数"""
    # 1. 采样超参数
    params = sample_hyperparams(trial)
    
    # 2. 创建trial输出目录
    trial_dir = f"trial_{trial.number:04d}"
    
    # 3. 保存超参数
    json.dump(params, "hyperparameters.json")
    
    # 4. 在多个数据集上评估
    results = run_multi_dataset_evaluation(params, datasets, trial_dir, cuda_device)
    
    # 5. 保存聚合结果
    json.dump(results, "aggregated_metrics.json")
    
    # 6. 返回优化目标（mean_val_pcc）
    return results["mean_val_pcc"]  # 最大化PCC
```
**关键**：返回`mean_val_pcc`作为优化目标，Optuna会最大化此值

### 3.7 run_worker函数（第192-221行）
```python
def run_worker(worker_id, gpu_id, study_name, storage_path, datasets, 
               n_trials_per_worker, daemon, max_trials):
    """单个Worker进程的执行逻辑"""
    # 1. 连接到共享的Optuna Study
    storage = JournalStorage(JournalFileBackend(storage_path))
    study = optuna.load_study(study_name=study_name, storage=storage)
    
    # 2. 设置回调（如果有max_trials限制）
    callbacks = [MaxTrialsCallback(max_trials)] if max_trials else []
    
    # 3. 执行优化
    study.optimize(
        lambda trial: objective(trial, datasets, gpu_id),
        n_trials=n_trials_per_worker,  # 此worker运行的trial数
        callbacks=callbacks,
    )
```
**重要**：多个Worker共享同一个Study，通过JournalStorage同步

### 3.8 run_parallel_optimization函数（第224-313行）
```python
def run_parallel_optimization(n_trials, n_datasets, cuda_devices, daemon):
    """主函数：启动并行优化"""
    
    # 1. 解析GPU配置
    # cuda_devices = "0:4,1:2,2:2,3:2" 
    # → worker_gpus = ['0','0','0','0','1','1','2','2','3','3']
    # → 10个workers
    
    # 2. 创建/加载Optuna Study
    study = optuna.create_study(
        study_name="marcp_hyperparam_optimization",
        storage=JournalStorage(...),
        direction="maximize",      # 最大化目标
        sampler=TPESampler(seed=42),  # TPE采样器
        load_if_exists=True,       # 支持断点续传
    )
    
    # 3. 【关键修复】计算每个worker的trial数
    current_complete = len([完成的trial])
    needed_trials = n_trials - current_complete  # 还需要多少
    n_trials_per_worker = ceil(needed_trials / n_workers)  # 每worker分配数
    
    # 4. 启动多进程
    for i, gpu_id in enumerate(worker_gpus):
        p = multiprocessing.Process(target=run_worker, args=(...))
        p.start()
    
    # 5. 等待所有进程完成
    for p in processes:
        p.join()
```

---

## 4. run_bayesian_rounds.py - 分阶段优化控制脚本

### 4.1 设计目的
控制多轮贝叶斯优化，每轮之间让TPE采样器充分利用新数据

### 4.2 get_completed_trials函数（第20-43行）
```python
def get_completed_trials():
    """获取当前已完成的trial数量"""
    # 读取Optuna journal
    # 统计COMPLETE状态的trial数
```

### 4.3 run_optimization_round函数（第46-90行）
```python
def run_optimization_round(target_trials, cuda_config, round_num):
    """运行一轮优化"""
    # 1. 计算本轮需要新增的trial数
    current = get_completed_trials()
    new_trials = target_trials - current
    
    # 2. 如果已达目标，跳过
    if new_trials <= 0:
        return
    
    # 3. 调用run_optimization_parallel.py
    subprocess.run(["python3", "run_optimization_parallel.py", 
                    "--n_trials", target_trials, ...])
```

### 4.4 main函数流程（第123-192行）
```python
def main():
    # 参数：
    # --initial_trials: 第1阶段完成数（当前是35）
    # --rounds: 额外轮数（当前还需要2轮）
    # --trials_per_round: 每轮trial数（10）
    
    for round_idx in range(args.rounds):
        target = args.initial_trials + (round_idx + 1) * args.trials_per_round
        run_optimization_round(target, cuda_config, round_num)
```

---

## 5. 数据流图

```
用户启动 run_bayesian_rounds.py
    │
    ▼
┌─────────────────────────────────┐
│ 第N轮优化                        │
│ target_trials = 35 + N*10       │
└───────────────┬─────────────────┘
                │
                ▼
    run_optimization_parallel.py
                │
    ┌───────────┴───────────┐
    │   解析GPU配置          │
    │   0:4,1:2,2:2,3:2     │
    │   → 10个workers       │
    └───────────┬───────────┘
                │
    ┌───────────┴───────────┐
    │   计算每worker的trial数│
    │   needed = target - current
    │   per_worker = ceil(needed/10)
    └───────────┬───────────┘
                │
    ┌───────────┴───────────┐
    │   启动10个Worker进程   │
    │   共享JournalStorage   │
    └───────────┬───────────┘
                │
        ┌───────┴───────┐
        ▼               ▼
    Worker 0-3      Worker 4-5     ...
    (GPU 0)         (GPU 1)
        │               │
        ▼               ▼
    study.optimize()
        │
        ▼
    objective(trial)
        │
        ├─→ sample_hyperparams() ← TPE采样
        │
        ├─→ run_multi_dataset_evaluation()
        │       │
        │       ├─→ run_single_training(dataset1)
        │       ├─→ run_single_training(dataset2)
        │       └─→ ... (共5个数据集)
        │
        └─→ return mean_val_pcc
                │
                ▼
        Trial完成，结果写入Journal
        TPE更新概率模型
```

---

## 6. 当前状态

### 6.1 已完成的Trial
- **Trial 0-24**: 第1阶段随机探索（25个）
- **Trial 25-34**: 第2阶段贝叶斯优化（10个，原编号50-59）
- **总计**: 35个COMPLETE

### 6.2 待运行的Trial
- **Trial 35-44**: 第3阶段（10个）
- **Trial 45-54**: 第4阶段（10个）
- **总计**: 还需20个，共55个

### 6.3 最佳结果
- **最佳Trial**: 31
- **最佳PCC**: 0.3349
- **最佳参数**: lr=0.000685, bsz=256, weight_decay=0.000362, alpha=1.59, beta=0.345, gamma=1.08, grad_clip=1.0, max_stopping_step=8

---

## 7. 启动命令

### 运行第3-4阶段
```bash
cd /data/CARBON_PROJECT/参数最优化

# 后台运行
nohup python3 run_bayesian_rounds.py \
    --initial_trials 35 \
    --rounds 2 \
    --trials_per_round 10 \
    --cuda "0:4,1:2,2:2,3:2" \
    > optimization_phase3-4.log 2>&1 &

# 监控
tail -f optimization_phase3-4.log
```

### 预期流程
```
第3轮: 目标45个 → 运行10个新trial (35-44)
第4轮: 目标55个 → 运行10个新trial (45-54)
```
