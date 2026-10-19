"""数据集选择器 - 从223个数据集中选择30个平衡分布的数据集"""
import os
import json
import random
from collections import defaultdict
from config import DATA_ROOT

SELECTED_FILE = "/data/CARBON_PROJECT/参数最优化/selected_datasets.json"


def get_all_datasets():
    """获取所有可用的5ch数据集"""
    datasets = set()
    for f in os.listdir(DATA_ROOT):
        if f.endswith("_5ch_train.h5"):
            name = f.replace("_5ch_train.h5", "")
            train_path = os.path.join(DATA_ROOT, f"{name}_5ch_train.h5")
            val_path = os.path.join(DATA_ROOT, f"{name}_5ch_val.h5")
            if os.path.exists(train_path) and os.path.exists(val_path):
                datasets.add(name)
    return sorted(datasets)


def parse_cell_line(name: str) -> str:
    """从数据集名称解析细胞系"""
    parts = name.split("_")
    if len(parts) >= 2:
        return parts[1]
    return "Unknown"


def select_balanced_datasets(n_select: int = 30, seed: int = 42):
    """选择平衡分布的数据集"""
    all_datasets = get_all_datasets()
    
    by_cell_line = defaultdict(list)
    for ds in all_datasets:
        cell = parse_cell_line(ds)
        by_cell_line[cell].append(ds)
    
    random.seed(seed)
    selected = []
    
    cell_lines = sorted(by_cell_line.keys())
    per_cell = n_select // len(cell_lines)
    remainder = n_select % len(cell_lines)
    
    for i, cell in enumerate(cell_lines):
        n = per_cell + (1 if i < remainder else 0)
        available = by_cell_line[cell]
        random.shuffle(available)
        selected.extend(available[:n])
    
    if len(selected) < n_select:
        remaining = [ds for ds in all_datasets if ds not in selected]
        random.shuffle(remaining)
        selected.extend(remaining[:n_select - len(selected)])
    
    return selected[:n_select]


def get_selected_datasets(n_select: int = 30):
    """获取已选择的数据集，如不存在则创建"""
    if os.path.exists(SELECTED_FILE):
        with open(SELECTED_FILE, "r") as f:
            data = json.load(f)
            return data.get("datasets", [])
    
    selected = select_balanced_datasets(n_select)
    with open(SELECTED_FILE, "w") as f:
        json.dump({"datasets": selected, "count": len(selected)}, f, indent=2)
    
    return selected


if __name__ == "__main__":
    datasets = get_selected_datasets()
    print(f"选择了 {len(datasets)} 个数据集:")
    for ds in datasets:
        print(f"  - {ds}")
