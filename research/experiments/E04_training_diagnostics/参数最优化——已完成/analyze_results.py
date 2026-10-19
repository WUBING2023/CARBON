"""
超参数优化结果分析和可视化脚本
"""

import os
import json
import argparse
from typing import Dict, List, Optional, Tuple
from datetime import datetime

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("警告: matplotlib未安装，图表功能不可用")

try:
    import optuna
    from optuna.visualization import (
        plot_optimization_history,
        plot_param_importances,
        plot_parallel_coordinate,
        plot_slice,
    )
    HAS_OPTUNA_VIZ = True
except ImportError:
    HAS_OPTUNA_VIZ = False

from config import OUTPUT_ROOT, OPTIMIZATION_CONFIG, SEARCH_SPACE


def load_study() -> Optional[optuna.Study]:
    """加载Optuna study"""
    from optuna.storages import JournalStorage
    try:
        from optuna.storages.journal import JournalFileBackend
    except ImportError:
        from optuna.storages import JournalFileStorage as JournalFileBackend
    
    storage_path = os.path.join(OUTPUT_ROOT, OPTIMIZATION_CONFIG["storage_name"])
    if not os.path.isfile(storage_path):
        print(f"Study文件不存在: {storage_path}")
        return None
    
    storage = JournalStorage(JournalFileBackend(storage_path))
    study = optuna.load_study(
        study_name=OPTIMIZATION_CONFIG["study_name"],
        storage=storage,
    )
    return study


def get_trials_dataframe(study: optuna.Study) -> pd.DataFrame:
    """将trials转换为DataFrame"""
    data = []
    for trial in study.trials:
        if trial.state == optuna.trial.TrialState.COMPLETE:
            row = {
                "trial_number": trial.number,
                "value": trial.value,
                "datetime_start": trial.datetime_start,
                "datetime_complete": trial.datetime_complete,
            }
            row.update(trial.params)
            data.append(row)
    
    df = pd.DataFrame(data)
    if "datetime_start" in df.columns and "datetime_complete" in df.columns:
        df["duration_minutes"] = (
            df["datetime_complete"] - df["datetime_start"]
        ).dt.total_seconds() / 60
    
    return df


def print_study_summary(study: optuna.Study):
    """打印study摘要"""
    print("\n" + "=" * 60)
    print("超参数优化结果摘要")
    print("=" * 60)
    
    completed_trials = [
        t for t in study.trials 
        if t.state == optuna.trial.TrialState.COMPLETE
    ]
    
    print(f"\n总试验次数: {len(study.trials)}")
    print(f"完成的试验: {len(completed_trials)}")
    print(f"优化方向: {study.direction}")
    
    if completed_trials:
        values = [t.value for t in completed_trials if t.value is not None]
        if values:
            print(f"\n性能统计 (val_pcc):")
            print(f"  最佳: {max(values):.4f}")
            print(f"  最差: {min(values):.4f}")
            print(f"  平均: {np.mean(values):.4f}")
            print(f"  标准差: {np.std(values):.4f}")
    
    print("\n最佳超参数:")
    print("-" * 40)
    for name, value in study.best_params.items():
        if isinstance(value, float):
            print(f"  {name}: {value:.6g}")
        else:
            print(f"  {name}: {value}")
    print(f"\n最佳验证PCC: {study.best_value:.4f}")
    print(f"最佳trial编号: {study.best_trial.number}")


def analyze_parameter_importance(study: optuna.Study) -> Dict[str, float]:
    """分析参数重要性"""
    try:
        importance = optuna.importance.get_param_importances(study)
        
        print("\n参数重要性排序:")
        print("-" * 40)
        for name, imp in sorted(importance.items(), key=lambda x: -x[1]):
            bar = "█" * int(imp * 30)
            print(f"  {name:20s} {imp:.3f} {bar}")
        
        return importance
    except Exception as e:
        print(f"参数重要性分析失败: {e}")
        return {}


def analyze_top_trials(
    study: optuna.Study, 
    n_top: int = 10,
) -> pd.DataFrame:
    """分析表现最好的trials"""
    df = get_trials_dataframe(study)
    
    if df.empty:
        print("没有完成的trials")
        return df
    
    df_sorted = df.sort_values("value", ascending=False).head(n_top)
    
    print(f"\n表现最好的 {n_top} 次试验:")
    print("-" * 60)
    
    param_cols = list(SEARCH_SPACE.keys())
    display_cols = ["trial_number", "value"] + param_cols
    display_cols = [c for c in display_cols if c in df_sorted.columns]
    
    print(df_sorted[display_cols].to_string(index=False))
    
    return df_sorted


def analyze_parameter_correlations(study: optuna.Study) -> pd.DataFrame:
    """分析参数与性能的相关性"""
    df = get_trials_dataframe(study)
    
    if df.empty or "value" not in df.columns:
        return pd.DataFrame()
    
    param_cols = [c for c in SEARCH_SPACE.keys() if c in df.columns]
    
    correlations = {}
    for col in param_cols:
        if df[col].dtype in [np.float64, np.int64]:
            corr = df[col].corr(df["value"])
            correlations[col] = corr
    
    print("\n参数与性能的相关性:")
    print("-" * 40)
    for name, corr in sorted(correlations.items(), key=lambda x: -abs(x[1])):
        direction = "↑" if corr > 0 else "↓"
        print(f"  {name:20s} {corr:+.3f} {direction}")
    
    return pd.DataFrame(correlations.items(), columns=["parameter", "correlation"])


def generate_optuna_visualizations(
    study: optuna.Study, 
    output_dir: str,
):
    """生成Optuna内置的可视化图表"""
    if not HAS_OPTUNA_VIZ:
        print("Optuna可视化功能不可用")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n生成可视化图表...")
    
    try:
        fig = plot_optimization_history(study)
        fig.write_html(os.path.join(output_dir, "optimization_history.html"))
        print("  - optimization_history.html")
    except Exception as e:
        print(f"  优化历史图失败: {e}")
    
    try:
        fig = plot_param_importances(study)
        fig.write_html(os.path.join(output_dir, "param_importances.html"))
        print("  - param_importances.html")
    except Exception as e:
        print(f"  参数重要性图失败: {e}")
    
    try:
        fig = plot_parallel_coordinate(study)
        fig.write_html(os.path.join(output_dir, "parallel_coordinate.html"))
        print("  - parallel_coordinate.html")
    except Exception as e:
        print(f"  平行坐标图失败: {e}")
    
    try:
        fig = plot_slice(study)
        fig.write_html(os.path.join(output_dir, "slice_plot.html"))
        print("  - slice_plot.html")
    except Exception as e:
        print(f"  切片图失败: {e}")


def generate_matplotlib_plots(
    study: optuna.Study,
    output_dir: str,
):
    """生成Matplotlib图表"""
    if not HAS_MATPLOTLIB:
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    df = get_trials_dataframe(study)
    if df.empty:
        return
    
    plt.style.use('seaborn-v0_8-whitegrid')
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(df["trial_number"], df["value"], 'b-', alpha=0.5, label="trial value")
    best_so_far = df["value"].cummax()
    ax.plot(df["trial_number"], best_so_far, 'r-', linewidth=2, label="best so far")
    ax.set_xlabel("Trial Number")
    ax.set_ylabel("Validation PCC")
    ax.set_title("Optimization Progress")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "optimization_progress.png"), dpi=150)
    plt.close()
    print("  - optimization_progress.png")
    
    param_cols = [c for c in SEARCH_SPACE.keys() if c in df.columns]
    n_params = len(param_cols)
    n_cols = min(3, n_params)
    n_rows = (n_params + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    axes = np.atleast_2d(axes).flatten()
    
    for idx, param in enumerate(param_cols):
        ax = axes[idx]
        scatter = ax.scatter(
            df[param], df["value"], 
            c=df["trial_number"], 
            cmap="viridis", 
            alpha=0.7
        )
        ax.set_xlabel(param)
        ax.set_ylabel("Validation PCC")
        ax.set_title(f"Effect of {param}")
    
    for idx in range(len(param_cols), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "parameter_effects.png"), dpi=150)
    plt.close()
    print("  - parameter_effects.png")


def generate_report(
    study: optuna.Study,
    output_dir: str,
) -> str:
    """生成完整的分析报告"""
    os.makedirs(output_dir, exist_ok=True)
    
    report_lines = []
    report_lines.append("# MARCP超参数优化分析报告")
    report_lines.append(f"\n生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    report_lines.append("\n## 1. 优化概况")
    report_lines.append(f"- Study名称: {study.study_name}")
    report_lines.append(f"- 总试验次数: {len(study.trials)}")
    completed = len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])
    report_lines.append(f"- 完成的试验: {completed}")
    
    report_lines.append("\n## 2. 最佳超参数")
    report_lines.append(f"\n**最佳验证PCC: {study.best_value:.4f}**")
    report_lines.append(f"\nTrial编号: {study.best_trial.number}")
    report_lines.append("\n| 参数 | 值 |")
    report_lines.append("|------|-----|")
    for name, value in study.best_params.items():
        if isinstance(value, float):
            report_lines.append(f"| {name} | {value:.6g} |")
        else:
            report_lines.append(f"| {name} | {value} |")
    
    report_lines.append("\n## 3. 推荐的训练命令")
    report_lines.append("\n```bash")
    cmd_parts = ["python -u marcp_train.py \\"]
    for name, value in study.best_params.items():
        if isinstance(value, float):
            cmd_parts.append(f"  --{name} {value:.6g} \\")
        else:
            cmd_parts.append(f"  --{name} {value} \\")
    cmd_parts.append("  --input_channels 5 \\")
    cmd_parts.append("  --modelFramework Treat \\")
    cmd_parts.append("  --lossM CombinedProfileLoss \\")
    cmd_parts.append("  -train <train_file> \\")
    cmd_parts.append("  -validation <val_file> \\")
    cmd_parts.append("  -test <test_file> \\")
    cmd_parts.append("  --output <output_dir> \\")
    cmd_parts.append("  --task <task_name>")
    report_lines.append("\n".join(cmd_parts))
    report_lines.append("```")
    
    try:
        importance = optuna.importance.get_param_importances(study)
        report_lines.append("\n## 4. 参数重要性")
        report_lines.append("\n| 参数 | 重要性 |")
        report_lines.append("|------|--------|")
        for name, imp in sorted(importance.items(), key=lambda x: -x[1]):
            report_lines.append(f"| {name} | {imp:.3f} |")
    except Exception:
        pass
    
    df = get_trials_dataframe(study)
    if not df.empty and "value" in df.columns:
        report_lines.append("\n## 5. 性能统计")
        values = df["value"].dropna()
        report_lines.append(f"- 最佳: {values.max():.4f}")
        report_lines.append(f"- 最差: {values.min():.4f}")
        report_lines.append(f"- 平均: {values.mean():.4f}")
        report_lines.append(f"- 中位数: {values.median():.4f}")
        report_lines.append(f"- 标准差: {values.std():.4f}")
    
    report_lines.append("\n## 6. 可视化图表")
    report_lines.append("\n图表文件保存在 `visualizations/` 目录下:")
    report_lines.append("- `optimization_history.html` - 优化历史")
    report_lines.append("- `param_importances.html` - 参数重要性")
    report_lines.append("- `parallel_coordinate.html` - 平行坐标图")
    report_lines.append("- `slice_plot.html` - 参数切片图")
    report_lines.append("- `optimization_progress.png` - 优化进度图")
    report_lines.append("- `parameter_effects.png` - 参数效果图")
    
    report_content = "\n".join(report_lines)
    report_file = os.path.join(output_dir, "analysis_report.md")
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(report_content)
    
    print(f"\n报告已保存到: {report_file}")
    return report_file


def export_best_params(study: optuna.Study, output_file: str = None) -> str:
    """导出最佳参数到JSON文件"""
    if output_file is None:
        output_file = os.path.join(OUTPUT_ROOT, "best_hyperparameters.json")
    
    data = {
        "best_params": study.best_params,
        "best_value": study.best_value,
        "best_trial_number": study.best_trial.number,
        "n_trials": len(study.trials),
        "n_completed": len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]),
        "timestamp": datetime.now().isoformat(),
    }
    
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"最佳参数已导出到: {output_file}")
    return output_file


def main():
    parser = argparse.ArgumentParser(description="分析超参数优化结果")
    parser.add_argument("--output", type=str, default=None,
                       help="输出目录 (默认: results/analysis)")
    parser.add_argument("--top", type=int, default=10,
                       help="显示表现最好的N次试验 (默认: 10)")
    parser.add_argument("--no-viz", action="store_true",
                       help="不生成可视化图表")
    parser.add_argument("--export", action="store_true",
                       help="导出最佳参数到JSON")
    
    args = parser.parse_args()
    
    study = load_study()
    if study is None:
        print("无法加载study，请先运行优化")
        return
    
    print_study_summary(study)
    
    analyze_parameter_importance(study)
    
    analyze_top_trials(study, n_top=args.top)
    
    analyze_parameter_correlations(study)
    
    output_dir = args.output or os.path.join(OUTPUT_ROOT, "analysis")
    
    if not args.no_viz:
        viz_dir = os.path.join(output_dir, "visualizations")
        generate_optuna_visualizations(study, viz_dir)
        generate_matplotlib_plots(study, viz_dir)
    
    generate_report(study, output_dir)
    
    if args.export:
        export_best_params(study)


if __name__ == "__main__":
    main()
