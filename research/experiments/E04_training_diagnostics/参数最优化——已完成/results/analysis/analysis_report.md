# MARCP超参数优化分析报告

生成时间: 2026-03-19 07:29:47

## 1. 优化概况
- Study名称: marcp_hyperparam_optimization
- 总试验次数: 65
- 完成的试验: 58

## 2. 最佳超参数

**最佳验证PCC: 0.3371**

Trial编号: 64

| 参数 | 值 |
|------|-----|
| lr | 0.000785101 |
| bsz | 256 |
| weight_decay | 0.000114342 |
| alpha | 1.41942 |
| beta | 0.540218 |
| gamma | 1.79834 |
| grad_clip | 2 |
| max_stopping_step | 9 |

## 3. 推荐的训练命令

```bash
python -u marcp_train.py \
  --lr 0.000785101 \
  --bsz 256 \
  --weight_decay 0.000114342 \
  --alpha 1.41942 \
  --beta 0.540218 \
  --gamma 1.79834 \
  --grad_clip 2 \
  --max_stopping_step 9 \
  --input_channels 5 \
  --modelFramework Treat \
  --lossM CombinedProfileLoss \
  -train <train_file> \
  -validation <val_file> \
  -test <test_file> \
  --output <output_dir> \
  --task <task_name>
```

## 4. 参数重要性

| 参数 | 重要性 |
|------|--------|
| bsz | 0.716 |
| lr | 0.240 |
| gamma | 0.019 |
| beta | 0.010 |
| grad_clip | 0.007 |
| weight_decay | 0.003 |
| alpha | 0.003 |
| max_stopping_step | 0.002 |

## 5. 性能统计
- 最佳: 0.3371
- 最差: 0.2642
- 平均: 0.3203
- 中位数: 0.3258
- 标准差: 0.0169

## 6. 可视化图表

图表文件保存在 `visualizations/` 目录下:
- `optimization_history.html` - 优化历史
- `param_importances.html` - 参数重要性
- `parallel_coordinate.html` - 平行坐标图
- `slice_plot.html` - 参数切片图
- `optimization_progress.png` - 优化进度图
- `parameter_effects.png` - 参数效果图