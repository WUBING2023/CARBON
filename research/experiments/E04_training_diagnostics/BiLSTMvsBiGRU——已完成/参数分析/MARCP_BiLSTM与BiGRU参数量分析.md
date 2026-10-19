# MARCP 消融对比：BiLSTM vs BiGRU 参数量分析

本文档记录在使用默认 MARCP 超参、输入 `(B, L, C) = (1, 101, 5)` 时，**仅将时序混合模块由 `ModernBiGRU` 换为对称实现的 `ModernBiLSTM`** 后的参数量对比结论。数据由 [`compare_marcp_rnn_params.py`](compare_marcp_rnn_params.py) 结合 PyTorch 统计与 `torchinfo.summary` 得到。

## 1. 模型设置（与官方 `MARCP` 默认一致）

| 项 | 取值 |
|----|------|
| `in_channels` | 5 |
| `gru_hidden`（单向隐藏维，双向输出维为 `2×`） | 64 |
| `motif_out_per_branch` | 32 |
| `motif_kernels` | (3, 5, 7) |
| `profile_hidden` | (64, 32) |
| `dropOut` | 0.2 |
| 消融开关 | 均未开启（`ablation_no_se/multiscale/gru` = False） |

BiLSTM 版在脚本中为 `MARCPBiLSTM`：除 `mixer` 内 RNN 由 `nn.GRU` 改为 `nn.LSTM` 外，`proj` / `residual_proj` / `LayerNorm` / `Dropout` 及前后模块与 BiGRU 版逐行对齐。

## 2. 总参数量对比

| 版本 | 总参数 | 可训练参数 |
|------|--------|------------|
| BiGRU（官方 `MARCP`） | **100,826** | **100,826** |
| BiLSTM（`MARCPBiLSTM`） | **117,466** | **117,466** |

- **绝对差**：BiLSTM 比 BiGRU 多 **16,640** 个参数（$\Delta = \text{LSTM} - \text{GRU}$）。
- **相对比例**：$\text{LSTM} / \text{GRU} \approx \mathbf{1.165}$，即 LSTM 版总参数约为 GRU 版的 **16.5%** 增量。

由于两版仅在 `mixer` 内的循环子模块不同，**全模型参数量之差应等于 BiLSTM 与 BiGRU 子模块参数之差**。下面专门对比这一块。

## 3. 时序模块（mixer）内：GRU vs LSTM 参数量

`torchinfo` 在相同深度下对 `ModernBiGRU` / `ModernBiLSTM` 汇报的 RNN 子模块参数为：

| 子模块 | 参数量 | 说明 |
|--------|--------|------|
| **GRU**（双向、单层，`hidden_size=64`） | **49,920** | 3 组门（reset/update/new），每组对 `input` 与 `hidden` 各一张权重矩阵及偏置 |
| **LSTM**（双向、单层，同上） | **66,560** | 4 组门（input/forget/cell/output），权重与偏置规模按 PyTorch 实现累计 |

- **RNN 部分参数量差**：$66{,}560 - 49{,}920 = \mathbf{16{,}640}$，与**全模型总差 16,640 完全一致**。
- **结论**：在你当前的 MARCP 默认配置下，**BiLSTM 相对 BiGRU 的额外参数全部来自单层双向 LSTM 相对单层双向 GRU 的容量增加**；卷积多尺度、`profile_head`、`count_head` 等其余部分参数量两版相同。

### 简要机理（直观理解）

在输入维度与隐藏维相同的前提下，标准实现中 **LSTM 的可训练权重规模通常大于 GRU**，主要因为：

- GRU 对每个方向、每个时间步用 **3 组** 线性变换（对上一隐藏与当前输入）组合门控；
- LSTM 使用 **4 组**（输入门、遗忘门、细胞候选、输出门），参数项更多。

因此，**不调整 `hidden_size` 或层数时，将 `ModernBiGRU` 换成对称的 `ModernBiLSTM`，参数量上升是预期行为**，量级由上表中的 **约 1.66 万** 刻画。

## 4. 计算量（torchinfo 估算，仅供参考）

| 版本 | Total mult-adds (M) |
|------|---------------------|
| BiGRU | ≈ 9.23 |
| BiLSTM | ≈ 10.91 |

LSTM 前向通常也比同规模的 GRU 更重，与参数量与门控结构一致；具体数值随序列长度与实现细节略有变化。

## 5. 如何复现

```bash
pip install -r /path/to/RBP/BiLSTMvsBiGRU/参数分析/requirements.txt
python3 /path/to/RBP/BiLSTMvsBiGRU/参数分析/compare_marcp_rnn_params.py
```

若开启 `--ablation_no_gru`，两侧均去掉 RNN 路径，总参数应一致（仅 residual + LayerNorm 等），可用于 Sanity check。

---

**摘要**：在默认 MARCP 结构下，**BiLSTM 总参数 **117,466** > BiGRU **100,826****，差额 **16,640** 与 **mixer 内 LSTM（66,560）与 GRU（49,920）** 的差一一对应；比较两个模块时，应聚焦这 **~16.6k** 参数差异及约 **1.165×** 的总参数量比例。
