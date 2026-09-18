"""
CARBON: RNA binding profile model (multi-scale motif + conservation).
Multi-scale motif + SE channel attention + BiGRU or BiLSTM; input (B, L, C) -> logits_profile (B, L), log_total_count (B, 1).
"""

from __future__ import absolute_import, division, print_function

import random
import os
import numpy as np
import torch
import torch.nn as nn
from einops import rearrange


def seed_everything(seed):
    if seed is None:
        seed = random.randint(1, 10000)
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.enabled = True


class MultiScaleMotifExtractor(nn.Module):
    """Parallel Conv1d branches (k=3,5,7), concat + SE channel attention, fuse + residual. Supports ablation no_se, no_multiscale."""

    def __init__(self, in_channels: int, out_per_branch: int = 32, kernels: tuple = (3, 5, 7), se_reduction: int = 4,
                 ablation_no_se: bool = False, ablation_no_multiscale: bool = False):
        super().__init__()
        self.ablation_no_se = ablation_no_se
        self.ablation_no_multiscale = ablation_no_multiscale
        if ablation_no_multiscale:
            kernels = (5,)
        self.branches = nn.ModuleList()
        for k in kernels:
            padding = k // 2
            self.branches.append(
                nn.Sequential(
                    nn.Conv1d(in_channels, out_per_branch, kernel_size=k, padding=padding),
                    nn.BatchNorm1d(out_per_branch),
                    nn.GELU(),
                )
            )
        self.out_channels = out_per_branch * len(self.branches)
        mid = max(self.out_channels // se_reduction, 8)
        if not ablation_no_se:
            self.se = nn.Sequential(
                nn.AdaptiveAvgPool1d(1),
                nn.Conv1d(self.out_channels, mid, 1),
                nn.ReLU(inplace=True),
                nn.Conv1d(mid, self.out_channels, 1),
                nn.Sigmoid(),
            )
        else:
            self.se = None
        self.fuse = nn.Conv1d(self.out_channels, self.out_channels, 1)
        self.shortcut = nn.Conv1d(in_channels, self.out_channels, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        shortcut = self.shortcut(x)
        outs = [b(x) for b in self.branches]
        out = torch.cat(outs, dim=1)
        if self.se is not None:
            se_weight = self.se(out)
            out = out * se_weight
        return self.fuse(out) + shortcut


class ModernBiGRU(nn.Module):
    """Single-layer BiGRU + residual + LayerNorm. Supports ablation no_gru (residual + LayerNorm only)."""

    def __init__(self, in_channels: int, hidden_size: int = 64, dropout: float = 0.2, ablation_no_gru: bool = False):
        super().__init__()
        self.hidden_size = hidden_size
        self.ablation_no_gru = ablation_no_gru
        if not ablation_no_gru:
            self.proj = nn.Conv1d(in_channels, hidden_size, 1)
            self.gru = nn.GRU(hidden_size, hidden_size, num_layers=1, batch_first=False, bidirectional=True, dropout=0)
        else:
            self.proj = None
            self.gru = None
        self.residual_proj = nn.Conv1d(in_channels, hidden_size * 2, 1)
        self.norm = nn.LayerNorm(hidden_size * 2)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B, D, L = x.shape
        residual = self.residual_proj(x)
        if self.gru is not None:
            x = self.proj(x)
            x = x.permute(2, 0, 1)
            y, _ = self.gru(x)
            y = y.permute(1, 2, 0)
        else:
            y = torch.zeros_like(residual, device=residual.device, dtype=residual.dtype)
        out = y + residual
        out = out.transpose(1, 2)
        out = self.norm(out)
        out = self.dropout(out)
        return out.transpose(1, 2)


class ModernBiLSTM(nn.Module):
    """Single-layer BiLSTM + residual + LayerNorm. Same layout as ModernBiGRU; supports ablation_no_gru (no RNN)."""

    def __init__(self, in_channels: int, hidden_size: int = 64, dropout: float = 0.2, ablation_no_gru: bool = False):
        super().__init__()
        self.hidden_size = hidden_size
        self.ablation_no_gru = ablation_no_gru
        if not ablation_no_gru:
            self.proj = nn.Conv1d(in_channels, hidden_size, 1)
            self.lstm = nn.LSTM(
                hidden_size, hidden_size, num_layers=1, batch_first=False, bidirectional=True, dropout=0
            )
        else:
            self.proj = None
            self.lstm = None
        self.residual_proj = nn.Conv1d(in_channels, hidden_size * 2, 1)
        self.norm = nn.LayerNorm(hidden_size * 2)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        residual = self.residual_proj(x)
        if self.lstm is not None:
            x = self.proj(x)
            x = x.permute(2, 0, 1)
            y, _ = self.lstm(x)
            y = y.permute(1, 2, 0)
        else:
            y = torch.zeros_like(residual, device=residual.device, dtype=residual.dtype)
        out = y + residual
        out = out.transpose(1, 2)
        out = self.norm(out)
        out = self.dropout(out)
        return out.transpose(1, 2)


class Conv1dProfileHead(nn.Module):
    """All Conv1d(kernel_size=1): (B, D, L) -> (B, 1, L), strict spatial alignment."""

    def __init__(self, in_channels: int, hidden: tuple = (64, 32), dropout: float = 0.2):
        super().__init__()
        layers = []
        c = in_channels
        for h in hidden:
            layers += [
                nn.Conv1d(c, h, 1),
                nn.BatchNorm1d(h),
                nn.GELU(),
                nn.Dropout(dropout),
            ]
            c = h
        layers.append(nn.Conv1d(c, 1, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class CARBON(nn.Module):
    """CARBON: multi-scale motif + SE + BiGRU or BiLSTM; outputs logits_profile and log_total_count. Supports ablations: no_se, no_multiscale, no_gru, lstm_instead_of_gru."""

    def __init__(
        self,
        dropOut: float = 0.2,
        control: bool = False,
        in_channels: int = 5,
        motif_out_per_branch: int = 32,
        motif_kernels: tuple = (3, 5, 7),
        gru_hidden: int = 64,
        profile_hidden: tuple = (64, 32),
        ablation_no_se: bool = False,
        ablation_no_multiscale: bool = False,
        ablation_no_gru: bool = False,
        ablation_lstm_instead_of_gru: bool = False,
    ):
        super().__init__()
        self.control = control
        self.in_channels = in_channels

        self.motif = MultiScaleMotifExtractor(
            in_channels, motif_out_per_branch, motif_kernels,
            ablation_no_se=ablation_no_se, ablation_no_multiscale=ablation_no_multiscale,
        )
        d = self.motif.out_channels
        if ablation_no_gru:
            self.mixer = ModernBiGRU(d, hidden_size=gru_hidden, dropout=dropOut, ablation_no_gru=True)
        elif ablation_lstm_instead_of_gru:
            self.mixer = ModernBiLSTM(d, hidden_size=gru_hidden, dropout=dropOut, ablation_no_gru=False)
        else:
            self.mixer = ModernBiGRU(d, hidden_size=gru_hidden, dropout=dropOut, ablation_no_gru=False)
        head_in = gru_hidden * 2
        self.profile_head = Conv1dProfileHead(head_in, profile_hidden, dropOut)
        self.count_head = nn.Sequential(
            nn.Linear(head_in, 32),
            nn.GELU(),
            nn.Dropout(dropOut),
            nn.Linear(32, 1),
        )

    def _backbone(self, x: torch.Tensor) -> torch.Tensor:
        x = rearrange(x, 'b n d -> b d n')
        x = self.motif(x)
        x = self.mixer(x)
        return x

    def forward(self, x: torch.Tensor) -> dict:
        x = self._backbone(x)
        logits_profile = self.profile_head(x).squeeze(1)
        x_pool = x.mean(dim=-1)
        log_total_count = self.count_head(x_pool)
        return {
            'logits_profile': logits_profile,
            'log_total_count': log_total_count,
        }
