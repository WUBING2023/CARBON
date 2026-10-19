"""
Compare MARCP total parameters: official BiGRU vs symmetric BiLSTM mixer.
Run from any cwd; adds ../.. / RBP/MARCP to sys.path for imports.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import torch
import torch.nn as nn
from einops import rearrange
from torchinfo import summary

# MARCP package lives in RBP/MARCP
_RBP_ROOT = Path(__file__).resolve().parent.parent.parent
_MARCP_ROOT = _RBP_ROOT / "MARCP"
if not _MARCP_ROOT.is_dir():
    raise FileNotFoundError(f"MARCP repo not found at {_MARCP_ROOT}")
sys.path.insert(0, str(_MARCP_ROOT))

from MARCP.models.marcp_model import (  # noqa: E402
    MARCP,
    MultiScaleMotifExtractor,
    Conv1dProfileHead,
)
from MARCP.operations.train_ops import get_parameter_number  # noqa: E402


class ModernBiLSTM(nn.Module):
    """Mirror of ModernBiGRU with nn.LSTM (same proj / residual / norm / dropout)."""

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


class MARCPBiLSTM(nn.Module):
    """Same as MARCP but mixer uses BiLSTM instead of BiGRU."""

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
    ):
        super().__init__()
        self.control = control
        self.in_channels = in_channels
        self.motif = MultiScaleMotifExtractor(
            in_channels,
            motif_out_per_branch,
            motif_kernels,
            ablation_no_se=ablation_no_se,
            ablation_no_multiscale=ablation_no_multiscale,
        )
        d = self.motif.out_channels
        self.mixer = ModernBiLSTM(d, hidden_size=gru_hidden, dropout=dropOut, ablation_no_gru=ablation_no_gru)
        head_in = gru_hidden * 2
        self.profile_head = Conv1dProfileHead(head_in, profile_hidden, dropOut)
        self.count_head = nn.Sequential(
            nn.Linear(head_in, 32),
            nn.GELU(),
            nn.Dropout(dropOut),
            nn.Linear(32, 1),
        )

    def _backbone(self, x: torch.Tensor) -> torch.Tensor:
        x = rearrange(x, "b n d -> b d n")
        x = self.motif(x)
        x = self.mixer(x)
        return x

    def forward(self, x: torch.Tensor) -> dict:
        x = self._backbone(x)
        logits_profile = self.profile_head(x).squeeze(1)
        x_pool = x.mean(dim=-1)
        log_total_count = self.count_head(x_pool)
        return {"logits_profile": logits_profile, "log_total_count": log_total_count}


class _SummaryWrapper(nn.Module):
    """torchinfo expects a tensor return; MARCP returns a dict."""

    def __init__(self, inner: nn.Module):
        super().__init__()
        self.inner = inner

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        out = self.inner(x)
        return out["logits_profile"]


def _total_trainable_params(model: nn.Module) -> tuple[int, int]:
    total = sum(p.numel() for p in model.parameters())
    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    return total, trainable


def build_common_kwargs(args: argparse.Namespace) -> dict:
    return {
        "dropOut": args.dropout,
        "control": args.control,
        "in_channels": args.in_channels,
        "motif_out_per_branch": args.motif_out_per_branch,
        "motif_kernels": tuple(args.motif_kernels),
        "gru_hidden": args.gru_hidden,
        "profile_hidden": tuple(args.profile_hidden),
        "ablation_no_se": args.ablation_no_se,
        "ablation_no_multiscale": args.ablation_no_multiscale,
        "ablation_no_gru": args.ablation_no_gru,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="MARCP BiGRU vs BiLSTM parameter count (PyTorch + torchinfo)")
    parser.add_argument("--batch_size", type=int, default=1, help="Dummy batch size for torchinfo")
    parser.add_argument("--seq_length", type=int, default=101, help="Sequence length L (B, L, C)")
    parser.add_argument("--in_channels", type=int, default=5)
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--control", action="store_true")
    parser.add_argument("--motif_out_per_branch", type=int, default=32)
    parser.add_argument("--motif_kernels", type=int, nargs="+", default=[3, 5, 7])
    parser.add_argument("--gru_hidden", type=int, default=64)
    parser.add_argument("--profile_hidden", type=int, nargs="+", default=[64, 32])
    parser.add_argument("--ablation_no_se", action="store_true")
    parser.add_argument("--ablation_no_multiscale", action="store_true")
    parser.add_argument("--ablation_no_gru", action="store_true")
    parser.add_argument("--summary_depth", type=int, default=5)
    args = parser.parse_args()

    if args.ablation_no_gru:
        print("[warn] ablation_no_gru=True removes both GRU and LSTM paths; recurrent param diff should be ~0.\n")

    kw = build_common_kwargs(args)
    model_gru = MARCP(**kw)
    model_lstm = MARCPBiLSTM(**kw)

    print("=== PyTorch parameter count (train_ops.get_parameter_number) ===")
    print("BiGRU (MARCP):", get_parameter_number(model_gru))
    print("BiLSTM:       ", get_parameter_number(model_lstm))

    n_gru, t_gru = _total_trainable_params(model_gru)
    n_lst, t_lst = _total_trainable_params(model_lstm)
    delta = n_lst - n_gru
    ratio = n_lst / n_gru if n_gru else float("nan")
    print("\n=== Total params (integers) ===")
    print(f"BiGRU total: {n_gru:,}  trainable: {t_gru:,}")
    print(f"BiLSTM total: {n_lst:,}  trainable: {t_lst:,}")
    print(f"Delta (LSTM - GRU): {delta:,}  ratio LSTM/GRU: {ratio:.6f}")

    bs, L, C = args.batch_size, args.seq_length, args.in_channels
    dummy = torch.randn(bs, L, C)

    print("\ntorchinfo.summary — BiGRU (MARCP)")
    print(
        summary(
            _SummaryWrapper(model_gru),
            input_size=(bs, L, C),
            depth=args.summary_depth,
            col_names=("output_size", "num_params"),
            verbose=0,
        )
    )

    print("\ntorchinfo.summary — BiLSTM (MARCPBiLSTM)")
    print(
        summary(
            _SummaryWrapper(model_lstm),
            input_size=(bs, L, C),
            depth=args.summary_depth,
            col_names=("output_size", "num_params"),
            verbose=0,
        )
    )


if __name__ == "__main__":
    main()
