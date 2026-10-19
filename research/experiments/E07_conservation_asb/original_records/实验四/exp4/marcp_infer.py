from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

from exp4.paths import marcp_root


def ensure_marcp_on_path() -> Path:
    root = marcp_root()
    s = str(root.resolve())
    if s not in sys.path:
        sys.path.insert(0, s)
    return root


def load_marcp(in_channels: int, weights: Path, device: torch.device) -> torch.nn.Module:
    ensure_marcp_on_path()
    from MARCP.models.marcp_model import MARCP

    model = MARCP(dropOut=0.2, control=False, in_channels=in_channels).to(device)
    try:
        state = torch.load(str(weights), map_location=device, weights_only=True)
    except TypeError:
        state = torch.load(str(weights), map_location=device)
    model.load_state_dict(state, strict=False)
    model.eval()
    return model


@torch.no_grad()
def predict_profile_batched(
    model: torch.nn.Module,
    x: np.ndarray,
    device: torch.device,
    *,
    batch_size: int = 4096,
) -> np.ndarray:
    """
    x: (N, L, C) float32，C=4 或 5。
    返回 (N, L) 预测 profile。
    """
    model.eval()
    n = int(x.shape[0])
    outs: list[np.ndarray] = []
    for i in range(0, n, batch_size):
        sl = slice(i, min(i + batch_size, n))
        t = torch.from_numpy(np.ascontiguousarray(x[sl])).float().to(device)
        o = model(t)
        lp = o["logits_profile"]
        ltc = o["log_total_count"]
        if ltc.dim() == 1:
            ltc = ltc.unsqueeze(-1)
        prob = F.softmax(lp, dim=-1)
        pred = prob * torch.exp(ltc)
        if pred.dim() == 3 and pred.shape[-1] == 1:
            pred = pred.squeeze(-1)
        outs.append(pred.detach().cpu().numpy().astype(np.float32, copy=False))
    return np.concatenate(outs, axis=0)
