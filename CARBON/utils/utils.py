import torch
import os
import random
import numpy as np


def pearsonr_tensor(input, target, batch_first=True):
    """
    Per-sample Pearson correlation, compatible with audtorch interface.
    input/target: (B, L) when batch_first=True; returns (B,) tensor, call .nanmean() for scalar.
    """
    if not batch_first:
        input = input.t()
        target = target.t()
    B = input.shape[0]
    device = input.device
    out = []
    for i in range(B):
        a = input[i].flatten().float()
        b = target[i].flatten().float()
        n = a.numel()
        if n < 2:
            out.append(torch.tensor(float("nan"), device=device))
            continue
        am, bm = a.mean(), b.mean()
        cov = ((a - am) * (b - bm)).sum() / (n - 1)
        std_a = ((a - am).pow(2).sum() / (n - 1)).clamp(min=1e-12).sqrt()
        std_b = ((b - bm).pow(2).sum() / (n - 1)).clamp(min=1e-12).sqrt()
        r = (cov / (std_a * std_b)).clamp(-1.0, 1.0)
        out.append(r)
    return torch.stack(out)


class epochAverageMeter(object):
    def __init__(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count


class Meter(object):
    def __init__(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count


def save_model(path, model, optimizer, name="model"):
    if not os.path.exists(path):
        os.mkdir(path)
    torch.save(model, f"{path}/{name}.pth")
    torch.save(model.state_dict(), f'{path}/{name}_params.pth')
    torch.save(optimizer.state_dict(), f'{path}/{name}_optimizer.pth')


def seed_everything(seed):
    """Set all relevant RNG seeds for reproducibility."""
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


def transfer_state_dict(pretrained_dict, model_dict):
    """
    Filter pretrained_dict by model_dict keys and shapes for transfer learning.
    """
    state_dict = {}
    for k, v in pretrained_dict.items():
        if k in model_dict.keys():
            if model_dict[k].shape == v.shape:
                state_dict[k] = v
            else:
                shape0, shape1 = model_dict[k].shape
                state_dict[k] = v[:shape0, :shape1]
    return state_dict


def transfer_model(pretrained_model, model):
    """Load matching parameters from pretrained_model into model."""
    pretrained_dict = pretrained_model.state_dict()
    model_dict = model.state_dict()
    pretrained_dict = transfer_state_dict(pretrained_dict, model_dict)
    model_dict.update(pretrained_dict)
    model.load_state_dict(model_dict)
    return model
