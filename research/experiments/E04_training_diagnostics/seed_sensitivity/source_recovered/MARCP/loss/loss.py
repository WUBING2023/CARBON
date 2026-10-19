"""Profile and count loss: CombinedProfileLoss (shape + jitter + count) and legacy loss_fun."""

import torch
import torch.nn.functional as F
from torch.distributions.multinomial import Multinomial
import torch.nn as nn
import math

try:
    from torcheval.metrics.functional import r2_score
except ImportError:
    def r2_score(input: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        ss_res = ((target - input).pow(2)).sum()
        ss_tot = ((target - target.mean()).pow(2)).sum().clamp(min=1e-8)
        return (1 - ss_res / ss_tot).unsqueeze(0)


class CombinedProfileLoss(nn.Module):
    """
    Combined loss: profile shape (Multinomial NLL), jitter-tolerant (CDF MSE), total count (MSE on log).
    Weights: alpha (shape), beta (jitter), gamma (count); total scaled by loss_scale.
    """

    def __init__(self, alpha: float = 1.0, beta: float = 1.0, gamma: float = 1.0, loss_scale: float = 1.0):
        super().__init__()
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.loss_scale = loss_scale

    def profile_shape_loss(self, logits_profile: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        log_probs = F.log_softmax(logits_profile, dim=-1)
        loss = -torch.sum(y_true * log_probs, dim=-1).mean()
        return loss

    def jitter_loss(self, logits_profile: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        prob_pred = F.softmax(logits_profile, dim=-1)
        y_sum = y_true.sum(dim=-1, keepdim=True)
        prob_true = y_true / (y_sum + 1e-6)
        cdf_pred = torch.cumsum(prob_pred, dim=-1)
        cdf_true = torch.cumsum(prob_true, dim=-1)
        loss = F.mse_loss(cdf_pred, cdf_true)
        return loss

    def count_loss(self, log_total_count: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
        true_total = y_true.sum(dim=-1, keepdim=True)
        log_true_total = torch.log(true_total + 1.0)
        loss = F.mse_loss(log_total_count, log_true_total)
        return loss

    def forward(
        self,
        logits_profile: torch.Tensor,
        log_total_count: torch.Tensor,
        y_true: torch.Tensor,
    ) -> dict:
        l_shape = self.profile_shape_loss(logits_profile, y_true)
        l_jitter = self.jitter_loss(logits_profile, y_true)
        l_count = self.count_loss(log_total_count, y_true)
        raw_total = self.alpha * l_shape + self.beta * l_jitter + self.gamma * l_count
        total = raw_total / self.loss_scale
        return {
            'loss': total,
            'loss_shape': l_shape,
            'loss_jitter': l_jitter,
            'loss_count': l_count,
        }


class loss_fun(object):
    def __init__(self):
        self.eps = 1e-8

    def kl_divergence_loss(self, prediction, target):
        kl_loss = torch.nn.KLDivLoss(reduction="batchmean", log_target=True)
        log_prediction = torch.nn.functional.log_softmax(prediction, dim=1)
        log_target = torch.nn.functional.log_softmax(target, dim=1)
        loss = kl_loss(log_prediction, log_target)
        return loss

    def mse(self, prediction, target):
        loss = torch.nn.MSELoss()(prediction, target)
        return loss

    def SmoothL1Loss(self, prediction, target):
        loss = torch.nn.SmoothL1Loss()(prediction, target)
        return loss

    def multinomial_ng(self, prediction, target):
        total_counts = int(torch.sum(target).item())
        dist = Multinomial(total_count=total_counts, logits=prediction)
        loss = -torch.mean(dist.log_prob(target))
        return loss

    def multinomialnll_batch(self, prediction, target):
        total_counts = int(torch.sum(target).item())
        dist = Multinomial(total_count=total_counts, logits=prediction)
        batchLen = torch.tensor(target.shape[0], dtype=torch.float32)
        loss = -torch.sum(dist.log_prob(target)) / batchLen
        return loss

    def multinomialnll_seq(self, prediction, target):
        total_counts = int(torch.sum(target).item())
        dist = Multinomial(total_count=total_counts, logits=prediction)
        seqLen = torch.tensor(target.shape[1], dtype=torch.float32)
        loss = -torch.sum(dist.log_prob(target)) / seqLen
        return loss

    @staticmethod
    def __log(t, eps=1e-20):
        return torch.log(t.clamp(min=eps))

    def poissonLoss(self, prediction, target):
        return (prediction - target * self.__log(prediction)).mean()

    def PoissonNLLLoss(self, prediction, target):
        loss = nn.PoissonNLLLoss(log_input=False, eps=1e-8)(prediction, target)
        return loss

    def mse_poisson(self, prediction, target):
        mse_loss = self.mse(prediction, target)
        poiss_loss = self.poissonLoss(prediction, target)
        total_loss = (2 * mse_loss * poiss_loss) / (mse_loss + poiss_loss)
        return total_loss

    def mse_poisson_alpha(self, prediction, target, alpha=1):
        mse_loss = self.mse(prediction, target)
        poiss_loss = self.poissonLoss(prediction, target)
        total_loss = alpha * poiss_loss + mse_loss
        return total_loss

    def mse_poisson_Mean(self, prediction, target):
        mse_loss = self.mse(prediction, target)
        poiss_loss = self.PoissonNLLLoss(prediction, target)
        total_loss = (2 * mse_loss * poiss_loss) / (mse_loss + poiss_loss)
        return total_loss

    def multinomialnll_mse(self, y_pred, y_true, alpha=1):
        mult_loss = self.multinomialnll_batch(y_pred, y_true)
        mse_loss = nn.MSELoss()(y_pred, y_true)
        total_loss = alpha * mult_loss + mse_loss
        return total_loss

    def l1_loss(self, y_pred, y_true):
        return nn.L1Loss()(y_pred, y_true)
