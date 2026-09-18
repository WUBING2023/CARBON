"""Training utilities: padding, parameter count, and batch metrics."""

import random
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.stats import pearsonr as pearsonr_array
from sklearn import metrics
import torch

from CARBON.utils.utils import seed_everything
from CARBON.utils.functions import find_peaks


def randomPadN(df, NLen=101, column="seq", signal="signal"):
    processedData = []
    for index, peakElement in df.iterrows():
        tempSignal = np.array(eval(peakElement[signal]), dtype=float)
        if random.choice(["left", "right"]) == "right":
            leftPadN = random.randint(0, NLen // 2)
            peakElement[column] = "N" * leftPadN + peakElement[column][leftPadN:(len(peakElement[column]))]
            tempSignal[0:leftPadN] = 0
        else:
            rightPadN = random.randint(0, NLen // 2)
            peakElement[column] = peakElement[column][0:((len(peakElement[column]) - rightPadN))] + "N" * rightPadN
            tempSignal[(len(tempSignal) - rightPadN):len(tempSignal)] = 0
        peakElement[signal] = str(tempSignal.tolist())
        processedData.append(peakElement.tolist())
    return pd.DataFrame(processedData, columns=df.columns)


def get_parameter_number(model_analyse):
    total_num = sum(p.numel() for p in model_analyse.parameters())
    trainable_num = sum(p.numel() for p in model_analyse.parameters() if p.requires_grad)
    return f'Total parameters: {total_num}, Trainable parameters: {trainable_num}'


def get_used_parameters(model, *inputs):
    used_params = set()

    def hook_fn(module, input, output):
        for name, param in module.named_parameters(recurse=False):
            if param.requires_grad:
                used_params.add(param)

    hooks = []
    for layer in model.children():
        hooks.append(layer.register_forward_hook(hook_fn))

    with torch.no_grad():
        model(*inputs)

    total_params = sum(p.numel() for p in used_params)
    for hook in hooks:
        hook.remove()

    return total_params


def calculate_batch_metrics_seq(outputs_np, targets_np):
    roc_auc_values = []
    pr_auc_values = []
    sp_cor_values = []

    for row_predict, row_target in zip(outputs_np, targets_np):
        spearman_corr = spearmanr(row_predict, row_target).correlation
        signal_peak_mark, _ = find_peaks(row_target)
        fpr, tpr, thresholds = metrics.roc_curve(y_true=signal_peak_mark, y_score=row_predict, pos_label=1)
        roc_auc = metrics.auc(fpr, tpr)
        pr_auc = metrics.average_precision_score(y_true=signal_peak_mark, y_score=row_predict, pos_label=1)

        roc_auc_values.append(roc_auc)
        pr_auc_values.append(pr_auc)
        sp_cor_values.append(spearman_corr)

    return roc_auc_values, pr_auc_values, sp_cor_values


def calculate_batch_metrics_whole_auc(outputs, targets):
    outputs_flat = outputs.flatten()
    targets_flat = targets.flatten()

    sp_whole = spearmanr(outputs_flat, targets_flat)[0]
    p_whole = pearsonr_array(outputs_flat, targets_flat)[0]

    signal_peak_mark, peak_count = find_peaks(targets_flat, signal_max_threshold=2, signal_sum_threshold=2)
    fpr, tpr, thresholds = metrics.roc_curve(y_true=signal_peak_mark, y_score=outputs_flat, pos_label=1)
    roc_auc_whole = metrics.auc(fpr, tpr)
    pr_auc_whole = metrics.average_precision_score(y_true=signal_peak_mark, y_score=outputs_flat, pos_label=1)

    return roc_auc_whole, pr_auc_whole, sp_whole, p_whole


def calculate_batch_metrics_whole(outputs, targets):
    outputs_flat = outputs.flatten()
    targets_flat = targets.flatten()

    sp_whole = spearmanr(outputs_flat, targets_flat)[0]
    p_whole = pearsonr_array(outputs_flat, targets_flat)[0]

    peak_mark, peak_count = find_peaks(targets_flat, signal_max_threshold=2, signal_sum_threshold=2)
    n_peak_pos = int(peak_mark.sum())
    n_peak_neg = len(peak_mark) - n_peak_pos
    if n_peak_pos > 0 and n_peak_neg > 0:
        fpr_p, tpr_p, _ = metrics.roc_curve(y_true=peak_mark, y_score=outputs_flat, pos_label=1)
        peak_auc = float(metrics.auc(fpr_p, tpr_p))
        peak_ap = float(metrics.average_precision_score(y_true=peak_mark, y_score=outputs_flat, pos_label=1))
    else:
        peak_auc = peak_ap = float("nan")

    med = np.median(targets_flat)
    signal_mark = (targets_flat > med).astype(np.int32)
    n_sig_pos = int(signal_mark.sum())
    n_sig_neg = len(signal_mark) - n_sig_pos
    if n_sig_pos > 0 and n_sig_neg > 0:
        fpr_s, tpr_s, _ = metrics.roc_curve(y_true=signal_mark, y_score=outputs_flat, pos_label=1)
        signal_auc = float(metrics.auc(fpr_s, tpr_s))
        signal_ap = float(metrics.average_precision_score(y_true=signal_mark, y_score=outputs_flat, pos_label=1))
    else:
        signal_auc = signal_ap = float("nan")

    return signal_auc, signal_ap, peak_auc, peak_ap, sp_whole, p_whole
