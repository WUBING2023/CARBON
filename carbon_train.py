"""
Training and evaluation script for the CARBON model.
Supports Treat/Control frameworks, CombinedProfileLoss, 4ch/5ch input, and ablation options.
"""

import os
import sys
import argparse
import gc
import time
from collections import Counter

import h5py
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from torch.utils.tensorboard import SummaryWriter
from tqdm import tqdm

sys.path.append('./')
try:
    from audtorch.metrics.functional import pearsonr as pearsonr_tensor
except ImportError:
    from CARBON.utils.utils import pearsonr_tensor
from CARBON.utils.utils import (
    epochAverageMeter,
    save_model,
    seed_everything,
)
from CARBON.utils.functions import onehot_encode
from CARBON.loss.loss import loss_fun, CombinedProfileLoss
from CARBON.operations.train_ops import (
    randomPadN,
    get_parameter_number,
    calculate_batch_metrics_seq,
    calculate_batch_metrics_whole,
)
from CARBON.models.carbon_model import CARBON

def parse_arguments():
    """Parse command-line arguments for training and evaluation."""
    parser = argparse.ArgumentParser(
        prog='carbon_train.py',
        description='Training script for CARBON model',
        epilog="Show help message."
    )
    
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    parser.add_argument("--bsz", type=int, default=256, help="Batch size")
    parser.add_argument("--epochs", type=int, default=100, help="Number of training epochs")
    
    parser.add_argument("--lr_decay_count_max", type=int, default=5, 
                       help="Max count for learning rate decay")
    parser.add_argument("--max_stopping_step", type=int, default=9, 
                       help="Early stopping patience")
    parser.add_argument("--lr", type=float, default=0.000785, help="Initial learning rate")
    parser.add_argument("--lr_min", type=float, default=1e-7, 
                       help="Minimum learning rate")
    
    parser.add_argument("--pretrained_model_path", type=str, 
                       help="Path to pretrained model")
    parser.add_argument("--output", dest="output", type=str, required=True,
                       help="Output directory path")
    parser.add_argument("--task", dest="task", type=str)
    parser.add_argument("--CUDA", type=str, default="0", 
                       help="CUDA device ID(s), e.g. 0,1,2,3 for multi-GPU")
    parser.add_argument("--seqLength", type=int, default=101, 
                       help="Sequence length (default 101 for 5ch data)")
    parser.add_argument("--randomPadNLen", type=int, default=101,
                       help="Random padding length")
    parser.add_argument("--encode", type=str, default="OneHot",
                       choices=["OneHot", "Embedding"], 
                       help="Sequence encoding method (default OneHot for 5ch data)")
    parser.add_argument("--lossM", type=str, dest="loss_fun_name",
                       required=True, help="Loss function name")
    parser.add_argument("--eval_only", action="store_true",
                       help="Only run test set evaluation (load model from --load_model)")
    parser.add_argument("--load_model", type=str, default=None,
                       help="Path to model_params.pth for --eval_only")
    parser.add_argument("--resume", action="store_true",
                       help="Resume training from --output (requires model_params.pth and checkpoint_state.txt)")
    parser.add_argument('-train', '--trainFile', required=False,
                       metavar="trainFile", type=str,
                       help='Path to training data file')
    parser.add_argument('-validation', '--validationFile', required=False,
                       metavar="validationFile", type=str,
                       help='Path to validation data file')
    parser.add_argument('-test', '--testFile', type=str, default=None,
                       metavar="testFile",
                       help='Path to test data file (optional, for final evaluation)')
    parser.add_argument('--log', action="store_true",
                       help="Apply log transformation to data")
    parser.add_argument('--modelFramework', required=True,
                       choices=["Treat", "Control"],
                       help="Model framework type")
    parser.add_argument('--shuff', action="store_true",
                       help="Enable sequence shuffling augmentation")
    parser.add_argument("--signalMax", type=int, default=0,
                       help="Maximum signal threshold for filtering")
    parser.add_argument('--filter', action="store_true",
                       help="Enable data filtering")
    parser.add_argument("--input_channels", type=int, default=5,
                       help="Input feature channels: 4=OneHot only, 5=OneHot+PhyloP (default 5 for _5ch.h5 data)")

    parser.add_argument("--alpha", type=float, default=1.42,
                       help="Weight for Profile Shape Loss (Multinomial NLL)")
    parser.add_argument("--beta", type=float, default=0.54,
                       help="Weight for Jitter-tolerant Loss (CDF-Wasserstein)")
    parser.add_argument("--gamma", type=float, default=1.80,
                       help="Weight for Total Count Loss (MSE on log-counts)")
    parser.add_argument("--loss_scale", type=float, default=1.0,
                       help="Divide total loss by this factor to match poissonLoss magnitude")
    parser.add_argument("--weight_decay", type=float, default=0.000114,
                       help="Weight decay for AdamW (0 to use Adam)")
    parser.add_argument("--grad_clip", type=float, default=2.0,
                       help="Gradient clipping max norm (0 to disable)")
    parser.add_argument("--use_cosine_lr", action="store_true",
                       help="Use CosineAnnealingLR (T_max=epochs, eta_min=1e-5)")
    parser.add_argument("--use_reduce_on_plateau", action="store_true",
                       help="Use ReduceLROnPlateau on val_prAUC (mode=max), factor=0.2, patience=5, min_lr=1e-6")
    parser.add_argument("--plateau_factor", type=float, default=0.2,
                       help="ReduceLROnPlateau factor (lr *= factor when no improvement)")
    parser.add_argument("--plateau_patience", type=int, default=5,
                       help="ReduceLROnPlateau patience (must be < max_stopping_step)")
    parser.add_argument("--plateau_min_lr", type=float, default=1e-6,
                       help="ReduceLROnPlateau minimum lr")

    parser.add_argument("--ablation_no_se", action="store_true", help="Ablation: remove SE channel attention")
    parser.add_argument("--ablation_no_multiscale", action="store_true", help="Ablation: single-scale (k=5 only)")
    parser.add_argument("--ablation_no_gru", action="store_true", help="Ablation: remove BiGRU (residual+LayerNorm only)")
    parser.add_argument("--ablation_lstm_instead_of_gru", action="store_true",
                       help="Use BiLSTM instead of BiGRU in the sequence mixer (mutually exclusive with --ablation_no_gru)")
    parser.add_argument("--ablation_no_combined_loss", action="store_true",
                       help="Ablation: use baseline loss instead of CombinedProfileLoss")
    parser.add_argument("--ablation_no_conservative", action="store_true",
                       help="Ablation: drop 5th (conservation) channel; 5ch input but in_channels=4")
    parser.add_argument("--ablation_baseline_loss", type=str, default="poissonLoss",
                       help="Baseline loss name when using --ablation_no_combined_loss (default: poissonLoss)")
    parser.add_argument("--ablation_suffix", type=str, default="", help="Optional suffix for output subdir")
    parser.add_argument("--resource_fraction", type=float, default=0.9, help="Resource fraction for run_info")
    parser.add_argument("--num_workers", type=int, default=0, help="DataLoader num_workers (0 = main process)")
    parser.add_argument("--prefetch_factor", type=int, default=2, help="DataLoader prefetch_factor")

    return parser.parse_args()


def _get_ablation_kw(args):
    """Return ablation kwargs for CARBON (no in_channels; no_conservative is handled via effective_in_channels)."""
    return {
        "ablation_no_se": getattr(args, "ablation_no_se", False),
        "ablation_no_multiscale": getattr(args, "ablation_no_multiscale", False),
        "ablation_no_gru": getattr(args, "ablation_no_gru", False),
        "ablation_lstm_instead_of_gru": getattr(args, "ablation_lstm_instead_of_gru", False),
    }


def _get_ablation_subdir(args):
    """Build output subdir name from ablation flags; returns None if no ablations."""
    parts = []
    if getattr(args, "ablation_no_se", False):
        parts.append("no_se")
    if getattr(args, "ablation_no_multiscale", False):
        parts.append("no_multiscale")
    if getattr(args, "ablation_no_gru", False):
        parts.append("no_gru")
    if getattr(args, "ablation_no_combined_loss", False):
        parts.append("no_combined_loss")
    if getattr(args, "ablation_no_conservative", False):
        parts.append("no_conservative")
    if getattr(args, "ablation_lstm_instead_of_gru", False):
        parts.append("lstm_instead_of_gru")
    if not parts:
        return None
    subdir = "_".join(parts)
    suffix = getattr(args, "ablation_suffix", "") or ""
    if suffix:
        subdir = subdir + "_" + suffix
    return subdir

def _read_data_single(bed_fa_name, signal_max, filter_flag, EncodeMode):
    """Load a single h5 file (key='df') and return encoded data and signal lists."""
    bed_fa = pd.read_hdf(bed_fa_name, key='df')
    bed_fa.columns = ["chr", "chromStart", "chromEnd", "peak", "score", "strand", "signal", "signal_ctl", "seq"]
    if filter_flag and signal_max > 0:
        bed_fa["mainPeakMaxSignal"] = np.max(np.array(bed_fa["signal"].tolist()), axis=1)
        bed_fa = bed_fa[bed_fa['mainPeakMaxSignal'] >= signal_max]
    bed_fa = bed_fa[~bed_fa['seq'].str.contains('N', case=False)]
    bed_fa = bed_fa.dropna()
    seqs = bed_fa["seq"].str.upper().str.replace('T', 'U').tolist()
    if EncodeMode == "OneHot":
        vocab = list("AUGC")
        train_data = onehot_encode(seqs, vocab, 4)
    else:
        train_data = None
    return train_data, bed_fa["signal"].tolist(), bed_fa["signal_ctl"].tolist()

def read_data(bed_fa_name, signal_max, filter_flag, EncodeMode):
    """Load and preprocess data; supports single path or comma-separated list of paths."""
    try:
        paths = [p.strip() for p in str(bed_fa_name).split(',') if p.strip()]
        all_data, all_signal, all_signal_ctl = [], [], []
        for p in paths:
            d, s, s_ctl = _read_data_single(p, signal_max, filter_flag, EncodeMode)
            all_data.append(np.array(d))
            all_signal.extend(s)
            all_signal_ctl.extend(s_ctl)
        train_data = np.concatenate(all_data, axis=0)
        print(f"Loaded {len(paths)} file(s), {len(all_signal)} sequences total")
        return train_data, all_signal, all_signal_ctl
    except Exception as e:
        print(f"Error loading data: {str(e)}")
        raise


def read_data_5ch(h5_path, signal_max, filter_flag):
    """
    Load 5-channel h5 data (OneHot + PhyloP).
    Supports single path or comma-separated paths; multiple files concatenated on axis=0.
    """
    raw = str(h5_path).replace('\n', ',').replace('\r', ',')
    paths = [p.strip() for p in raw.split(',') if p.strip()]
    all_X, all_Y, all_Y_ctl = [], [], []
    for i, p in enumerate(paths):
        if not os.path.isfile(p):
            raise FileNotFoundError(
                f"read_data_5ch: path {i+1}/{len(paths)} missing or unreadable: {p!r}"
            )
    for p in paths:
        with h5py.File(p, 'r') as hf:
            X = np.array(hf['X'], dtype=np.float32)
            Y = np.array(hf['Y'], dtype=np.float32)
            Y_ctl = np.array(hf['Y_ctl'], dtype=np.float32)
        n = X.shape[0]
        if filter_flag and signal_max > 0:
            main_peak_max = np.max(Y, axis=1)
            keep = main_peak_max >= signal_max
            X, Y, Y_ctl = X[keep], Y[keep], Y_ctl[keep]
        all_X.append(X)
        all_Y.append(Y)
        all_Y_ctl.append(Y_ctl)
    X = np.concatenate(all_X, axis=0)
    Y = np.concatenate(all_Y, axis=0)
    Y_ctl = np.concatenate(all_Y_ctl, axis=0)
    n = X.shape[0]
    signal_list = [Y[i] for i in range(n)]
    signal_ctl_list = [Y_ctl[i] for i in range(n)]
    print(f"5ch data: {len(paths)} file(s), X shape {X.shape}, {len(signal_list)} signals")
    return X, signal_list, signal_ctl_list


if __name__ == "__main__":
    seed_everything(0)
    args = parse_arguments()
    if getattr(args, "ablation_no_gru", False) and getattr(args, "ablation_lstm_instead_of_gru", False):
        raise SystemExit("--ablation_no_gru and --ablation_lstm_instead_of_gru are mutually exclusive")
    if getattr(args, 'eval_only', False):
        if not args.load_model or not args.testFile or not args.output or not args.task:
            raise SystemExit("--eval_only requires --load_model, -test, --output, --task")
    else:
        if not args.trainFile or not args.validationFile:
            raise SystemExit("Training requires -train and -validation")

    task = args.task
    protein,cell,accession = task.split("_")[0], task.split("_")[1],task.split("_")[2]

    # Set encode mode and data type
    EncodeMode = args.encode

    # Learning rate
    lr = args.lr

    # Filter or not
    filter = args.filter

    # Output directory; append ablation subdir and write ablation_config.txt if needed
    outputDir = args.output
    ablation_subdir = _get_ablation_subdir(args)
    if ablation_subdir:
        outputDir = os.path.join(args.output, ablation_subdir)

    # Device; DataParallel when multiple GPUs (e.g. --CUDA 0,1,2,3)
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.CUDA).strip()
    args.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {args.device}")
    args.multi_gpu = torch.cuda.is_available() and torch.cuda.device_count() > 1
    if args.multi_gpu:
        print(f"Multi-GPU: DataParallel on {torch.cuda.device_count()} GPU(s)")

    def create_folder_if_not_exists(folder_path):
        """Create directory if it does not exist."""
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Folder '{folder_path}' created.")
        else:
            print(f"Folder '{folder_path}' already exists.")
    
    # Set output directory
    create_folder_if_not_exists(outputDir)
    if ablation_subdir:
        ablation_config_path = os.path.join(outputDir, "ablation_config.txt")
        with open(ablation_config_path, "w", encoding="utf-8") as f:
            f.write(f"ablation_subdir\t{ablation_subdir}\n")
            for k in ["ablation_no_se", "ablation_no_multiscale", "ablation_no_gru",
                      "ablation_lstm_instead_of_gru", "ablation_no_combined_loss", "ablation_no_conservative"]:
                v = 1 if getattr(args, k, False) else 0
                f.write(f"{k}\t{v}\n")
            if getattr(args, "ablation_no_combined_loss", False):
                f.write(f"ablation_baseline_loss\t{getattr(args, 'ablation_baseline_loss', 'poissonLoss')}\n")
        print(f"Ablation config written to {ablation_config_path}")

    # Effective input channels (4 when no_conservative)
    effective_in_channels = 4 if getattr(args, "ablation_no_conservative", False) else args.input_channels

    def _test_files_exist(test_file_arg):
        if not test_file_arg or not str(test_file_arg).strip():
            return False
        return any(os.path.isfile(p.strip()) for p in str(test_file_arg).split(',') if p.strip())

    if getattr(args, 'eval_only', False):
        if not _test_files_exist(args.testFile):
            raise SystemExit(f"Test file missing or unreadable: {args.testFile}")
        create_folder_if_not_exists(outputDir + '/tensorboard')
        writer = SummaryWriter(outputDir + '/tensorboard')
        from CARBON.models.carbon_model import CARBON
        eff_ch = 4 if getattr(args, "ablation_no_conservative", False) else args.input_channels
        ab_kw = _get_ablation_kw(args)
        if args.modelFramework == "Treat":
            model = CARBON(dropOut=0.2, control=False, in_channels=eff_ch, **ab_kw).to(args.device)
        else:
            model = CARBON(dropOut=0.2, control=True, in_channels=eff_ch, **ab_kw).to(args.device)
        if not os.path.isfile(args.load_model):
            raise SystemExit(f"Model file not found: {args.load_model}")
        try:
            state = torch.load(args.load_model, map_location=args.device, weights_only=True)
        except TypeError:
            state = torch.load(args.load_model, map_location=args.device)
        model.load_state_dict(state)
        model.eval()
        if args.input_channels == 5:
            test_data, test_label, test_label_ctl = read_data_5ch(args.testFile, args.signalMax, args.filter)
            if getattr(args, "ablation_no_conservative", False):
                test_data = test_data[..., :4]
        else:
            test_data, test_label, test_label_ctl = read_data(args.testFile, args.signalMax, args.filter, EncodeMode)
        if args.log:
            test_label = np.log(np.asarray(test_label) + 1)
            test_label_ctl = np.log(np.asarray(test_label_ctl) + 1)
        else:
            test_label = np.asarray(test_label)
            test_label_ctl = np.asarray(test_label_ctl)
        test_data = torch.from_numpy(np.asarray(test_data))
        test_label = torch.from_numpy(test_label).float()
        test_label_ctl = torch.from_numpy(test_label_ctl).float()
        test_dataset = torch.utils.data.TensorDataset(test_data, test_label, test_label_ctl)
        testLoader = torch.utils.data.DataLoader(test_dataset, batch_size=args.bsz, shuffle=False, drop_last=False)
        del test_data, test_label, test_label_ctl
        gc.collect()
        test_start_time = time.perf_counter()
        all_preds, all_targets = [], []
        with torch.no_grad():
            for (inputs, targets, targets_ctl) in testLoader:
                inputs = inputs.to(args.device)
                targets = targets.to(args.device)
                targets_ctl = targets_ctl.to(args.device)
                if EncodeMode == "OneHot" or args.input_channels == 5:
                    outputs = model(inputs.float())
                if args.modelFramework == "Control":
                    w = outputs[2].unsqueeze(1)
                    target_predict, ctl_predict = outputs[0], outputs[1]
                    if args.log:
                        target_predict = torch.exp(target_predict) - 1
                        ctl_predict = torch.exp(ctl_predict) - 1
                        targets = torch.exp(targets) - 1
                    treat_predict = (target_predict * w) + (1 - w) * ctl_predict
                else:
                    if args.log:
                        targets = torch.exp(targets) - 1
                    logits_profile = outputs['logits_profile']
                    log_total_count = outputs['log_total_count']
                    prob = torch.softmax(logits_profile, dim=-1)
                    treat_predict = prob * torch.exp(log_total_count)
                all_preds.append(treat_predict.cpu().numpy())
                all_targets.append(targets.cpu().numpy())
        test_duration_sec = time.perf_counter() - test_start_time
        outputs_np = np.concatenate(all_preds, axis=0)
        targets_np = np.concatenate(all_targets, axis=0)
        signal_auc, signal_ap, peak_auc, peak_ap, spearman, mean_pcc = calculate_batch_metrics_whole(outputs_np, targets_np)
        mean_pcc_val = "nan" if np.isnan(mean_pcc) else float(mean_pcc)
        run_info_path = f"{outputDir}/run_info.txt"
        train_duration_sec = 0.0
        if os.path.isfile(run_info_path):
            with open(run_info_path, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("training_time_sec\t"):
                        try:
                            train_duration_sec = float(line.strip().split("\t")[1])
                        except Exception:
                            pass
                        break
        with open(f"{outputDir}/test_metrics_summary.txt", "w", encoding="utf-8") as f:
            f.write("Metric\tValue\n")
            f.write(f"Mean_PCC\t{mean_pcc_val}\n")
            f.write(f"Spearman\t{spearman}\n")
            f.write(f"Signal_AUC\t{signal_auc}\n")
            f.write(f"Signal_AP\t{signal_ap}\n")
            f.write(f"Peak_AUC\t{peak_auc}\n")
            f.write(f"Peak_AP\t{peak_ap}\n")
            f.write(f"Training_time_sec\t{train_duration_sec}\n")
            f.write(f"Testing_time_sec\t{test_duration_sec}\n")
        with open(f"{outputDir}/test_metric.txt", "w", encoding="utf-8") as f:
            f.write("Mean_PCC\tSpearman\tSignal_AUC\tSignal_AP\tPeak_AUC\tPeak_AP\tTraining_time_sec\tTesting_time_sec\n")
            f.write(f"{mean_pcc_val}\t{spearman}\t{signal_auc}\t{signal_ap}\t{peak_auc}\t{peak_ap}\t{train_duration_sec}\t{test_duration_sec}\n")
        with open(run_info_path, "a", encoding="utf-8") as f:
            f.write(f"testing_time_sec\t{test_duration_sec}\n")
            f.write(f"testing_time_min\t{test_duration_sec/60}\n")
        print(f"\n=== eval_only done; metrics written to {outputDir} ===")
        sys.exit(0)

    create_folder_if_not_exists(outputDir + '/tensorboard')
    writer = SummaryWriter(outputDir + '/tensorboard')

    # Train data
    if args.input_channels == 5:
        train_data, train_label, train_label_ctl = read_data_5ch(args.trainFile, args.signalMax, args.filter)
        if getattr(args, "ablation_no_conservative", False):
            train_data = train_data[..., :4]
    else:
        train_data, train_label, train_label_ctl = read_data(args.trainFile, args.signalMax, args.filter, EncodeMode)

    # Apply log transform if specified
    if args.log:
        print("Applying log transform to training data")
        train_label = np.log(np.asarray(train_label) + 1)
        train_label_ctl = np.log(np.asarray(train_label_ctl) + 1)
    else:
        print("No log transform applied to training data")
        train_label = np.asarray(train_label)
        train_label_ctl = np.asarray(train_label_ctl)

    # Convert to PyTorch tensors
    train_data = torch.from_numpy(np.asarray(train_data))
    train_label = torch.from_numpy(train_label).float()
    train_label_ctl = torch.from_numpy(train_label_ctl).float()

    print(f"Training data shape: {train_data.shape}")
    print(f"Training label shape: {train_label.shape}")
    print(f"Training control label shape: {train_label_ctl.shape}")

    # Create DataLoader for training data
    train_dataset = torch.utils.data.TensorDataset(train_data, train_label, train_label_ctl)
    del train_data, train_label, train_label_ctl
    gc.collect()

    trainLoader = torch.utils.data.DataLoader(
        train_dataset, batch_size=args.bsz, shuffle=True, drop_last=True,
        num_workers=getattr(args, "num_workers", 0),
        prefetch_factor=getattr(args, "prefetch_factor", 2) if getattr(args, "num_workers", 0) > 0 else None,
    )
    del train_dataset
    gc.collect()

    # Validation data
    if args.input_channels == 5:
        val_data, val_label, val_label_ctl = read_data_5ch(args.validationFile, args.signalMax, args.filter)
        if getattr(args, "ablation_no_conservative", False):
            val_data = val_data[..., :4]
    else:
        val_data, val_label, val_label_ctl = read_data(args.validationFile, args.signalMax, args.filter, EncodeMode)

    # Apply log transform if specified
    if args.log:
        print("Applying log transform to validation data")
        val_label = np.log(np.asarray(val_label) + 1)
        val_label_ctl = np.log(np.asarray(val_label_ctl) + 1)
    else:
        print("No log transform applied to validation data")
        val_label = np.asarray(val_label)
        val_label_ctl = np.asarray(val_label_ctl)

    # Convert to PyTorch tensors
    val_data = torch.from_numpy(np.asarray(val_data))
    val_label = torch.from_numpy(val_label).float()
    val_label_ctl = torch.from_numpy(val_label_ctl).float()

    print(f"Validation data shape: {val_data.shape}")
    print(f"Validation label shape: {val_label.shape}")
    print(f"Validation control label shape: {val_label_ctl.shape}")

    val_dataset = torch.utils.data.TensorDataset(val_data, val_label, val_label_ctl)
    del val_data, val_label, val_label_ctl
    gc.collect()

    valLoader = torch.utils.data.DataLoader(
        val_dataset, batch_size=args.bsz, shuffle=True, drop_last=True,
        num_workers=getattr(args, "num_workers", 0),
        prefetch_factor=getattr(args, "prefetch_factor", 2) if getattr(args, "num_workers", 0) > 0 else None,
    )
    del val_dataset
    gc.collect()

    print(f"Initializing {args.modelFramework} model...")
    from CARBON.models.carbon_model import CARBON

    ab_kw = _get_ablation_kw(args)
    if args.modelFramework == "Treat":
        model = CARBON(dropOut=0.2, control=False, in_channels=effective_in_channels, **ab_kw).to(args.device)
    elif args.modelFramework == "Control":
        model = CARBON(dropOut=0.2, control=True, in_channels=effective_in_channels, **ab_kw).to(args.device)

    if getattr(args, 'multi_gpu', False):
        model = torch.nn.DataParallel(model)

    print(f"Model parameter count: {get_parameter_number(model)}")
    
    # CombinedProfileLoss (Treat only when --lossM CombinedProfileLoss)
    _use_combined_loss = not getattr(args, "ablation_no_combined_loss", False)
    _effective_baseline_loss = getattr(args, "ablation_baseline_loss", "poissonLoss")
    if args.modelFramework == "Treat" and args.loss_fun_name == "CombinedProfileLoss":
        criterion = CombinedProfileLoss(
            alpha=args.alpha, beta=args.beta, gamma=args.gamma,
            loss_scale=args.loss_scale
        ).to(args.device)
        if _use_combined_loss:
            print(f"Using CombinedProfileLoss: alpha={args.alpha}, beta={args.beta}, gamma={args.gamma}, loss_scale={args.loss_scale}")
        else:
            print(f"Ablation: using baseline loss '{_effective_baseline_loss}' instead of CombinedProfileLoss")
    
    # Initialize optimizer (AdamW with weight_decay for better generalization)
    if getattr(args, 'weight_decay', 0) > 0:
        optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
        print(f"Using AdamW, lr={args.lr}, weight_decay={args.weight_decay}")
    else:
        optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
        print(f"Using Adam, lr={args.lr}")

    if getattr(args, 'use_cosine_lr', False):
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs, eta_min=1e-5)
        print(f"Using CosineAnnealingLR, T_max={args.epochs}, eta_min=1e-5")
    elif getattr(args, 'use_reduce_on_plateau', False):
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='max', factor=getattr(args, 'plateau_factor', 0.2),
            patience=getattr(args, 'plateau_patience', 5), min_lr=getattr(args, 'plateau_min_lr', 1e-6)
        )
        print(f"Using ReduceLROnPlateau on val_prAUC (mode=max), factor={getattr(args, 'plateau_factor', 0.2)}, patience={getattr(args, 'plateau_patience', 5)}, min_lr={getattr(args, 'plateau_min_lr', 1e-6)}")
    else:
        scheduler = None

    train_loss_best = np.inf
    val_loss_best = np.inf
    val_prAUC_best = 0
    train_metric_epoch = []
    val_metric_epoch = []
    start_epoch = 0
    
    # Resume from checkpoint in output dir
    if getattr(args, 'resume', False):
        ckpt_params = os.path.join(outputDir, "model_params.pth")
        ckpt_optimizer = os.path.join(outputDir, "model_optimizer.pth")
        ckpt_state = os.path.join(outputDir, "checkpoint_state.txt")
        if os.path.isfile(ckpt_params) and os.path.isfile(ckpt_state):
            try:
                state = torch.load(ckpt_params, map_location=args.device, weights_only=True)
            except TypeError:
                state = torch.load(ckpt_params, map_location=args.device)
            # DataParallel saves with 'module.' prefix; strip for single-module load
            if state and next(iter(state.keys()), "").startswith("module."):
                state = {k.replace("module.", "", 1): v for k, v in state.items()}
            model_to_load = model.module if getattr(model, 'module', None) else model
            model_to_load.load_state_dict(state, strict=True)
            if os.path.isfile(ckpt_optimizer):
                try:
                    opt_state = torch.load(ckpt_optimizer, map_location=args.device, weights_only=True)
                except TypeError:
                    opt_state = torch.load(ckpt_optimizer, map_location=args.device)
                optimizer.load_state_dict(opt_state)
            with open(ckpt_state, "r", encoding="utf-8") as f:
                line = f.read().strip().split("\t")
                start_epoch = int(line[0]) + 1
                if len(line) > 1:
                    val_prAUC_best = float(line[1])
            print(f"Resume: from epoch {start_epoch}, val_prAUC_best={val_prAUC_best}")
    else:
        print("--resume: model_params.pth or checkpoint_state.txt not found; training from scratch")

    # Training
    max_stopping_step = args.max_stopping_step
    stopping_step = 0

    print("Training epochs: ", args.epochs)
    train_start_time = time.perf_counter()
    for epoch in range(start_epoch, args.epochs):
        lr = optimizer.param_groups[0]['lr']
        # Training phase
        model.train()
        train_loss_meter = epochAverageMeter()
        train_p_meter = epochAverageMeter()
    
        with tqdm(total=(len(trainLoader))) as t:   
            t.set_description('epoch:{}/{}'.format(epoch, args.epochs))
            for (inputs, targets, targets_ctl) in trainLoader:
                inputs = inputs.to(args.device)
                targets = targets.to(args.device)
                targets_ctl = targets_ctl.to(args.device)
    
                if EncodeMode == "OneHot" or args.input_channels == 5:
                    outputs = model(inputs.float())
    
                if (args.modelFramework == "Control"): 
                    w = outputs[2]
                    w = w.unsqueeze(axis=1)
                    target_predict = outputs[0]
                    ctl_predict = outputs[1]
    
                    if(args.log):
                        target_predict = torch.exp(target_predict) - 1
                        ctl_predict = torch.exp(ctl_predict) - 1
                        targets = torch.exp(targets) - 1
                
                    treat_predict = (target_predict * w) + (1-w) * ctl_predict
                    loss_treat = getattr(loss_fun(), args.loss_fun_name)(treat_predict, targets)
                    loss_ctl = getattr(loss_fun(), args.loss_fun_name)(ctl_predict, targets_ctl)
    
                    loss = loss_treat + loss_ctl
                    p = pearsonr_tensor(treat_predict, targets, batch_first=True).nanmean()
    
                elif(args.modelFramework == "Treat"):
                    if args.loss_fun_name == "CombinedProfileLoss":
                        logits_profile = outputs['logits_profile']    # (B, L)
                        log_total_count = outputs['log_total_count']  # (B, 1)
    
                        if args.log:
                            targets = torch.exp(targets) - 1
    
                        if _use_combined_loss:
                            loss_dict = criterion(logits_profile, log_total_count, targets)
                            loss = loss_dict['loss']
                        else:
                            prob = torch.softmax(logits_profile, dim=-1)
                            reconstructed = prob * torch.exp(log_total_count)
                            loss = getattr(loss_fun(), _effective_baseline_loss)(reconstructed, targets)
    
                        with torch.no_grad():
                            prob = torch.softmax(logits_profile, dim=-1)
                            reconstructed = prob * torch.exp(log_total_count)
                            p = pearsonr_tensor(reconstructed, targets, batch_first=True).nanmean()
                    else:
                        logits_profile = outputs['logits_profile']
                        outputs_compat = F.softplus(logits_profile)
                        if args.log:
                            outputs_compat = torch.exp(outputs_compat) - 1
                            targets = torch.exp(targets) - 1
                        loss = getattr(loss_fun(), args.loss_fun_name)(outputs_compat, targets)
                        p = pearsonr_tensor(outputs_compat, targets, batch_first=True).nanmean()
    
                train_loss_meter.update(loss.item())
    
                if (not torch.isnan(p)): 
                    train_p_meter.update(p.item())
    
                optimizer.zero_grad()
                loss.backward()
                if getattr(args, 'grad_clip', 0) > 0:
                    torch.nn.utils.clip_grad_norm_(model.parameters(), args.grad_clip)
                optimizer.step()
    
                if (args.modelFramework == "Treat"):
                    if args.loss_fun_name == "CombinedProfileLoss" and _use_combined_loss:
                        t.set_postfix(
                            loss='{:.3f}'.format(train_loss_meter.avg),
                            shape='{:.3f}'.format(loss_dict['loss_shape'].item()),
                            jitter='{:.3f}'.format(loss_dict['loss_jitter'].item()),
                            count='{:.3f}'.format(loss_dict['loss_count'].item()),
                            p='{:.3f}'.format(train_p_meter.avg),
                            lr='{:.5f}'.format(lr))
                    else:
                        t.set_postfix(loss='{:.3f}'.format(train_loss_meter.avg), p='{:.3f}'.format(train_p_meter.avg), lr='{:.5f}'.format(lr))
                elif(args.modelFramework == "Control"):  
                    t.set_postfix(loss_treat='{:.3f}'.format(loss_treat), loss_ctl='{:.3f}'.format(loss_ctl), weight='{:.3f}'.format(torch.mean(w)), loss='{:.3f}'.format(train_loss_meter.avg), p='{:.3f}'.format(train_p_meter.avg))
                t.update(1)
    
        val_mse_meter = epochAverageMeter()
        val_loss_meter = epochAverageMeter()
        val_p_meter = epochAverageMeter()
        val_sp_meter = epochAverageMeter()
        val_AUC_meter = epochAverageMeter()
        val_prAUC_meter = epochAverageMeter()
    
        # validation data
        model.eval()
        for (inputs, targets, targets_ctl) in valLoader:
            inputs = inputs.to(args.device)
            targets = targets.to(args.device)
            targets_ctl = targets_ctl.to(args.device)
    
            if EncodeMode == "OneHot" or args.input_channels == 5:
                outputs = model(inputs.float())
        
            if (args.modelFramework == "Control"): 
                w = outputs[2]
                w = w.unsqueeze(axis=1)
                target_predict = outputs[0]
                ctl_predict = outputs[1]
    
                if(args.log):
                    target_predict = torch.exp(target_predict) - 1
                    ctl_predict = torch.exp(ctl_predict) - 1
                    targets = torch.exp(targets) - 1
            
                treat_predict = (target_predict * w) + (1-w) * ctl_predict
                loss_treat = getattr(loss_fun(), args.loss_fun_name)(treat_predict, targets)
                loss_ctl = getattr(loss_fun(), args.loss_fun_name)(ctl_predict, targets_ctl)
    
                valLoss = loss_treat + loss_ctl
                p = pearsonr_tensor(treat_predict, targets, batch_first=True).nanmean()
    
                outputs_np = treat_predict.cpu().detach().numpy()
                targets_np = targets.cpu().detach().numpy()
    
            elif(args.modelFramework == "Treat"):
                logits_profile = outputs['logits_profile']
                log_total_count = outputs['log_total_count']
    
                if args.log:
                    targets = torch.exp(targets) - 1
    
                if args.loss_fun_name == "CombinedProfileLoss":
                    if _use_combined_loss:
                        loss_dict = criterion(logits_profile, log_total_count, targets)
                        valLoss = loss_dict['loss']
                    else:
                        prob = torch.softmax(logits_profile, dim=-1)
                        reconstructed = prob * torch.exp(log_total_count)
                        valLoss = getattr(loss_fun(), _effective_baseline_loss)(reconstructed, targets)
                else:
                    outputs_compat = F.softplus(logits_profile)
                    valLoss = getattr(loss_fun(), args.loss_fun_name)(outputs_compat, targets)
    
                # Reconstruct profile for metrics
                prob = torch.softmax(logits_profile, dim=-1)
                reconstructed = prob * torch.exp(log_total_count)
                p = pearsonr_tensor(reconstructed, targets, batch_first=True).nanmean()
    
                outputs_np = reconstructed.cpu().detach().numpy()
                targets_np = targets.cpu().detach().numpy()
        
            precise = False
            if precise:
                p_cor = pearsonr_tensor(outputs, targets).nanmean()
                roc_auc_values, pr_auc_values, sp_cor_values = calculate_batch_metrics_seq(outputs_np, targets_np)
    
                # Mse
                mse_values = np.nanmean((outputs_np - targets_np)**2, axis=None)
            
                sp_cor = np.nanmean(sp_cor_values)
                roc_auc = np.nanmean(roc_auc_values)
                pr_auc = np.nanmean(pr_auc_values)
    
            else:
                _sig_auc, _sig_ap, peak_auc, peak_ap, sp_cor, p_cor = calculate_batch_metrics_whole(outputs_np, targets_np)
                roc_auc, pr_auc = peak_auc, peak_ap
    
                # Mse
                mse_values = np.nanmean((outputs_np - targets_np)**2, axis = None)
    
            val_loss_meter.update(valLoss.item())
        
            if (not np.isnan(p_cor)):
                val_p_meter.update(p_cor.item())
            val_sp_meter.update(sp_cor)
            val_AUC_meter.update(roc_auc)
            if (not np.isnan(pr_auc)):
                val_prAUC_meter.update(pr_auc)
            val_mse_meter.update(mse_values)
    
        # lr decay (skip when using CosineAnnealingLR)
        if scheduler is None:
            if val_loss_best > val_loss_meter.avg:
                decay_count = 0
                val_loss_best = val_loss_meter.avg
            else:
                decay_count += 1
                if decay_count >= args.lr_decay_count_max: 
                    lr = optimizer.param_groups[0]['lr'] / 2  #
                    for param_group in optimizer.param_groups:
                        param_group['lr'] = lr
                    if lr < args.lr_min:    # stop training
                        print(f"Early stoping by lr, epoch={epoch+1}, ")
                        break
                    decay_count = 0
        else:
            val_loss_best = min(val_loss_best, val_loss_meter.avg)
    
        # Early stop
        if (val_prAUC_best < val_prAUC_meter.avg): 
            save_model(outputDir, model, optimizer, name="model")
            with open(os.path.join(outputDir, "checkpoint_state.txt"), "w", encoding="utf-8") as f:
                f.write(f"{epoch}\t{val_prAUC_meter.avg}\n")
            stopping_step = 0
            val_prAUC_best = val_prAUC_meter.avg
        else:
            stopping_step += 1
            if stopping_step >= max_stopping_step:
                print(f"Early stoping by val_loss, epoch={epoch+1}")
                break
    
        # writer
        writer.add_scalar('train_loss', train_loss_meter.avg, epoch)
        writer.add_scalar('train_p', train_p_meter.avg, epoch)
    
        writer.add_scalar('val_loss', val_loss_meter.avg, epoch)
        writer.add_scalar('val_p', val_p_meter.avg, epoch)
        writer.add_scalar('val_sp', val_sp_meter.avg, epoch)
        writer.add_scalar('val_AUC', val_AUC_meter.avg, epoch)
        writer.add_scalar('val_prAUC', val_prAUC_meter.avg, epoch)
        writer.add_scalar('val_mse', val_mse_meter.avg, epoch)
        writer.add_scalar('lr', lr, epoch)
    
        val_metric_epoch.append([val_loss_meter.avg, val_p_meter.avg, val_sp_meter.avg, val_AUC_meter.avg, val_prAUC_meter.avg, val_mse_meter.avg])
        train_metric_epoch.append([train_loss_meter.avg, train_p_meter.avg])
    
        if scheduler is not None:
            if getattr(args, "use_reduce_on_plateau", False):
                scheduler.step(val_prAUC_meter.avg)
            else:
                scheduler.step()
    
        train_end_time = time.perf_counter()
        train_duration_sec = train_end_time - train_start_time
    
        val_metric_epoch_df = pd.DataFrame(val_metric_epoch)
        train_metric_epoch_df = pd.DataFrame(train_metric_epoch)
    
        val_metric_epoch_df.columns = ["loss", "pearson", "spearman", "ROC_AUC", "PR_AUC", "mse"]
        train_metric_epoch_df.columns = ["loss", "pearson"]
        val_metric_epoch_df.to_csv(f"{outputDir}/val_metric.txt", sep="\t", index =False)
        train_metric_epoch_df.to_csv(f"{outputDir}/train_metric.txt", sep="\t", index =False)
    
    if args.log:
        print(f"{protein}, {cell}, {args.loss_fun_name}, log, loss={val_loss_meter.avg}, pearson={val_p_meter.avg}, spearman={val_sp_meter.avg}, auc={val_AUC_meter.avg}, prauc={val_prAUC_meter.avg}, mse={val_mse_meter.avg}")
    else:
        print(f"{protein}, {cell}, {args.loss_fun_name}, non-log, loss={val_loss_meter.avg}, pearson={val_p_meter.avg}, spearman={val_sp_meter.avg}, auc={val_AUC_meter.avg}, prauc={val_prAUC_meter.avg}, mse={val_mse_meter.avg}")

    print(f"\n=== Training time: {train_duration_sec:.2f}s ({train_duration_sec/60:.2f} min) ===")

    # Ablation pipelines check outputDir/model_params.pth for success; if val never beat val_prAUC_best,
    # save_model may never run, yet test still uses the in-memory model for test_metrics_summary.txt.
    # Persist the final-epoch weights here for testing and downstream copy.
    _mp = f"{outputDir}/model_params.pth"
    if not os.path.isfile(_mp):
        save_model(outputDir, model, optimizer, name="model")
        print(f"[checkpoint] Wrote {_mp} (val-based save never triggered; using last epoch for test/downstream)")

    # Test set evaluation (if --testFile; supports comma-separated paths)
    test_duration_sec = None
    if args.testFile and _test_files_exist(args.testFile):
        print(f"\n=== Evaluating on test set: {args.testFile} ===")
        test_start_time = time.perf_counter()
        if args.input_channels == 5:
            test_data, test_label, test_label_ctl = read_data_5ch(args.testFile, args.signalMax, args.filter)
        else:
            test_data, test_label, test_label_ctl = read_data(args.testFile, args.signalMax, args.filter, EncodeMode)
        if getattr(args, "ablation_no_conservative", False):
            test_data = test_data[..., :4]
        if args.log:
            test_label = np.log(np.asarray(test_label) + 1)
            test_label_ctl = np.log(np.asarray(test_label_ctl) + 1)
        else:
            test_label = np.asarray(test_label)
            test_label_ctl = np.asarray(test_label_ctl)
        test_data = torch.from_numpy(np.asarray(test_data))
        test_label = torch.from_numpy(test_label).float()
        test_label_ctl = torch.from_numpy(test_label_ctl).float()
        test_dataset = torch.utils.data.TensorDataset(test_data, test_label, test_label_ctl)
        testLoader = torch.utils.data.DataLoader(test_dataset, batch_size=args.bsz, shuffle=False, drop_last=False)
        del test_data, test_label, test_label_ctl
        gc.collect()

        # Load best model for test evaluation
        model_path = f"{outputDir}/model_params.pth"
        if os.path.isfile(model_path):
            try:
                state = torch.load(model_path, map_location=args.device, weights_only=True)
            except TypeError:
                state = torch.load(model_path, map_location=args.device)
            model.load_state_dict(state)
        model.eval()

        all_preds = []
        all_targets = []
        with torch.no_grad():
            for (inputs, targets, targets_ctl) in testLoader:
                inputs = inputs.to(args.device)
                targets = targets.to(args.device)
                targets_ctl = targets_ctl.to(args.device)
                if EncodeMode == "OneHot" or args.input_channels == 5:
                    outputs = model(inputs.float())
                if args.modelFramework == "Control":
                    w = outputs[2].unsqueeze(1)
                    target_predict = outputs[0]
                    ctl_predict = outputs[1]
                    if args.log:
                        target_predict = torch.exp(target_predict) - 1
                        ctl_predict = torch.exp(ctl_predict) - 1
                        targets = torch.exp(targets) - 1
                    treat_predict = (target_predict * w) + (1 - w) * ctl_predict
                else:
                    if args.log:
                        targets = torch.exp(targets) - 1
                    logits_profile = outputs['logits_profile']
                    log_total_count = outputs['log_total_count']
                    prob = torch.softmax(logits_profile, dim=-1)
                    treat_predict = prob * torch.exp(log_total_count)
                all_preds.append(treat_predict.cpu().numpy())
                all_targets.append(targets.cpu().numpy())
        test_end_time = time.perf_counter()
        test_duration_sec = test_end_time - test_start_time

        outputs_np = np.concatenate(all_preds, axis=0)
        targets_np = np.concatenate(all_targets, axis=0)
        signal_auc, signal_ap, peak_auc, peak_ap, spearman, mean_pcc = calculate_batch_metrics_whole(outputs_np, targets_np)
        mean_pcc_val = "nan" if np.isnan(mean_pcc) else float(mean_pcc)

        print("\n" + "=" * 60)
        print("Test Set Metrics")
        print("=" * 60)
        print(f"  Mean PCC (Pearson):     {mean_pcc_val if mean_pcc_val == 'nan' else f'{float(mean_pcc_val):.6f}'}")
        print(f"  Spearman correlation:  {spearman:.6f}")
        print(f"  Signal AUC:             {signal_auc:.6f}")
        print(f"  Signal AP:              {signal_ap:.6f}")
        print(f"  Peak AUC:               {peak_auc:.6f}")
        print(f"  Peak AP:                {peak_ap:.6f}")
        print("=" * 60)
        print(f"Test time: {test_duration_sec:.2f}s ({test_duration_sec/60:.2f} min)")
        print("=" * 60)
        summary_path = f"{outputDir}/test_metrics_summary.txt"
        with open(summary_path, "w", encoding="utf-8") as f:
            f.write("Metric\tValue\n")
            f.write(f"Mean_PCC\t{mean_pcc_val}\n")
            f.write(f"Spearman\t{spearman}\n")
            f.write(f"Signal_AUC\t{signal_auc}\n")
            f.write(f"Signal_AP\t{signal_ap}\n")
            f.write(f"Peak_AUC\t{peak_auc}\n")
            f.write(f"Peak_AP\t{peak_ap}\n")
            f.write(f"Training_time_sec\t{train_duration_sec}\n")
            f.write(f"Testing_time_sec\t{test_duration_sec}\n")
        print(f"\nMetrics and timing saved to: {summary_path}")
        test_metric_path = f"{outputDir}/test_metric.txt"
        with open(test_metric_path, "w", encoding="utf-8") as f:
            f.write("Mean_PCC\tSpearman\tSignal_AUC\tSignal_AP\tPeak_AUC\tPeak_AP\tTraining_time_sec\tTesting_time_sec\n")
            f.write(f"{mean_pcc_val}\t{spearman}\t{signal_auc}\t{signal_ap}\t{peak_auc}\t{peak_ap}\t{train_duration_sec}\t{test_duration_sec}\n")
        with open(f"{outputDir}/run_info.txt", "a", encoding="utf-8") as f:
            f.write(f"testing_time_sec\t{test_duration_sec}\n")
            f.write(f"testing_time_min\t{test_duration_sec/60}\n")
    elif args.testFile:
        print(f"\nTest file not found: {args.testFile}; skipping test evaluation.")


