#!/usr/bin/env bash
# Evaluation only: load model_params.pth and run metrics on test h5 (carbon_train.py --eval_only)
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

TASK="AARS_K562_ENCSR825SVO"
TEST="data/h5data/${TASK}/_5ch/${TASK}_chr3_5ch.h5"
MODEL="${1:-carbon_trained_model/AARS_K562_ENCSR825SVO_demo/model_params.pth}"
OUT="${2:-carbon_trained_model/AARS_K562_ENCSR825SVO_demo_eval_only}"

if [[ ! -f "$MODEL" ]]; then
  echo "Not found: $MODEL" >&2
  exit 1
fi

python3 -u carbon_train.py \
  --eval_only \
  --load_model "$MODEL" \
  --CUDA 0 \
  --signalMax 0 \
  --seqLength 101 \
  --bsz 64 \
  --encode OneHot \
  --modelFramework Treat \
  --input_channels 5 \
  -test "$TEST" \
  --task "$TASK" \
  --lossM CombinedProfileLoss \
  --alpha 1.42 --beta 0.54 --gamma 1.80 \
  --output "$OUT"
