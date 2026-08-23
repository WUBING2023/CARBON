#!/usr/bin/env bash
# Train CARBON on minimal demo data (chr1/2/3) from configs/demo.yaml; output to carbon_trained_model/AARS_K562_ENCSR825SVO_demo/
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="${1:-carbon_trained_model/AARS_K562_ENCSR825SVO_demo}"
EPOCHS="${DEMO_EPOCHS:-3}"
TASK="AARS_K562_ENCSR825SVO"
TRAIN="data/h5data/${TASK}/_5ch/${TASK}_chr1_5ch.h5"
VAL="data/h5data/${TASK}/_5ch/${TASK}_chr2_5ch.h5"
TEST="data/h5data/${TASK}/_5ch/${TASK}_chr3_5ch.h5"

for f in "$TRAIN" "$VAL" "$TEST"; do
  if [[ ! -f "$f" ]]; then
    echo "Missing data file: $f" >&2
    exit 1
  fi
done

python3 -u carbon_train.py \
  --lr 0.000785 \
  --CUDA 0 \
  --signalMax 0 \
  --seqLength 101 \
  --bsz 64 \
  --encode OneHot \
  --modelFramework Treat \
  --epochs "$EPOCHS" \
  --input_channels 5 \
  -train "$TRAIN" \
  -validation "$VAL" \
  -test "$TEST" \
  --task "$TASK" \
  --lossM CombinedProfileLoss \
  --alpha 1.42 --beta 0.54 --gamma 1.80 \
  --weight_decay 0.000114 \
  --grad_clip 2.0 \
  --output "$OUT"

echo "Done. Weights: $OUT/model_params.pth"
