#!/usr/bin/env bash
# Use weights from scripts/run_demo_train.sh; generate Figure-2-style four-track demo plot.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
mkdir -p figures/outputs

W="${1:-$ROOT/carbon_trained_model/AARS_K562_ENCSR825SVO_demo/model_params.pth}"
OUT="${2:-$ROOT/figures/outputs/demo_signal_track_4rows.png}"

if [[ ! -f "$W" ]]; then
  echo "Weights not found: $W — run: bash scripts/run_demo_train.sh" >&2
  exit 1
fi

if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
  python3 figures/plot_demo_signal_track.py --weights "$W" --out "$OUT"
else
  python3 figures/plot_demo_signal_track.py --weights "$W" --out "$OUT" --device cpu
fi

echo "Wrote: $OUT"
