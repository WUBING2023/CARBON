#!/usr/bin/env bash
# Compatibility entry: included AARS CPU saved-metric check, not MARCP equivalence.
set -euo pipefail
CARBON_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$CARBON_ROOT"
python3 scripts/verify_demo_cpu.py "$@"
