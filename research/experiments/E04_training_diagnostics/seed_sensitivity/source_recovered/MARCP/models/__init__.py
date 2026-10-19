"""MARCP model variants: main (multi-scale motif + BiGRU) and RNAFM-based."""

from .marcp_model import MARCP, MultiScaleMotifExtractor, ModernBiGRU, ModernBiLSTM, Conv1dProfileHead

__all__ = ["MARCP", "MultiScaleMotifExtractor", "ModernBiGRU", "ModernBiLSTM", "Conv1dProfileHead"]
