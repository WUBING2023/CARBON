"""CARBON model variants: main (multi-scale motif + BiGRU) and RNAFM-based."""

from .carbon_model import CARBON, MultiScaleMotifExtractor, ModernBiGRU, ModernBiLSTM, Conv1dProfileHead

__all__ = ["CARBON", "MultiScaleMotifExtractor", "ModernBiGRU", "ModernBiLSTM", "Conv1dProfileHead"]
