"""Gradient-based attribution and sequence-level interpretation."""

from .marcp_attribution import (
    AttributionMethod,
    GradientBasedMethod,
    Saliency,
    GradientXInput,
    IntegratedGradients,
    MARCP_integrated_gradients_fast,
    MARCP_integrated_gradients_fast_interpret_target,
    MARCP_integrated_gradients_fast_interpret,
    attributionFromSeq,
    interpret1HeadFromSeq,
)

__all__ = [
    "AttributionMethod",
    "GradientBasedMethod",
    "Saliency",
    "GradientXInput",
    "IntegratedGradients",
    "MARCP_integrated_gradients_fast",
    "MARCP_integrated_gradients_fast_interpret_target",
    "MARCP_integrated_gradients_fast_interpret",
    "attributionFromSeq",
    "interpret1HeadFromSeq",
]
