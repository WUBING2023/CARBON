"""Gradient-based attribution and sequence-level interpretation."""

from .carbon_attribution import (
    AttributionMethod,
    GradientBasedMethod,
    Saliency,
    GradientXInput,
    IntegratedGradients,
    CARBON_integrated_gradients_fast,
    CARBON_integrated_gradients_fast_interpret_target,
    CARBON_integrated_gradients_fast_interpret,
    attributionFromSeq,
    interpret1HeadFromSeq,
)

__all__ = [
    "AttributionMethod",
    "GradientBasedMethod",
    "Saliency",
    "GradientXInput",
    "IntegratedGradients",
    "CARBON_integrated_gradients_fast",
    "CARBON_integrated_gradients_fast_interpret_target",
    "CARBON_integrated_gradients_fast_interpret",
    "attributionFromSeq",
    "interpret1HeadFromSeq",
]
