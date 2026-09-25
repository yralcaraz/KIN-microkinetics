"""Microkinetics engine: Transition State Theory, BEP, Marcus, and Arrhenius regression."""

try:
    from kinetics.microkinetics.rate_constants import calculate_rate_constants, bep_eyring, marcus_eyring
    from kinetics.microkinetics.arrhenius import fit_modified_arrhenius, generate_arrhenius_summary
except ImportError:
    from calculate_rate_constants import calculate_rate_constants, bep_eyring, marcus_eyring
    from fit_modified_arrhenius import fit_modified_arrhenius, generate_arrhenius_summary

__all__ = [
    "calculate_rate_constants",
    "bep_eyring",
    "marcus_eyring",
    "fit_modified_arrhenius",
    "generate_arrhenius_summary",
]
