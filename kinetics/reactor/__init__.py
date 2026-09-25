"""Reactor dynamics: Homogeneous batch and multi-stage temperature-programmed ODE engines."""

from kinetics.reactor.batch_reactor import simulate_tank_reactor, build_stoichiometric_matrix
from kinetics.reactor.protocol_reactor import (
    simulate_protocol_reactor,
    compute_recipe_molarities,
    build_default_protocol_schedule
)

__all__ = [
    "simulate_tank_reactor",
    "build_stoichiometric_matrix",
    "simulate_protocol_reactor",
    "compute_recipe_molarities",
    "build_default_protocol_schedule",
]
