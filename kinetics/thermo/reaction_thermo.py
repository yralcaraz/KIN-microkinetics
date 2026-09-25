# ==============================================================================

# BLOCK 4: REACTION THERMODYNAMICS & EQUILIBRIUM CONSTANT

# This module calculates the reaction thermodynamics (ΔH_rxn, ΔS_rxn, ΔG_rxn) 
# and equilibrium constant (K_eq) for a given reaction at a specified temperature.

# The main function `calculate_reaction_thermo` takes a reaction ID, temperature,
# and optional reaction network and species database, and returns a dictionary
# containing the computed thermodynamic properties and equilibrium constant. 

# ==============================================================================
import inspect
import numpy as np

try:
    from kinetics.thermo.solution_gibbs import calculate_solution_gibbs
    from kinetics.thermo.gas_thermo import R_SI, EV_TO_KJ_MOL
except ImportError:
    from calculate_solution_gibbs import calculate_solution_gibbs
    from calculate_gas_thermo import R_SI, EV_TO_KJ_MOL


DEFAULT_REACTIONS_NETWORK = {
    'R1': {'reactants': {'TMSPA': 1, 'H2O': 1},     'products': {'BMSPA': 1, 'TMSOH': 1},   'class': 'hydrolysis'},
    'R2': {'reactants': {'BMSPA': 1, 'H2O': 1},     'products': {'MMSPA': 1, 'TMSOH': 1},   'class': 'hydrolysis'},
    'R3': {'reactants': {'MMSPA': 1, 'H2O': 1},     'products': {'H3PO4': 1, 'TMSOH': 1},   'class': 'hydrolysis'},
    'R4': {'reactants': {'TMSOH': 2},               'products': {'HMDSO': 1, 'H2O': 1},     'class': 'condensation'},
    'R5': {'reactants': {'TMSPA': 1, 'TMSOH': 1},   'products': {'BMSPA': 1, 'HMDSO': 1},   'class': 'transfer'},
    'R6': {'reactants': {'BMSPA': 1, 'TMSOH': 1},   'products': {'MMSPA': 1, 'HMDSO': 1},   'class': 'transfer'},
    'R7': {'reactants': {'MMSPA': 1, 'TMSOH': 1},   'products': {'H3PO4': 1, 'HMDSO': 1},   'class': 'transfer'},
    'R8': {'reactants': {'EC': 1, 'TMSOH': 1},      'products': {'TMSOEG': 1, 'CO2': 1},     'class': 'solvent_attack'},
    'R9': {'reactants': {'EC': 1, 'TMSOEG': 1},     'products': {'TMSOdiEG': 1, 'CO2': 1},   'class': 'solvent_attack'},
}

def _get_reactions_network(reactions_net=None):
    """Retrieves reactions network from arguments, globals, calling frame, or default network."""
    if reactions_net is not None:
        return reactions_net
    if 'reactions_network' in globals():
        return globals()['reactions_network']
    frame = inspect.currentframe().f_back
    while frame:
        if 'reactions_network' in frame.f_globals:
            return frame.f_globals['reactions_network']
        if 'reactions_network' in frame.f_locals:
            return frame.f_locals['reactions_network']
        frame = frame.f_back
    return DEFAULT_REACTIONS_NETWORK.copy()

def calculate_reaction_thermo(rxn_id: str, T_K: float, reactions_net: dict = None, species_db: dict = None, mode: str = 'qRRHO', **kwargs) -> dict:
    """
    Computes ΔH_rxn, ΔS_rxn, ΔG_rxn, and K_eq at temperature T_K for reaction rxn_id.
    """
    net = _get_reactions_network(reactions_net)
    rxn = net[rxn_id]
    dG_rxn_kJ_mol = 0.0
    dH_rxn_kJ_mol = 0.0
    dS_rxn_J_mol_K = 0.0
    
    # Products (+)
    for p, coeff in rxn['products'].items():
        th = calculate_solution_gibbs(p, T_K, species_db=species_db, mode=mode, **kwargs)
        dG_rxn_kJ_mol += coeff * th['G_sol_kJ_mol']
        dH_rxn_kJ_mol += coeff * th['H_sol_kJ_mol']
        dS_rxn_J_mol_K += coeff * th['S_sol_J_mol_K']
        
    # Reactants (-)
    for r, coeff in rxn['reactants'].items():
        th = calculate_solution_gibbs(r, T_K, species_db=species_db, mode=mode, **kwargs)
        dG_rxn_kJ_mol -= coeff * th['G_sol_kJ_mol']
        dH_rxn_kJ_mol -= coeff * th['H_sol_kJ_mol']
        dS_rxn_J_mol_K -= coeff * th['S_sol_J_mol_K']
        
    # Equilibrium constant K_eq = exp(-ΔG / RT)
    dG_J_mol = dG_rxn_kJ_mol * 1000.0
    K_eq = np.exp(-dG_J_mol / (R_SI * T_K))
    
    return {
        'rxn_id': rxn_id,
        'class': rxn['class'],
        'dH_rxn_kJ_mol': dH_rxn_kJ_mol,
        'dS_rxn_J_mol_K': dS_rxn_J_mol_K,
        'dG_rxn_kJ_mol': dG_rxn_kJ_mol,
        'dG_rxn_eV': dG_rxn_kJ_mol / EV_TO_KJ_MOL,
        'K_eq': K_eq
    }