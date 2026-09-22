# ==============================================================================
# BLOCK 5: KINETIC RATE CONSTANTS (BEP + EYRING + DETAILED BALANCE)
# ==============================================================================
import inspect
import numpy as np
from calculate_reaction_thermo import calculate_reaction_thermo
from calculate_gas_thermo import R_SI, KB_SI, H_SI, EV_TO_KJ_MOL

def _get_bep_parameters(bep_params=None):
    """Retrieves BEP parameters from arguments, globals, or calling frame."""
    if bep_params is not None:
        return bep_params
    if 'family_bep_parameters' in globals():
        return globals()['family_bep_parameters']
    frame = inspect.currentframe().f_back
    while frame:
        if 'family_bep_parameters' in frame.f_globals:
            return frame.f_globals['family_bep_parameters']
        if 'family_bep_parameters' in frame.f_locals:
            return frame.f_locals['family_bep_parameters']
        frame = frame.f_back
    raise NameError("family_bep_parameters is not defined. Please define it in your notebook or pass bep_params explicitly.")

def calculate_rate_constants(
    rxn_id: str,
    T_K: float,
    bep_params: dict = None,
    reactions_net: dict = None,
    species_db: dict = None,
    mode: str = 'qRRHO',
    **kwargs
) -> dict:
    """
    Computes k_f and k_r at temperature T_K strictly enforcing detailed balance.
    
    ΔG‡_f = max(E0, E0 + α * ΔG_rxn)
    k_f = (kB * T / h) * exp(-ΔG‡_f / RT)
    k_r = k_f / K_eq
    """
    thermo = calculate_reaction_thermo(
        rxn_id, T_K,
        reactions_net=reactions_net,
        species_db=species_db,
        mode=mode,
        **kwargs
    )
    rxn_class = thermo['class']
    all_bep = _get_bep_parameters(bep_params)
    bep = all_bep.get(rxn_class, all_bep.get('default', {'E0_eV': 0.80, 'alpha': 0.50}))
    
    E0_kJ_mol = bep['E0_eV'] * EV_TO_KJ_MOL
    alpha = bep['alpha']
    dG_rxn_kJ_mol = thermo['dG_rxn_kJ_mol']
    
    # 1. Forward Barrier: ΔG‡_f = max(E0, E0 + α * ΔG_rxn)
    dG_barrier_f_kJ_mol = max(E0_kJ_mol, E0_kJ_mol + alpha * dG_rxn_kJ_mol)
    
    # 2. Eyring Forward Rate: k_f = (kB*T / h) * exp(-ΔG‡ / RT)
    eyring_prefactor = (KB_SI * T_K) / H_SI
    k_f = eyring_prefactor * np.exp(-(dG_barrier_f_kJ_mol * 1000.0) / (R_SI * T_K))
    
    # 3. Microscopic Reversibility: k_r = k_f / K_eq
    K_eq = thermo['K_eq']
    k_r = k_f / K_eq if K_eq > 1e-300 else 0.0
    
    # Reverse barrier from detailed balance
    dG_barrier_r_kJ_mol = dG_barrier_f_kJ_mol - dG_rxn_kJ_mol
    
    return {
        'rxn_id': rxn_id,
        'T_K': T_K,
        'dG_barrier_f_kJ_mol': dG_barrier_f_kJ_mol,
        'dG_barrier_f_eV': dG_barrier_f_kJ_mol / EV_TO_KJ_MOL,
        'dG_barrier_r_kJ_mol': dG_barrier_r_kJ_mol,
        'dG_barrier_r_eV': dG_barrier_r_kJ_mol / EV_TO_KJ_MOL,
        'k_f': k_f,
        'k_r': k_r,
        'K_eq': K_eq
    }