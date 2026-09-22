# ==============================================================================
# BLOCK 6: ARRHENIUS REGRESSION & PARAMETER EXTRACTION
# ==============================================================================
import numpy as np
import pandas as pd
from calculate_rate_constants import calculate_rate_constants
from calculate_gas_thermo import R_SI, EV_TO_KJ_MOL
from calculate_reaction_thermo import _get_reactions_network

def fit_modified_arrhenius(
    rxn_id: str,
    T_grid: np.ndarray,
    direction: str = 'f',
    bep_params: dict = None,
    reactions_net: dict = None,
    species_db: dict = None,
    mode: str = 'qRRHO',
    **kwargs
) -> dict:
    """
    Fits k(T) to the Modified Arrhenius equation:
        k(T) = A * (T / T0)^beta * exp(-Ea / RT)
    via Ordinary Least Squares in log-space:
        ln k(T) = ln A + beta * ln(T / T0) - (Ea / R) * (1 / T)
    """
    T0 = 298.15
    k_values = []
    
    for T in T_grid:
        rates = calculate_rate_constants(
            rxn_id, T,
            bep_params=bep_params,
            reactions_net=reactions_net,
            species_db=species_db,
            mode=mode,
            **kwargs
        )
        k_val = rates['k_f'] if direction == 'f' else rates['k_r']
        k_values.append(k_val)
        
    k_values = np.array(k_values)
    
    # Linear system: ln(k) = ln(A) + beta * ln(T/T0) - (Ea/R) * (1/T)
    Y = np.log(k_values)
    X = np.column_stack([
        np.ones_like(T_grid),
        np.log(T_grid / T0),
        -1.0 / (R_SI * T_grid)
    ])
    
    params, residuals, rank, s = np.linalg.lstsq(X, Y, rcond=None)
    ln_A, beta, Ea_J_mol = params
    
    # Compute R^2 determination coefficient
    Y_pred = X @ params
    ss_tot = np.sum((Y - np.mean(Y))**2)
    ss_res = np.sum((Y - Y_pred)**2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0
    
    return {
        'rxn_id': rxn_id,
        'direction': direction,
        'A': np.exp(ln_A),
        'beta': beta,
        'Ea_kJ_mol': Ea_J_mol / 1000.0,
        'Ea_eV': (Ea_J_mol / 1000.0) / EV_TO_KJ_MOL,
        'R2': r2
    }

def generate_arrhenius_summary(
    T_grid: np.ndarray,
    reactions_net: dict = None,
    bep_params: dict = None,
    species_db: dict = None,
    mode: str = 'qRRHO',
    **kwargs
) -> list:
    """
    Generates summary list of forward and reverse Arrhenius parameters across all reactions in the network.
    """
    net = _get_reactions_network(reactions_net)
    summary = []
    for rxn_id in net.keys():
        fit_f = fit_modified_arrhenius(
            rxn_id, T_grid, direction='f',
            bep_params=bep_params, reactions_net=net,
            species_db=species_db, mode=mode,
            **kwargs
        )
        fit_r = fit_modified_arrhenius(
            rxn_id, T_grid, direction='r',
            bep_params=bep_params, reactions_net=net,
            species_db=species_db, mode=mode,
            **kwargs
        )
        summary.append({
            'Reaction': rxn_id,
            'Class': net[rxn_id]['class'],
            'A_f (s^-1 or M^-n s^-1)': f"{fit_f['A']:.3e}",
            'beta_f': f"{fit_f['beta']:.2f}",
            'Ea_f (kJ/mol)': f"{fit_f['Ea_kJ_mol']:.2f}",
            'Ea_f (eV)': f"{fit_f['Ea_eV']:.3f}",
            'A_r': f"{fit_r['A']:.3e}",
            'beta_r': f"{fit_r['beta']:.2f}",
            'Ea_r (kJ/mol)': f"{fit_r['Ea_kJ_mol']:.2f}",
            'Ea_r (eV)': f"{fit_r['Ea_eV']:.3f}",
            'R2': f"{fit_f['R2']:.4f}"
        })
    return summary
