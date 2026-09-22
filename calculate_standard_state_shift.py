# ==============================================================================
# BLOCK 3A: GAS-TO-SOLUTION STANDARD-STATE COMPRESSION CORRECTION
# ==============================================================================
import numpy as np

# Universal gas constant [J / (mol · K)]
R_SI = 8.314462618

def calculate_standard_state_shift(T_K: float, P_ref_Pa: float = 1.0e5, **kwargs) -> float:
    """
    Computes standard-state free energy shift ΔG°→*(T) in kJ/mol.
    
    ΔG°→*(T) = R * T * ln(C*_soln / C°_gas(T))
    where:
        C°_gas(T) = P° / (R * T)  in mol / m^3
        C*_soln   = 1000.0 mol / m^3 (1.0 M)
    """
    C_gas_standard = P_ref_Pa / (R_SI * T_K)  # mol / m^3
    C_soln_standard = 1000.0                   # 1 mol / L = 1000 mol / m^3
    
    # Free energy correction in J/mol:
    dG_shift_J_mol = R_SI * T_K * np.log(C_soln_standard / C_gas_standard)
    
    return dG_shift_J_mol / 1000.0  # kJ/mol