# ==============================================================================
# BLOCK 8: VIRTUAL OPERANDO SPECTROMETER (SYNTHETIC 29Si / 31P NMR ENGINE)
# ==============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Standard 29Si NMR data from Peter's DFT GIAO calculations (tank_model.ipynb)
# Chemical shifts (ppm) relative to TMS (delta = 0 ppm), and Si multiplicity
DEFAULT_NMR_29SI = {
    'TMSPA':    (-19.86, 3),   # Three equivalent trimethylsilyl ester groups
    'BMSPA':    (-18.01, 2),   # Two equivalent trimethylsilyl ester groups
    'MMSPA':    (-17.58, 1),   # One trimethylsilyl group
    'TMSOH':    (-15.57, 1),   # Free silanol
    'siloxyl':  (-6.92, 2),    # Hexamethyldisiloxane (bridged dimer)
    'TMSOEG':   (-18.01, 1),   # Silylated ethylene glycol mono-adduct
    'TMSOdiEG': (-19.11, 1),   # Silylated diethylene glycol ether
    'TMSOCH3':  (-17.58, 1),   # Silylated methanol
}

def lorentzian(x: np.ndarray, x0: float, fwhm: float = 0.8) -> np.ndarray:
    """
    Normalized Lorentzian line shape.
    L(x) = (gamma^2) / ((x - x0)^2 + gamma^2), where gamma = fwhm / 2.
    """
    gamma = fwhm / 2.0
    return (gamma**2) / ((x - x0)**2 + gamma**2)

def simulate_virtual_nmr(
    t_s: np.ndarray,
    C_M: np.ndarray,
    species_list: list,
    nmr_data: dict = None,
    delta_range: tuple = (-35.0, 5.0),
    n_delta: int = 2000,
    lw: float = 0.8,
    n_snapshots: int = 50,
    **kwargs
) -> dict:
    """
    Convolutes time-resolved chemical concentrations into virtual NMR spectra.
    """
    if nmr_data is None:
        nmr_data = DEFAULT_NMR_29SI
        
    idx = {s: i for i, s in enumerate(species_list)}
    delta = np.linspace(delta_range[0], delta_range[1], n_delta)
    
    # Pick log-spaced snapshot indices
    nt = len(t_s)
    snap_idx = np.unique(np.round(np.logspace(0, np.log10(nt - 1), n_snapshots)).astype(int))
    snap_idx = snap_idx[snap_idx < nt]
    t_snap_h = t_s[snap_idx] / 3600.0
    
    # Build 2D intensity matrix (n_snapshots x n_delta) in mM·Si
    spectra_2d = np.zeros((len(snap_idx), len(delta)))
    
    for si, ti in enumerate(snap_idx):
        for sp, (shift, n_nuc) in nmr_data.items():
            if sp not in idx:
                continue
            conc_mM = C_M[idx[sp], ti] * 1000.0
            intensity = conc_mM * n_nuc
            spectra_2d[si] += intensity * lorentzian(delta, shift, lw)
            
    return {
        'delta_ppm': delta,
        'snap_idx': snap_idx,
        't_snap_h': t_snap_h,
        'spectra_2d': spectra_2d,
        'nmr_data': nmr_data,
        'lw': lw
    }
