"""Multi-nuclear Operando NMR Simulator: 29Si, 31P, 13C, 1H with Fast Proton Exchange.

Converts time-dependent species concentration matrices into realistic operando
multi-nuclear NMR spectra. Implements fast chemical exchange kinetics for labile
hydroxyl protons (-OH), automatic chemical shift window clustering (broken axes),
continuous peak tracking maps, and rigorous mass-balance water determination.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from kinetics.spectroscopy.molecular_symmetry import (
    DEFAULT_NUCLEI,
    build_referenced_nmr_sites
)

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
RULE = "#c3c2b7"
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def lorentzian(x: np.ndarray, x0: float | np.ndarray, fwhm: float) -> np.ndarray:
    """Normalized Lorentzian resonance lineshape function."""
    gamma = fwhm / 2.0
    return (gamma ** 2) / ((x - x0) ** 2 + gamma ** 2)


def site_peaks(
    element: str,
    C_mat: np.ndarray,
    species_idx: dict,
    nmr_sites_catalog: dict = None,
    solvent: tuple = ("EC",)
) -> list[tuple[str, np.ndarray, np.ndarray]]:
    """Computes peak centers and intensities (in mM of nuclei) for each resolved resonance."""
    catalog = nmr_sites_catalog if nmr_sites_catalog is not None else build_referenced_nmr_sites()
    element_sites = catalog.get(element, {})
    
    peaks = []
    labile_weight = np.zeros(C_mat.shape[1], dtype=float)
    labile_weighted_shift = np.zeros(C_mat.shape[1], dtype=float)
    
    for sp, sites in element_sites.items():
        if sp not in species_idx or sp in solvent:
            continue
            
        conc_mM = C_mat[species_idx[sp]] * 1000.0
        for shift, multiplicity, is_labile in sites:
            intensity = conc_mM * multiplicity
            if is_labile:
                labile_weight += intensity
                labile_weighted_shift += intensity * shift
            else:
                shift_vec = np.full(conc_mM.shape, shift, dtype=float)
                peaks.append((sp, shift_vec, intensity))
                
    if np.any(labile_weight > 0):
        denom = np.where(labile_weight > 0, labile_weight, 1.0)
        mean_labile_shift = labile_weighted_shift / denom
        peaks.append(("OH (exch.)", mean_labile_shift, labile_weight))
        
    return peaks


def auto_regions(shifts: list[float], pad: float = 4.0, max_gap: float = 40.0) -> list[tuple[float, float]]:
    """Clusters chemical shifts into disjoint windows [(high_ppm, low_ppm)]."""
    if not shifts:
        return [(10.0, -10.0)]
        
    s_sorted = sorted(shifts, reverse=True)
    regions = []
    current_high = s_sorted[0]
    
    for a, b in zip(s_sorted, s_sorted[1:]):
        if (a - b) > max_gap:
            regions.append((current_high + pad, a - pad))
            current_high = b
            
    regions.append((current_high + pad, s_sorted[-1] - pad))
    return regions


def simulate_multinuclear_spectra(
    element: str,
    C_mat: np.ndarray,
    species_idx: dict,
    x_grid: np.ndarray,
    lw: float = 0.4,
    nmr_sites_catalog: dict = None,
    solvent: tuple = ("EC",)
) -> np.ndarray:
    """Convolutes Lorentzian peaks across a frequency grid."""
    peaks = site_peaks(element, C_mat, species_idx, nmr_sites_catalog=nmr_sites_catalog, solvent=solvent)
    n_times = C_mat.shape[1]
    out_spectra = np.zeros((n_times, len(x_grid)), dtype=float)
    
    for _, shift_vec, inten_vec in peaks:
        out_spectra += inten_vec[:, None] * lorentzian(x_grid[None, :], shift_vec[:, None], lw)
        
    return out_spectra


def compute_water_mass_balance(
    C_mat: np.ndarray,
    species_idx: dict,
    initial_h2o_molarity: float,
    nmr_sites_catalog: dict = None
) -> dict:
    """Computes exact water concentration and OH carrier distribution over time."""
    catalog = nmr_sites_catalog if nmr_sites_catalog is not None else build_referenced_nmr_sites()
    h_sites = catalog.get("H", {})
    
    n_oh_per_species = {}
    for sp, sites in h_sites.items():
        if sp in species_idx:
            n_oh = sum(mult for _, mult, is_lab in sites if is_lab)
            if n_oh > 0:
                n_oh_per_species[sp] = n_oh
                
    n_times = C_mat.shape[1]
    non_water_oh_M = np.zeros(n_times, dtype=float)
    total_oh_M = np.zeros(n_times, dtype=float)
    
    carriers_dict = {}
    for sp, n_oh in n_oh_per_species.items():
        c_sp = C_mat[species_idx[sp]]
        oh_contrib = n_oh * c_sp
        total_oh_M += oh_contrib
        carriers_dict[f"{sp} ({n_oh} OH)"] = oh_contrib * 1000.0
        if sp != "H2O":
            non_water_oh_M += oh_contrib
            
    h2o_deduced_M = initial_h2o_molarity - 0.5 * non_water_oh_M
    h2o_direct_M = C_mat[species_idx["H2O"]] if "H2O" in species_idx else np.zeros(n_times)
    
    return {
        "h2o_deduced_M": h2o_deduced_M,
        "h2o_direct_M": h2o_direct_M,
        "total_oh_M": total_oh_M,
        "carriers_df": pd.DataFrame(carriers_dict)
    }


# Standard 29Si NMR data from Peter's DFT GIAO calculations
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
    """Convolutes time-resolved chemical concentrations into virtual 1D NMR spectra (Block 8 engine)."""
    if nmr_data is None:
        nmr_data = DEFAULT_NMR_29SI
        
    idx = {s: i for i, s in enumerate(species_list)}
    delta = np.linspace(delta_range[0], delta_range[1], n_delta)
    
    nt = len(t_s)
    snap_idx = np.unique(np.round(np.logspace(0, np.log10(nt - 1), n_snapshots)).astype(int))
    snap_idx = snap_idx[snap_idx < nt]
    t_snap_h = t_s[snap_idx] / 3600.0
    
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

