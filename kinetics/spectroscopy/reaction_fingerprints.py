"""Chemometric Reaction Fingerprints: Multi-nuclear Feature Subspaces & Extent Inversion.

Constructs multi-nuclear stoichiometric reaction fingerprints F = S · P_pure,
evaluates reaction identifiability via Singular Value Decomposition (SVD),
computes Net Analyte Signal (NAS) selectivity matrices per nucleus,
and performs inverse deconvolution of per-step reaction extents (Δξ) from
noisy or experimental operando spectra via pseudoinverse projection.
"""

import numpy as np
import pandas as pd

from kinetics.spectroscopy.molecular_symmetry import (
    DEFAULT_NUCLEI,
    build_referenced_nmr_sites
)
from kinetics.spectroscopy.multinuclear_nmr import auto_regions, lorentzian

DEFAULT_FP_LW = {"Si": 0.4, "P": 0.5, "C": 0.4, "H": 0.02}
DEFAULT_FP_PAD = {"Si": 4.0, "P": 4.0, "C": 4.0, "H": 0.4}
DEFAULT_FP_GAP = {"Si": 40.0, "P": 40.0, "C": 40.0, "H": 2.0}


def build_multinuclear_feature_space(
    visible_species: list[str],
    nmr_sites_catalog: dict = None,
    nuclei: list[str] = None,
    fp_lw: dict = None,
    fp_pad: dict = None,
    fp_gap: dict = None
) -> tuple[dict, dict, dict, int]:
    """Constructs concatenated multi-nuclear frequency grid space."""
    catalog = nmr_sites_catalog if nmr_sites_catalog is not None else build_referenced_nmr_sites()
    active_nuclei = [el for el in (nuclei or ["Si", "P", "C", "H"]) if el in catalog]
    lw_map = fp_lw or DEFAULT_FP_LW
    pad_map = fp_pad or DEFAULT_FP_PAD
    gap_map = fp_gap or DEFAULT_FP_GAP
    
    grids = {}
    for el in active_nuclei:
        shifts = [
            d for sp in visible_species
            for d, _, is_lab in catalog[el].get(sp, [])
            if not is_lab
        ]
        if not shifts:
            continue
            
        region_bounds = auto_regions(shifts, pad=pad_map[el], max_gap=gap_map[el])
        grids[el] = []
        for hi, lo in region_bounds:
            n_pts = int(np.clip(8.0 * (hi - lo) / lw_map[el], 600, 6000))
            grids[el].append((hi, lo, np.linspace(hi, lo, n_pts)))
            
    x_concat = {el: np.concatenate([x for *_, x in grids[el]]) for el in grids}
    
    block_slices = {}
    start_idx = 0
    for el in grids:
        n_feat_el = len(x_concat[el])
        block_slices[el] = slice(start_idx, start_idx + n_feat_el)
        start_idx += n_feat_el
        
    total_features = start_idx
    return grids, x_concat, block_slices, total_features


def build_pure_component_matrix(
    visible_species: list[str],
    x_concat: dict,
    block_slices: dict,
    total_features: int,
    nmr_sites_catalog: dict = None,
    fp_lw: dict = None
) -> np.ndarray:
    """Constructs pure-species spectral matrix P_pure of shape (N_species, N_features) per mM."""
    catalog = nmr_sites_catalog if nmr_sites_catalog is not None else build_referenced_nmr_sites()
    lw_map = fp_lw or DEFAULT_FP_LW
    
    p_pure_list = []
    for sp in visible_species:
        spectrum_vec = np.zeros(total_features, dtype=float)
        for el, sl in block_slices.items():
            for d, mult, is_lab in catalog.get(el, {}).get(sp, []):
                if not is_lab:
                    spectrum_vec[sl] += mult * lorentzian(x_concat[el], d, lw_map[el])
        p_pure_list.append(spectrum_vec)
        
    return np.array(p_pure_list, dtype=float)


def build_reaction_fingerprints(
    reactions_list: list[tuple[str, list[str], list[str]]],
    visible_species: list[str],
    p_pure: np.ndarray,
    block_slices: dict
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Constructs raw and block-variance normalized multi-nuclear reaction fingerprints."""
    s_vis = np.zeros((len(reactions_list), len(visible_species)), dtype=float)
    for r_idx, (_, reac, prod) in enumerate(reactions_list):
        for sp_idx, sp in enumerate(visible_species):
            s_vis[r_idx, sp_idx] = prod.count(sp) - reac.count(sp)
            
    f_raw = s_vis @ p_pure
    
    total_features = p_pure.shape[1]
    weights = np.ones(total_features, dtype=float)
    for el, sl in block_slices.items():
        norm_el = np.linalg.norm(f_raw[:, sl])
        if norm_el > 1e-12:
            weights[sl] = 1.0 / norm_el
            
    f_norm = f_raw * weights
    return f_norm, f_raw, weights, s_vis


def analyze_reaction_identifiability(
    f_norm: np.ndarray,
    f_raw: np.ndarray,
    rxn_labels: list[str],
    block_slices: dict,
    nuclei_dict: dict = None
) -> dict:
    """Performs SVD rank analysis, greedy basis selection, Net Analyte Signal (NAS), and selectivity matrix."""
    nuclei = nuclei_dict or DEFAULT_NUCLEI
    
    def matrix_rank(m, tol=1e-8):
        s = np.linalg.svd(m, compute_uv=False)
        return int((s > tol * s[0]).sum())
        
    basis_indices = []
    for i in range(len(f_norm)):
        trial_basis = basis_indices + [i]
        if matrix_rank(f_norm[trial_basis]) > len(basis_indices):
            basis_indices.append(i)
            
    fb_norm = f_norm[basis_indices]
    fb_raw = f_raw[basis_indices]
    
    a_proj = np.linalg.lstsq(fb_norm.T, f_norm.T, rcond=None)[0].T
    a_proj[np.abs(a_proj) < 1e-8] = 0.0
    a_proj = np.round(a_proj, 8)
    
    lumped_labels = []
    for j, b_idx in enumerate(basis_indices):
        members = [
            rxn_labels[r] for r in range(len(f_norm))
            if r != b_idx and abs(a_proj[r, j]) > 1e-6
        ]
        if members:
            lumped_labels.append(f"{rxn_labels[b_idx]} ⊕ {','.join(members)}")
        else:
            lumped_labels.append(rxn_labels[b_idx])
            
    def net_analyte_signal(f_mat, j):
        others = np.delete(f_mat, j, axis=0)
        coef = np.linalg.lstsq(others.T, f_mat[j], rcond=None)[0]
        return f_mat[j] - others.T @ coef
        
    active_nuclei = list(block_slices.keys())
    eval_cols = active_nuclei + ["all"]
    sel_matrix = np.full((len(basis_indices), len(eval_cols)), np.nan)
    
    for c_idx, col in enumerate(eval_cols):
        fb_sub = fb_norm if col == "all" else fb_norm[:, block_slices[col]]
        scale = np.abs(fb_sub).max()
        for j in range(len(basis_indices)):
            norm_j = np.linalg.norm(fb_sub[j])
            if norm_j > 1e-6 * scale:
                nas_j = net_analyte_signal(fb_sub, j)
                sel_matrix[j, c_idx] = np.linalg.norm(nas_j) / norm_j
                
    col_names = [nuclei.get(e, (e, e))[1] if e != "all" else "all nuclei" for e in eval_cols]
    selectivity_df = pd.DataFrame(sel_matrix, index=lumped_labels, columns=col_names)
    
    return {
        "basis_indices": basis_indices,
        "lumped_labels": lumped_labels,
        "A_projection": a_proj,
        "F_basis_norm": fb_norm,
        "F_basis_raw": fb_raw,
        "selectivity_df": selectivity_df,
        "effective_rank": len(basis_indices)
    }


def recover_reaction_extents(
    d_measured: np.ndarray,
    f_basis_norm: np.ndarray,
    weights: np.ndarray
) -> np.ndarray:
    """Recovers per-step reaction extents (Δξ in mM) via pseudoinverse projection."""
    f_basis_pinv = np.linalg.pinv(f_basis_norm)
    xi_hat = (d_measured * weights) @ f_basis_pinv
    return xi_hat
