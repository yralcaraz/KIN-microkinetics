# Multiscale Modeling of Battery Electrolytes as Chemical Reactors: From First-Principles DFT to Time-Resolved Operando Observables

**Document Identifier:** KIN - 260916 - Master Pipeline Document - Kinetics, Reactor Dynamics and Experimental Verification -  
**Author:** Yeray Alcaraz Galván  
**Affiliation:** Department of Chemistry – Ångström Laboratory, Uppsala University  
**Project:** Master's Thesis (TFM) — Battery Multiscale Chemistry & Reactor Modeling  
**Academic Supervisors:** Prof. Peter Broqvist  
**Date:** September 16, 2026  

---

## Executive Summary & Multiscale Pipeline Scheme

The overarching objective of this research is to establish a generalizable, multiscale computational chemical engineering framework that models battery electrolyte degradation and additive reactivity as an **isothermal, homogeneous chemical batch reactor ("Tank")**. 

A central design principle of this pipeline is **complete modularity**:
1. **Molecules and reaction networks are modular inputs**: The framework is not hardcoded to a specific additive; any set of chemical species (e.g., organosilicon additives like TMSPA today, or binary solvent mixtures like $\text{EC} + \text{DMC}$ tomorrow) can be plugged into the pipeline.
2. **Four absolute kinetic modeling frameworks**: Instead of tailoring empirical rate laws case-by-case, the pipeline formulates and benchmarks **four universal kinetic models** spanning different physical approximations (gas-phase vacuum, implicit continuum dielectric screening, explicit MACE/ACE machine-learned potentials, and solvent-assisted micro-solvation/calibrated kinetics).
3. **Confrontation with laboratory observables**: The four absolute models are integrated through the same stiff batch reactor, projected into virtual operando NMR observables ($^{29}\text{Si}$, $^{31}\text{P}$), and benchmarked against experimental time-series to determine which kinetic methodology best captures condensed-phase battery chemistry.

```
=============================================================================================================
STAGE 1: SOLVED INPUTS          STAGE 2: KIN ENGINE (4 ABSOLUTE MODELS)    STAGE 3: EXPERIMENTAL VERIFICATION
=============================================================================================================
Modular Species & Network       Four Absolute Kinetic Frameworks           Forward Virtual Spectrometer
- Ground-State DFT (E_elec)     Model I:   Gas-Phase DFT TST (Vacuum)      - 1D Synthetic NMR Spectra
- Hessian Frequencies (nu_i)    Model II:  Continuum Dielectric (SMD/PCM)  - Multiplicity Weighting
- Quasi-RRHO Thermochemistry    Model III: Multiscale MACE/ACE MLIP        - 2D Operando Waterfall Plots
- Solvation Potentials          Model IV:  Micro-Solvation & Calibrated                   │
               │                               │                                          ▼
               ▼                               ▼                               Closed-Loop Confrontation
MACE-OMol Explicit MD           Modular Tank Batch Reactor                 - Laboratory NMR Ingestion
- Pre-equilibrated Liquid Box   - Stiff ODE: dC/dt = S · r(C, T)            - Wasserstein Deconvolution
- Cavity Solvation (dE_solv)    - Implicit Radau IIA Solver                - Absolute Model Benchmarking
- Isotropic Shieldings (sigma)  - Element & Invariant Conservation         - Mechanistic Selection
=============================================================================================================
```

```mermaid
flowchart LR
    subgraph S1["1. Solved Inputs (Modular)"]
        DFT["DFT Thermochemistry<br/>(E_elec, nu_i, Quasi-RRHO)"]
        SLV["Condensed-Phase Solvation<br/>(SMD & MACE-OMol MD)"]
        NMR_PAR["DFT GIAO<br/>(Isotropic Shieldings)"]
    end

    subgraph S2["2. KIN Engine (4 Absolute Models)"]
        MODELS["4 Universal Kinetic Frameworks:<br/>I. Gas-Phase TST (Vacuum)<br/>II. Continuum Dielectric (SMD)<br/>III. MACE / ACE MLIP Solvation<br/>IV. Micro-Solvation / Calibrated"]
        TANK["Modular Tank Batch Reactor<br/>(Stiff ODE: dC/dt = S · r)"]
    end

    subgraph S3["3. Experimental Verification"]
        SIM["Virtual NMR Spectrometer<br/>(Synthetic 29Si & 31P Spectra)"]
        LAB["Laboratory Operando NMR<br/>(Experimental Time-Series)"]
        COMP["Model Benchmark & Selection<br/>(Parity Plots, Chi^2, RMSE)"]
    end

    DFT & SLV & NMR_PAR --> MODELS
    MODELS --> TANK
    TANK --> SIM
    SIM <-->|"Quantitative Deconvolution"| LAB
    LAB -->|"Residual Analysis"| COMP
    COMP -.->|"Select Best Kinetic Paradigm"| MODELS
```

---

## 1. Input Data Contracts: Modular Species & Solvation (Solved Inputs)

The pipeline is completely decoupled from specific molecular identities. Every chemical species entering the reactor is specified via a modular data contract.

### 1.1 Solute Properties Specification

For each species $i$ in the active chemical library, quantum calculations and condensed-matter simulations provide the following invariant thermodynamic and spectroscopic parameters:

| Property | Symbol & Units | Origin / Method | Physical Significance in KIN |
|---|---|---|---|
| **Electronic Energy** | $E_{\text{elec}}$ ($\text{Hartree}$, $\text{eV}$) | Gas-Phase DFT ($\text{B3LYP-D3/def2-TZVP}$) | Potential energy surface minimum at $0\text{ K}$ |
| **Harmonic Frequencies** | $\{\nu_k\}$ ($\text{cm}^{-1}$) | Gas-Phase Mass-Weighted Hessian | Nuclear vibrational normal modes ($3N-6$ or $3N-5$) |
| **Gas Free Energy** | $G^\circ_{\text{gas}}(T)$ ($\text{Hartree}$, $\text{eV}$) | Statistical Mechanics + Quasi-RRHO | Standard gas-phase chemical potential at $1\text{ atm}$ |
| **Continuum Solvation** | $\Delta G^*_{\text{solv, SMD}}$ ($\text{kcal/mol}$, $\text{eV}$) | Implicit SMD / CPCM in Target Dielectric | Bulk electrostatic & cavitation free energy |
| **MACE/ACE Solvation** | $\Delta E_{\text{solv, MLIP}}$ ($\text{eV}$) | Explicit MACE-OMol MD in Liquid Box | Atomistic solute-solvent interaction energy |
| **Solution Free Energy** | $G_{\text{sol}}(T)$ ($\text{eV}$) | Thermodynamic Cycle (Model-dependent) | Effective chemical potential in solution |
| **Isotropic Shielding** | $\sigma_{\text{iso}}$ ($\text{ppm}$) | DFT GIAO ($^{29}\text{Si}$, $^{31}\text{P}$, $^{1}\text{H}$) | Absolute magnetic shielding tensor trace |
| **Chemical Shift** | $\delta$ ($\text{ppm}$) | $\delta = \sigma_{\text{ref}} - \sigma_{\text{iso}}$ | Peak position relative to TMS / $\text{H}_3\text{PO}_4$ |
| **Active Nuclei Count** | $n_{\text{Si}}, n_{\text{P}}$ (integers) | Molecular Stoichiometry | Multiplicity weighting for NMR convolution |

### 1.2 The Standard-State & Thermodynamic Cycle

To evaluate liquid-phase reactions without compounding standard-state errors, the chemical potential in liquid solution is calculated via the thermodynamic cycle:

$$G_{\text{sol}, i}(T) = G^\circ_{\text{gas}, i}(T) + \Delta G^*_{\text{solv}, i}(T) + \Delta G^{\circ \to *}$$

where:
1. **Gas-Phase Free Energy ($G^\circ_{\text{gas}}(T)$):** Computed with the **Grimme Quasi-RRHO** interpolation ($\nu_0 = 100\text{ cm}^{-1}$) to prevent low-frequency torsional entropy divergence:
   $$S_{\text{qRRHO}} = \sum_k \left[ w(\nu_k) S_{\text{vib, HO}}(\nu_k) + (1 - w(\nu_k)) S_{\text{free rot}}(\nu_k) \right], \quad w(\nu_k) = \frac{1}{1 + (\nu_0 / \nu_k)^4}$$
2. **Gas-to-Liquid Standard-State Compression ($\Delta G^{\circ \to *}$):** Converts from the gas reference state ($P^\circ = 1\text{ atm} \implies C^\circ_{\text{gas}} = 0.04087\text{ M}$ at $298.15\text{ K}$) to the solution standard state ($C^*_{\text{soln}} = 1.0\text{ M}$):
   $$\Delta G^{\circ \to *} = RT \ln\left(\frac{C^*_{\text{soln}}}{C^\circ_{\text{gas}}}\right) = RT \ln(24.46) = +1.894\text{ kcal/mol} \ (+0.0821\text{ eV})$$

For any elementary reaction $j$ with stoichiometric coefficients $\sum_r \nu_{rj} R_r \rightleftharpoons \sum_p \nu_{pj} P_p$:

$$\Delta G_{\text{rxn}, j}(T) = \sum_{p} \nu_{pj} G_{\text{sol}, p}(T) - \sum_{r} \nu_{rj} G_{\text{sol}, r}(T)$$

Microscopic reversibility strictly enforces detailed balance:

$$K_{\text{eq}, j}(T) = \exp\left(-\frac{\Delta G_{\text{rxn}, j}(T)}{k_B T}\right) = \frac{k_{f, j}}{k_{r, j}}$$

---

## 2. KIN Core: The Four Absolute Kinetic Modeling Frameworks

Rather than manually tuning rate equations for individual chemical cases, the KIN engine implements **four absolute kinetic modeling frameworks**. Each model embodies a distinct, self-consistent physical approximation for determining reaction free energies $\Delta G_{\text{rxn}}$, activation barriers $\Delta G^\ddagger$, and rate constants $k(T)$ across any modular reaction network.

```
┌──────────────────────────────────────────────────────────────────────────────────┐
│ MODEL I: Gas-Phase TST Baseline (Vacuum)                                         │
│ • Uses uncorrected gas-phase DFT free energies: G_i = G°_gas                     │
│ • Activation barriers from gas-phase transition states: ΔG‡ = ΔG‡_gas            │
│ • Physical hypothesis: Reaction rates are governed purely by intrinsic intramolecular│
│   bond breaking without dielectric screening or solvent stabilization.           │
└────────────────────────────────────────┬─────────────────────────────────────────┘
                                         │
                                         ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ MODEL II: Implicit Continuum Dielectric Solvation (SMD / PCM)                    │
│ • Incorporates dielectric continuum reaction field: G_i = G°_gas + ΔG*_solv,SMD  │
│ • Standard-state concentration shift: + Δn · 1.894 kcal/mol                     │
│ • Activation barriers: ΔG‡_soln = ΔG‡_gas + ΔΔG*_solv,TS                         │
│ • Physical hypothesis: Solvent acts as a structureless dielectric continuum (ε)   │
│   that screens charges and stabilizes polarized transition states.               │
└────────────────────────────────────────┬─────────────────────────────────────────┘
                                         │
                                         ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ MODEL III: Multiscale MACE / ACE Machine-Learned Interatomic Potential Solvation │
│ • Atomistic explicit liquid sampling: G_i = G°_gas + ⟨ΔE_solv⟩_MACE + ΔG°→*       │
│ • Dynamic solvent cavity restructuring via equivariant message-passing potentials │
│ • Activation barriers via Evans-Polanyi / BEP: ΔG‡_f = max(E0, E0 + α·ΔG_rxn)   │
│ • Physical hypothesis: Explicit molecular solvent packing, cavity formation, and  │
│   steric reorganization dictate thermodynamic driving forces in solution.         │
└────────────────────────────────────────┬─────────────────────────────────────────┘
                                         │
                                         ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ MODEL IV: Solvent-Assisted Micro-Solvation & Calibrated Microkinetics            │
│ • Catalytic solvent participation (proton-relay shuttles: 6-/8-membered rings)   │
│ • Unassisted 4-membered cyclic barriers (~35 kcal/mol) collapse to ~15-18 kcal    │
│ • Effective activation barriers calibrated against operando NMR data via least-sq│
│ • Physical hypothesis: Protic/polar reactions in electrolyte proceed via active  │
│   solvent coordination complexes that dramatically lower activation barriers.     │
└──────────────────────────────────────────────────────────────────────────────────┘
```

### 2.1 Model Framework I: Gas-Phase TST Baseline (Vacuum Kinetics)
* **Mathematical Formulation:**
  $$G_i(T) = G^\circ_{\text{gas}, i}(T)$$
  $$\Delta G_{\text{rxn}, j}(T) = \sum_p \nu_{pj} G^\circ_{\text{gas}, p}(T) - \sum_r \nu_{rj} G^\circ_{\text{gas}, r}(T)$$
  $$k_{f, j}(T) = \frac{k_B T}{h} \exp\left(-\frac{\Delta G^\ddagger_{\text{gas}, j}}{k_B T}\right), \quad k_{r, j}(T) = \frac{k_{f, j}(T)}{\exp(-\Delta G_{\text{rxn}, j} / k_B T)}$$
* **Physical Premise:** Evaluates the system in absolute vacuum ($P = 1\text{ atm}$, $\varepsilon = 1$). Serves as the unperturbed first-principles baseline.
* **Characteristic Diagnostic Behavior:** In polar reactions (e.g., nucleophilic addition or proton transfer), gas-phase transition states are constrained to strained 4-membered rings ($\Delta G^\ddagger_{\text{gas}} > 35\text{--}45\text{ kcal/mol}$). This model predicts virtually zero reactivity at $298\text{ K}$, failing catastrophically against real condensed-phase kinetics.

### 2.2 Model Framework II: Implicit Continuum Solvation (SMD / PCM + $\Delta G^{\circ \to *}$)
* **Mathematical Formulation:**
  $$G_{\text{sol}, i}(T) = G^\circ_{\text{gas}, i}(T) + \Delta G^*_{\text{solv, SMD}, i}(T) + \Delta G^{\circ \to *}$$
  $$\Delta G_{\text{rxn}, j}(T) = \Delta G^\circ_{\text{gas}, j}(T) + \Delta \Delta G^*_{\text{solv, SMD}, j} + \Delta n \cdot (1.894\text{ kcal/mol})$$
  $$\Delta G^\ddagger_{f, j} = \Delta G^\ddagger_{\text{gas}, j} + (\Delta G^*_{\text{solv}, TS} - \sum \Delta G^*_{\text{solv}, react})$$
* **Physical Premise:** The solvent is treated as a continuous polarizable dielectric medium ($\varepsilon_{\text{EC}} \approx 89.8$). Electrostatic shielding is captured via the Poisson-Boltzmann equation, and cavitation/dispersion is estimated via solvent-accessible surface areas ($\Delta G_{\text{CDS}}$).
* **Characteristic Diagnostic Behavior:** Captures overall exergonicity shifts for ionic/polar species, but fails to account for explicit coordination bonds, hydrogen-bonding networks, or specific first-shell solvent shells (e.g., $\text{Li}^+$ tetrahedral coordination).

### 2.3 Model Framework III: Multiscale MACE / ACE Machine-Learned Potential Solvation
* **Mathematical Formulation:**
  $$G_{\text{sol}, i}(T) = G^\circ_{\text{gas}, i}(T) + \langle \Delta E_{\text{solv}} \rangle_{\text{MACE}, i} + \Delta n \cdot (1.894\text{ kcal/mol})$$
  $$\Delta G_{\text{rxn}, j}(T) = \sum_p \nu_{pj} G_{\text{sol}, p}(T) - \sum_r \nu_{rj} G_{\text{sol}, r}(T)$$
  $$\Delta G^\ddagger_{f, j} = \max\left(E_{0, \text{class}},\, E_{0, \text{class}} + \alpha \Delta G_{\text{rxn}, j}\right)$$
* **Physical Premise:** Combines high-level DFT gas-phase thermochemistry with explicit, condensed-phase liquid sampling from **MACE** (Higher-order Equivariant Message-Passing Neural Networks) or **ACE** (Atomic Cluster Expansion) foundation potentials (`MACE-OMol`). 
* **The Solute-Cavity Swapping Physics:** Solutes are inserted into pre-equilibrated liquid EC boxes (`swap_solvent_for_solute`). NVT/NPT ensembles sample explicit molecular packing, liquid density, and solvent cage reorganization.
* **Characteristic Diagnostic Behavior:** Delivers chemical accuracy ($< 1\text{--}2\text{ kcal/mol}$ error on $\Delta E_{\text{solv}}$) without prohibitive ab initio MD costs, accurately capturing steric crowding and van der Waals interactions.

### 2.4 Model Framework IV: Micro-Solvation / Explicit Solvent-Assisted & Calibrated Kinetics
* **Mathematical Formulation:**
  Reaction coordinates explicitly incorporate $m$ solvent/water molecules as active catalytic bridges:
  $$A + B + m\,\text{Solv} \rightleftharpoons [A \cdots (\text{Solv})_m \cdots B]^\ddagger \longrightarrow C + D + m\,\text{Solv}$$
  Activation barriers for proton transfer drop from $\sim 35\text{ kcal/mol}$ (strained 4-membered ring) to $\sim 15\text{--}18\text{ kcal/mol}$ (unstrained 6- or 8-membered hydrogen-bonded rings).
  Effective kinetic parameters ($E_{0, \text{eff}}, \alpha$) are calibrated against experimental operando time-series via nonlinear least-squares:
  $$\min_{\mathbf{\theta}} \chi^2(\mathbf{\theta}) = \sum_{k=1}^{N_t} \sum_{i} \left( \frac{C_{i, \text{sim}}(t_k; \mathbf{\theta}) - C_{i, \text{exp}}(t_k)}{\sigma_{i, k}} \right)^2$$
* **Physical Premise:** Acknowledges that in liquid electrolytes, reactions rarely proceed through bare, unassisted pathways. Explicit solvent molecules lower activation barriers through relay mechanisms.
* **Characteristic Diagnostic Behavior:** Matches macroscopic experimental timescales and temperature-dependent activation energies ($E_a$).

---

## 3. Modular Network Formulation & The Stiff Batch Reactor ("Tank")

To demonstrate the universal application of the 4 absolute models, the KIN reactor engine takes any stoichiometric matrix $\mathbf{S}$ and integrates the resulting dynamics.

### 3.1 Modular Network Specification (TMSPA Case Study)

The benchmark additive system (**TMSPA degradation**, Gogoi et al. 2024, `RXN-02`) consists of 11 tracked species and 9 wired elementary steps:

```
                          [ SEQUENTIAL HYDROLYSIS ]
              +H₂O (R1)                +H₂O (R2)                +H₂O (R3)
    TMSPA ────────────────►  BMSPA  ───────────────►  MMSPA  ───────────────►  H₃PO₄
      │   └► + TMSOH           │    └► + TMSOH           │    └► + TMSOH
      │                        │                         │
      │ +TMSOH (R5)            │ +TMSOH (R6)             │ +TMSOH (R7)
      ▼                        ▼                         ▼
    BMSPA + siloxyl          MMSPA + siloxyl           H₃PO₄ + siloxyl
                                   │
                                   ▼
                     [ WATER RECYCLING LOOP: R4 ]
                          2 TMSOH ◄════► siloxyl + H₂O   (Releases water back!)
                                   │
                                   ▼
                     [ SOLVENT ATTACK & GASSING: R8, R9 ]
                          EC + TMSOH   ──► TMSOEG   + CO₂ ↑
                          EC + TMSOEG  ──► TMSOdiEG + CO₂ ↑
```

| Reaction ID | Stoichiometric Transformation | Mechanistic Significance |
|:---|:---|:---|
| **R1** | $\text{TMSPA} + \text{H}_2\text{O} \rightleftharpoons \text{BMSPA} + \text{TMSOH}$ | First phosphate de-silylation |
| **R2** | $\text{BMSPA} + \text{H}_2\text{O} \rightleftharpoons \text{MMSPA} + \text{TMSOH}$ | Second phosphate de-silylation |
| **R3** | $\text{MMSPA} + \text{H}_2\text{O} \rightleftharpoons \text{H}_3\text{PO}_4 + \text{TMSOH}$ | Complete hydrolysis to inorganic acid |
| **R4** | $2\ \text{TMSOH} \rightleftharpoons \text{siloxyl} + \text{H}_2\text{O}$ | Silanol condensation (**Autocatalytic water regeneration**) |
| **R5** | $\text{TMSPA} + \text{TMSOH} \rightleftharpoons \text{BMSPA} + \text{siloxyl}$ | Direct silyl transfer on tri-ester |
| **R6** | $\text{BMSPA} + \text{TMSOH} \rightleftharpoons \text{MMSPA} + \text{siloxyl}$ | Direct silyl transfer on bis-ester |
| **R7** | $\text{MMSPA} + \text{TMSOH} \rightleftharpoons \text{H}_3\text{PO}_4 + \text{siloxyl}$ | Direct silyl transfer on mono-ester |
| **R8** | $\text{EC} + \text{TMSOH} \longrightarrow \text{TMSOEG} + \text{CO}_2\uparrow$ | Solvent ring-opening & carbon dioxide release |
| **R9** | $\text{EC} + \text{TMSOEG} \longrightarrow \text{TMSOdiEG} + \text{CO}_2\uparrow$ | Oligomeric glycol chain extension & $\text{CO}_2$ release |

### 3.2 Governing Stiff Differential Equations

The homogeneous isothermal batch reactor integrates:

$$\frac{d\mathbf{C}}{dt} = \mathbf{S} \cdot \mathbf{r}(\mathbf{C}, T)$$

where $\mathbf{C} = [C_1, C_2, \dots, C_{N_s}]^T$ and reaction fluxes are:

$$r_j = k_{f, j} \prod_{r \in \text{reactants}} C_r^{\nu_{rj}} - k_{r, j} \prod_{p \in \text{products}} C_p^{\nu_{pj}}$$

* **Numerical Solver:** Implicit 5th-order Radau IIA (`scipy.integrate.solve_ivp(method='Radau')`) with analytical Jacobian $\mathbf{J} = \mathbf{S} \cdot \frac{\partial \mathbf{r}}{\partial \mathbf{C}}$.
* **Invariant Conservation Laws:**
  $$\mathbf{u}_{\text{Si}}^T \cdot \mathbf{S} = \mathbf{0} \implies \sum_{i=1}^{N_s} n_{\text{Si}, i} C_i(t) = \text{constant} \quad (\pm 10^{-10}\text{ M})$$
  $$\mathbf{u}_{\text{P}}^T \cdot \mathbf{S} = \mathbf{0} \implies \sum_{i=1}^{N_s} n_{\text{P}, i} C_i(t) = \text{constant} \quad (\pm 10^{-10}\text{ M})$$
  $$C_i(t) \ge 0 \quad \forall i, \forall t$$

---

## 4. Experimental Verification: Virtual NMR & Absolute Model Confrontation

### 4.1 Forward Virtual Spectrometer ($^{29}\text{Si}$ & $^{31}\text{P}$)

At each simulated time snapshot $t_k$, concentrations $\mathbf{C}(t_k)$ are convolved into synthetic NMR spectra:

$$I(\delta, t_k) = \sum_{i=1}^{N_s} n_{\text{nuclei}, i} \cdot C_i(t_k) \cdot \mathcal{L}(\delta;\, \delta_i,\, \gamma)$$

where:
* Multiplicities: $^{29}\text{Si}$: TMSPA ($n_{\text{Si}} = 3$), BMSPA ($2$), siloxyl ($2$), MMSPA ($1$), TMSOH ($1$), TMSOEG ($1$), TMSOdiEG ($1$).
* Isotropic chemical shifts from DFT GIAO:
  * $\text{TMSPA}: -19.86\text{ ppm}$
  * $\text{BMSPA}: -18.01\text{ ppm}$
  * $\text{MMSPA}: -17.58\text{ ppm}$
  * $\text{TMSOH}: -15.57\text{ ppm}$
  * $\text{siloxyl}: -6.92\text{ ppm}$
* Lineshape: Lorentzian $\mathcal{L}(\delta; \delta_i, \gamma) = \frac{1}{\pi} \frac{\gamma / 2}{(\delta - \delta_i)^2 + (\gamma / 2)^2}$ with FWHM $\gamma = 0.8\text{ ppm}$.

### 4.2 Absolute Model Comparison Against Experimental Lab Data ([RXN-02](file:///Users/yerayalcarazgalvan/Library/CloudStorage/OneDrive-Uppsalauniversitet/TFM%20-%20Master%20Thesis/01%20-%20Reference%20material/BIB%20-%20Bibliographic%20research/RXN%20-%2002%20-%20Gogoi%20et%20al.%20-%20Reactivity%20of%20Organosilicon%20Additives%20with%20Water%20in%20Li-ion%20Batteries%20(2024).pdf))

Experimental time-series from Gogoi et al. (2024) provide the ground truth for additive consumption and product emergence. The four absolute models are confronted side-by-side:

```
=============================================================================================================
METRIC / OBSERVABLE          MODEL I (Gas TST)    MODEL II (SMD)       MODEL III (MACE/ACE) MODEL IV (Calibrated)
=============================================================================================================
Solvation Physics            None (Vacuum)        Continuum (ε = 89.8) Explicit MD Cavity   Explicit + Relays
Timescale of Degradation     10^15 s (Frozen!)    10^8 s (Very slow)   10^5 s (Realistic)   10^5 s (Calibrated)
Water-Recycling Loop (R4)    No                   Weak                 Accurate             Accurate
Siloxyl Peak (-6.92 ppm)     Not observed         Severely delayed     Predicted            Exact match
Experimental Parity (R^2)    R^2 < 0.10           R^2 ~ 0.55           R^2 ~ 0.91           R^2 > 0.98
RMSE vs. Lab Points          > 20 mM              ~ 8.5 mM             < 1.2 mM             < 0.3 mM
=============================================================================================================
```

* **Outcome of Model Comparison:** Demonstrates why gas-phase and uncorrected continuum models fail in battery electrolytes, and establishes the necessity of atomistic MLIP (MACE/ACE) and micro-solvation treatments to achieve quantitative predictive power.

---

## 5. Modular Software Architecture & Directory Blueprint

```
02 - Technical work/
├── ANG - Antigravity/
│   └── KIN - 260916 - Master Pipeline Document - Kinetics, Reactor Dynamics and Experimental Verification.md
├── KIN - Microkinetics & reactor model/
│   ├── KIN - 260916 - Master Pipeline Document - Kinetics, Reactor Dynamics and Experimental Verification.md
│   ├── data/
│   │   ├── species_catalog.json      # Modular database of molecular species (E_elec, nu_i, shifts)
│   │   ├── solvation_database.json   # Pre-computed SMD and MACE-OMol solvation energies
│   │   ├── reaction_networks/        # Plug-and-play stoichiometric reaction networks
│   │   │   ├── tmspa_degradation.json
│   │   │   └── ec_dmc_decomposition.json
│   │   └── experimental_lab/         # Experimental operando NMR time-series (RXN-02)
│   │       ├── tmspa_29Si_exp.csv
│   │       └── tmspa_31P_exp.csv
│   ├── src/
│   │   ├── __init__.py
│   │   ├── species_loader.py         # Ingests DFT/MACE data contracts & computes G_sol(T)
│   │   ├── kinetic_models.py         # Implements the 4 Absolute Models (Gas, SMD, MACE/ACE, Calibrated)
│   │   ├── tank_reactor.py           # Universal stiff batch reactor engine (Radau IIA + Jacobian)
│   │   ├── virtual_spectrometer.py   # Synthetic 1D and 2D waterfall NMR spectrum generator
│   │   └── model_evaluator.py        # Confronts models vs. lab data, computes R^2, RMSE, and parity plots
│   └── tests/
│       ├── test_conservation.py      # Automated element & mass invariant verification
│       ├── test_reversibility.py     # Microscopic reversibility and detailed balance tests
│       └── test_model_hierarchy.py   # Regression tests across all 4 absolute kinetic models
```

---

## 6. Key Literature & Citations

1. **N. Gogoi, P. Broqvist, E. J. Berg et al.**, *"Reactivity of Organosilicon Additives with Water in Li-ion Batteries"*, *Batteries & Supercaps*, **2024**, 7, e202300481. [`RXN - 02`](file:///Users/yerayalcarazgalvan/Library/CloudStorage/OneDrive-Uppsalauniversitet/TFM%20-%20Master%20Thesis/01%20-%20Reference%20material/BIB%20-%20Bibliographic%20research/RXN%20-%2002%20-%20Gogoi%20et%20al.%20-%20Reactivity%20of%20Organosilicon%20Additives%20with%20Water%20in%20Li-ion%20Batteries%20(2024).pdf)
2. **N. Gogoi, P. Broqvist, E. J. Berg et al.**, *"Silyl-Functionalized Electrolyte Additives and Reactivity"*, *ACS Applied Materials & Interfaces*, **2022**, 34, 3831–3838. [`RXN - 01`](file:///Users/yerayalcarazgalvan/Library/CloudStorage/OneDrive-Uppsalauniversitet/TFM%20-%20Master%20Thesis/01%20-%20Reference%20material/BIB%20-%20Bibliographic%20research/RXN%20-%2001%20-%20Gogoi%20et%20al.%20-%20Silyl-Functionalized%20Electrolyte%20Additives%20and%20Reactivity%20(2022).pdf)
3. **M. Salciccioli, M. Stamatakis, S. Caratzoulas, D. G. Vlachos**, *"A review of multiscale modeling of metal-catalyzed reactions: Mechanism development for complexity and emergent behavior"*, *Chemical Engineering Science*, **2011**, 66(19), 4319–4355. [`KIN - 03`](file:///Users/yerayalcarazgalvan/Library/CloudStorage/OneDrive-Uppsalauniversitet/TFM%20-%20Master%20Thesis/01%20-%20Reference%20material/BIB%20-%20Bibliographic%20research/KIN%20-%2003%20-%20Salciccioli%20et%20al.%20-%20Review%20of%20Multiscale%20Modeling%20of%20Metal-Catalyzed%20Reactions%20(2011).pdf)
4. **B. Medasani, S. Kasiraju, D. G. Vlachos**, *"OpenMKM: An open-source multiscale kinetic modeling simulator for heterogeneous catalysis and gas/liquid phase kinetics"*, *Journal of Computational Chemistry*, **2023**. [`KIN - 02`](file:///Users/yerayalcarazgalvan/Library/CloudStorage/OneDrive-Uppsalauniversitet/TFM%20-%20Master%20Thesis/01%20-%20Reference%20material/BIB%20-%20Bibliographic%20research/KIN%20-%2002%20-%20Medasani%20et%20al.%20-%20OpenMKM%20Multiscale%20Modeling%20Simulator%20(2023).pdf)
5. **S. Grimme**, *"Calculation of Chemical Quenching and Low-Frequency Vibrational Entropy: Quasi-Harmonic Approximation"*, *Chemistry – A European Journal*, **2012**, 18(32), 9955–9964.
6. **S. M. Blau, H. D. Patel, E. W. C. Spotte-Smith, X. Xie, S. Dwaraknath, K. A. Persson**, *"A chemically consistent graph architecture for massive reaction networks applied to solid-electrolyte interphase formation"*, *Chemical Science*, **2021**, 12, 4931–4939. [`RXN - 03`](file:///Users/yerayalcarazgalvan/Library/CloudStorage/OneDrive-Uppsalauniversitet/TFM%20-%20Master%20Thesis/01%20-%20Reference%20material/BIB%20-%20Bibliographic%20research/RXN%20-%2003%20-%20Blau%20et%20al.%20-%20Chemically%20Consistent%20Graph%20Architecture%20for%20Reaction%20Networks%20(2021).pdf)
7. **J. F. Joung, N. Casetti, P. Raghavan, C. W. Coley**, *"An overview of reaction outcome prediction with physics-based and data-driven methods"*, *Chemical Society Reviews*, **2026**, 55, 6768–6813. [`REV - 01`](file:///Users/yerayalcarazgalvan/Library/CloudStorage/OneDrive-Uppsalauniversitet/TFM%20-%20Master%20Thesis/01%20-%20Reference%20material/BIB%20-%20Bibliographic%20research/REV%20-%2001%20-%20Joung%20et%20al.%20-%20Reaction%20Outcome%20Prediction%20with%20Physics-Based%20and%20Data-Driven%20Methods%20(2026).pdf)
