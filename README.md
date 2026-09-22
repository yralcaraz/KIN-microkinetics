# KIN - Microkinetics & Reactor Dynamics

Microkinetic modeling, solution thermochemistry, and reactor dynamics simulation suite for the **BatteryAsTank** master's thesis research at Uppsala University (Department of Chemistry – Ångström Laboratory).

## Overview

This repository couples ab initio thermochemical parameters, solvation corrections, and transition state theory with macroscopic reactor dynamics to model chemical degradation and reaction kinetics in liquid-phase battery electrolyte systems.

## Project Structure

```
.
├── KIN - 260916 - Master Pipeline Document - Kinetics, Reactor Dynamics and Experimental Verification.md
├── KIN - 260918 - Architectural Comparison and Kinetic Monte Carlo Evaluation.md
├── KIN - 260918 - Kinetcs model - rev 1.ipynb       # Integrated kinetics & reactor simulation notebook
├── KIN - DRAFT - Initial modelling.ipynb           # Draft exploration notebook
├── calculate_gas_thermo.py                         # Gas-phase thermochemical evaluations
├── calculate_rate_constants.py                     # Eyring / Arrhenius rate constant derivation
├── calculate_reaction_thermo.py                    # Reaction Gibbs free energy calculations
├── calculate_solution_gibbs.py                     # Solution-phase Gibbs energy with solvation corrections
├── calculate_standard_state_shift.py               # Standard state shift (1 atm gas -> 1 M solution)
├── fit_modified_arrhenius.py                       # Modified Arrhenius parameter fitting
├── simulate_tank_reactor.py                        # Tank reactor dynamics (ODE system integration)
├── simulate_virtual_nmr.py                         # In silico NMR spectroscopy degradation tracking
├── degradation_timeseries.csv                      # Kinetics degradation time-series data
└── README.md
```

## Workflows & Methodologies

1. **Thermochemistry & Rate Parameter Generation:**
   - Calculation of standard-state corrections (\(1\text{ atm} \rightarrow 1\text{ M}\)).
   - Eyring transition-state formulation to derive forward and backward rate constants \(k(T)\).
2. **Reactor Dynamics:**
   - ODE-based non-linear microkinetic solver simulating transient concentrations of reactants, intermediates, and degradation products.
3. **Analytical Verification:**
   - Virtual NMR spectra prediction matching time-resolved species concentrations against experimental tracking data.

## Author & Academic Affiliation

- **Author:** Yeray Alcaraz Galván
- **Supervision:** Prof. Peter Broqvist
- **Institution:** Department of Chemistry – Ångström Laboratory, Uppsala University
