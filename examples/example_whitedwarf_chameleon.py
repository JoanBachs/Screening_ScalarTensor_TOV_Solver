# example_whitedwarf_chameleon.py
"""
Example: White Dwarf in Scalar–Tensor Gravity (Chameleon Model)

This script demonstrates how to:
1. Load a tabulated equation of state
2. Initialize the solver
3. Solve a scalar–tensor configuration
4. Print mass and radius

Run with:
    python example_white_dwarf_chameleon.py
"""

import numpy as np
import mpmath as mp
from pathlib import Path

from scalar_tov_solver.solver import RK4
from scalar_tov_solver.constants import g_cm3_to_code_units, Mp


def main():

    print("\nScreening_ScalarTensor_TOV_Solver")
    print("White Dwarf – Chameleon Model Example\n")

    # ------------------------------------------------------------------
    # Load Equation of State (relative to project root)
    # ------------------------------------------------------------------
    eos_path = Path(__file__).parent / "eos_data" / "WhiteDwarfEoS" / "Chandra_eos_RK4.dat"

    eos_data = np.genfromtxt(eos_path, names=True)
    rho_arr = eos_data["energy_density"]
    p_arr = eos_data["pressure"]

    # ------------------------------------------------------------------
    # Initialize Solver
    # ------------------------------------------------------------------
    solver = RK4(rho_arr, p_arr, add_crust=False, plot_eos=False)

    # ------------------------------------------------------------------
    # Model Parameters
    # ------------------------------------------------------------------
    central_density = 1e9 * g_cm3_to_code_units

    n = 1
    alpha = mp.mpf("1e-2") * Mp
    Lambda = mp.mpf("1e-3") * Mp

    scalar_par = (n, alpha, Lambda)

    # ------------------------------------------------------------------
    # Integration Parameters
    # ------------------------------------------------------------------
    max_r = mp.mpf("1e12")
    min_r = mp.mpf("1e5")
    eos_atm = mp.mpf("1e-20")

    phi_tol = mp.mpf("1e-8")
    phi_var = mp.mpf("1e-2") * Mp
    max_itr = 100

    # Initial scalar guess
    phi_inf = solver.phi_ini(
        solver.chameleon,
        scalar_par,
        1,
        phi_var,
        phi_var,
        eos_atm,
        None,
        None
    )

    # ------------------------------------------------------------------
    # Run Solver
    # ------------------------------------------------------------------
    print("Running integration...\n")

    result = solver.scalar_star_MR(
        solver.scalar_nwt_ode,
        solver.chameleon,
        scalar_par,
        np.array([central_density]),
        max_r,
        min_r,
        eos_atm,
        None,
        mp.mpf("1e-12"),
        mp.mpf("1e-12"),
        1,
        mp.mpf("0.0"),
        phi_inf,
        phi_tol,
        phi_var,
        None,
        max_itr,
        1,
        None,
        None,
        None,
        None
    )

    print("Integration finished.\n")
    print("Check output files or returned structure for results.")


if __name__ == "__main__":
    main()