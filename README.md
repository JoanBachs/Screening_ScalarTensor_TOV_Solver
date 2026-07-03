# Screening_ScalarTensor_TOV_Solver

Python package for solving the Tolman–Oppenheimer–Volkoff (TOV) equations in General Relativity and scalar–tensor theories with screening mechanisms.

---

## Overview

This package provides a modular framework for computing equilibrium configurations of relativistic compact stars (neutron stars and white dwarfs) in:

- General Relativity  
- Scalar–Tensor theories with screening mechanisms  
- Newtonian limits of the above systems  

The solver supports:

- Arbitrary precision integration via `mpmath`
- Multiple equations of state (tabulated and polytropic)
- Piecewise and simple polytropes
- Shooting methods for scalar field boundary conditions
- Modular implementation of screening mechanisms (currently three implemented)
- Tabulated EoS support via external data files

The code is structured as a reusable Python package and is suitable for research-level applications.

---

## Scientific Background

The code solves extensions of the Tolman–Oppenheimer–Volkoff system originally introduced in:

- R. C. Tolman, *Phys. Rev.* **55**, 364 (1939)  
- J. R. Oppenheimer and G. M. Volkoff, *Phys. Rev.* **55**, 374 (1939)

The scalar–tensor extensions implemented here allow the study of screened modified gravity models in compact stars, including chameleon-type screening mechanisms and related scalar field models.

Both relativistic and Newtonian limits are implemented, enabling controlled comparisons between gravity regimes.

---

## Project History

This project is based on the original TOV solver developed by Anton Motornenko. The original implementation is available at <https://github.com/amotornenko/TOVsolver>.

The present version substantially extends and restructures the original implementation by:

- Rewriting the integrator to support arbitrary precision arithmetic via `mpmath`
- Implementing scalar–tensor ODE systems
- Adding Newtonian counterparts
- Introducing screening mechanisms
- Implementing a new shooting method
- Supporting additional polytropic and tabulated equations of state
- Refactoring the code into a modular multi-file architecture

Only limited components related to the initialization of tabulated equations of state and crust matching are retained from the original implementation; the remainder of the codebase has been substantially rewritten and extended.

---

## Installation

Requires **Python ≥ 3.9**

Clone the repository:

```bash
git clone https://github.com/yourusername/Screening_ScalarTensor_TOV_Solver.git
cd Screening_ScalarTensor_TOV_Solver
pip install .
```

It is recommended to install inside a virtual environment.

---

## Example Usage

A minimal runnable example is provided to test the solver immediately after installation.

If the example script is located in the `examples/` directory, run:

```bash
python examples/example_white_dwarf_chameleon.py
```

This example:

- Loads a tabulated white dwarf equation of state  
- Initializes the Runge–Kutta 4 solver  
- Solves a scalar–Newtonian configuration with a chameleon scalar field  
- Demonstrates the shooting method for determining the scalar field at infinity  
- Produces output files in a structured output directory  

The example is intended as a starting point for building custom neutron star or white dwarf configurations in both GR and scalar–tensor gravity.

Users are encouraged to modify the script parameters to explore different central densities, scalar couplings, and potential parameters.

---

## Directory Structure

```
Screening_ScalarTensor_TOV_Solver/
│
├── scalar_tov_solver/        # Core solver implementation
├── eos_data/                 # Tabulated equations of state
├── examples/                 # Runnable example scripts
├── CITATION.cff
├── LICENSE
├── README.md
├── setup.py
└── setup.cfg
```

---

## Citation

If you use this software in your research, please cite it via the GitHub citation button or the provided `CITATION.cff` file.

You may also reference the original TOV formulation:

Tolman (1939)  
Oppenheimer & Volkoff (1939)

---

## License

This project is licensed under the GNU General Public License v3.0 or later (GPLv3+).

See the `LICENSE` file for details.

---

## Author

**Joan Bachs Esteban**  
Centro de Astrofísica e Gravitação (CENTRA)  
Instituto Superior Técnico – Universidade de Lisboa  
joan.bachs@tecnico.ulisboa.pt
