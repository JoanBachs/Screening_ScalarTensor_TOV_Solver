# solver.py
from .core import mpmath, numpy, pkg_resources, plt, set_precision, ZERO, ONE, TWO, MeV_fm3_to_code_units, g_cm3_to_code_units, dyn_cm2_to_code_units
from .general_relativity import tov_ode
from .scalar_tensor import chameleon, symmetron, dilaton, scalar_tov_ode
from .shooting import shoot_mtd


class ScalarTOVSolver:
    """
    High-level interface for solving Tolman-Oppenheimer-Volkoff equations
    with or without scalar-tensor screening mechanisms.
    """

    def __init__(self, en_arr, p_arr, add_crust=True, precision=100):
        """
        Initialize solver with EOS data.

        Parameters
        ----------
        en_arr : list or array
            Energy density array in MeV/fm^3.
        p_arr : list or array
            Pressure array in MeV/fm^3.
        add_crust : bool
            Whether to merge with crust EOS.
        precision : int
            Decimal precision for mpmath.
        """
        set_precision(precision)

        # Convert to code units
        self._en_arr = [mpmath.fmul(mpmath.mpf(e), MeV_fm3_to_code_units) for e in en_arr]
        self._p_arr = [mpmath.fmul(mpmath.mpf(p), MeV_fm3_to_code_units) for p in p_arr]

        # EOS interpolation
        self.en_dens = self._cubic_spline(self._p_arr, self._en_arr)
        self.press = self._cubic_spline(self._en_arr, self._p_arr)

        # Density/pressure range
        self.min_dens, self.max_dens = min(self._en_arr), max(self._en_arr)
        self.min_p, self.max_p = min(self._p_arr), max(self._p_arr)

        if add_crust:
            self._add_crust()

    # -------------------------
    # Public High-Level Methods
    # -------------------------
    def solve_gr(self, central_density, rmax, dr, atm_den):
        """
        Solve GR TOV equations for a given central density.
        """
        from .general_relativity import ode_sol
        return ode_sol(self, tov_ode, central_density, rmax, dr, atm_den)

    def solve_scalar(self, central_density, model, parameters, rmax, dr, atm_den):
        """
        Solve scalar-tensor TOV equations with a screening mechanism.
        """
        scalar_func = {"chameleon": chameleon,
                       "symmetron": symmetron,
                       "dilaton": dilaton}[model]

        return shoot_mtd(self, scalar_tov_ode, scalar_func, parameters,
                         central_density, rmax, dr, atm_den)

    def mass_radius_curve(self, densities, scalar=False, **kwargs):
        """
        Generate mass-radius curve over a list of central densities.
        """
        results = []
        for rho_c in densities:
            rho_c = mpmath.mpf(rho_c)
            if scalar:
                sol = self.solve_scalar(rho_c, **kwargs)
            else:
                sol = self.solve_gr(rho_c, **kwargs)
            results.append(sol)
        return results

    def plot_eos(self):
        """
        Plot EOS (energy density vs pressure).
        """
        plt.figure(figsize=(8, 6))
        plt.plot([float(x) for x in self._en_arr], [float(x) for x in self._p_arr], label="Original EOS")
        plt.xlabel("Energy Density")
        plt.ylabel("Pressure")
        plt.title("Equation of State")
        plt.legend()
        plt.show()

    # -------------------------
    # Public Utility Methods
    # -------------------------
    def check_density(self, dens):
        dens = mpmath.mpf(dens)
        if dens < self.min_dens:
            return self.min_dens
        elif dens > self.max_dens:
            raise ValueError(f"Density {dens} above EOS range ({self.max_dens})")
        return dens

    def check_pressure(self, P):
        P = mpmath.mpf(P)
        if P < self.min_p:
            return self.min_p
        elif P > self.max_p:
            raise ValueError(f"Pressure {P} above EOS range ({self.max_p})")
        return P

    # -------------------------
    # Internal Methods
    # -------------------------
    def _cubic_spline(self, x_data, y_data):
        n = len(x_data)
        if n < 2 or len(y_data) < 2:
            raise ValueError("Interpolation requires at least two data points.")

        h = [mpmath.fsub(x_data[i+1], x_data[i]) for i in range(n-1)]
        alpha = [mpmath.fsub(
                    mpmath.fmul(3/h[i], mpmath.fsub(y_data[i+1], y_data[i])),
                    mpmath.fmul(3/h[i-1], mpmath.fsub(y_data[i], y_data[i-1])))
                 for i in range(1, n-1)]

        l = [mpmath.mpf(1)] * n
        mu = [mpmath.mpf(0)] * n
        z = [mpmath.mpf(0)] * n

        for i in range(1, n-1):
            l[i] = mpmath.fsub(mpmath.fmul(TWO, mpmath.fsub(x_data[i+1], x_data[i-1])), mpmath.fmul(h[i-1], mu[i-1]))
            mu[i] = mpmath.fdiv(h[i], l[i])
            z[i] = mpmath.fdiv(mpmath.fsub(alpha[i-1], mpmath.fmul(h[i-1], z[i-1])), l[i])

        b = [ZERO] * n
        c = [ZERO] * n
        d = [ZERO] * n

        for j in range(n-2, -1, -1):
            c[j] = mpmath.fsub(z[j], mpmath.fmul(mu[j], c[j+1]))
            b[j] = mpmath.fsub(mpmath.fdiv(mpmath.fsub(y_data[j+1], y_data[j]), h[j]),
                           mpmath.fdiv(mpmath.fmul(h[j], mpmath.fadd(c[j+1], 2*c[j])), 3))
            d[j] = mpmath.fdiv(mpmath.fsub(c[j+1], c[j]), 3*h[j])

        def spline_func(x):
            x = mpmath.mpf(x) if not isinstance(x, mpmath.mpf) else x
            for i in range(n-1):
                if x_data[i] <= x <= x_data[i+1]:
                    dx = mpmath.fsub(x, x_data[i])
                    return mpmath.fsum([y_data[i], mpmath.fmul(b[i], dx),
                                    mpmath.fmul(c[i], mpmath.power(dx, 2)), mpmath.fmul(d[i], mpmath.power(dx, 3))])
            raise ValueError("Value out of interpolation range.")

        return spline_func

    def _add_crust(self):
        crust_file = pkg_resources.resource_filename("scalar_tov_solver.eos_data","Baym_eos.dat")
        crust_data = numpy.genfromtxt(crust_file, dtype=str, skip_header=1, names=["en", "p", "nB"])

        for col in ["en", "p", "nB"]:
            crust_data[col] = numpy.array([mpmath.mpf(x) for x in crust_data[col]], dtype=object)

        crust_data["en"] = numpy.array([mpmath.fmul(mpmath.mpf(en), g_cm3_to_code_units) for en in crust_data["en"]], dtype=object)
        crust_data["p"] = numpy.array([mpmath.fmul(mpmath.mpf(p), dyn_cm2_to_code_units) for p in crust_data["p"]], dtype=object)

        P_crust = self._cubic_spline(crust_data["en"], crust_data["p"])

        def intersection(n):
            return mpmath.fsub(P_crust(n), self.press(n))

        n_glue = mpmath.findroot(intersection, mpmath.mpf("1e-4"))

        en_arr, p_arr = [], []

        for en, p in zip(crust_data["en"], crust_data["p"]):
            if en < n_glue:
                en_arr.append(en)
                p_arr.append(p)

        for en, p in zip(self._en_arr, self._p_arr):
            if en >= n_glue:
                en_arr.append(en)
                p_arr.append(p)

        self.min_dens, self.min_p = min(en_arr), min(p_arr)
        self._en_arr, self._p_arr = en_arr, p_arr
        self.en_dens = self._cubic_spline(p_arr, en_arr)
        self.press = self._cubic_spline(en_arr, p_arr)