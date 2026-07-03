# shooting.py
from .core import *

"""
Shooting Method Algorithms
"""


def phi_ini(self, scalar, scalar_par, method, phi_var, phi_inf, c_den, EoS_param, file):
    eden = self.E_of_rho_general(c_den, EoS_param)
    P = self.P_of_rho_general(c_den, EoS_param)
    # Energy-momentum tensor trace
    T = mpmath.fsub(mpmath.fmul(THREE, mpmath.mpf(P)), mpmath.mpf(eden))

    if method == 0:
        phi_try = mpmath.mpf(phi_var)

    elif method == 1:
        if scalar == self.chameleon:
            n, a, L = scalar_par

            if T < 0:
                phi_try = mpmath.power(mpmath.fdiv(mpmath.fmul(n, mpmath.power(L, mpmath.fadd(
                    n, FOUR))), mpmath.fneg(mpmath.fmul(a, T))), mpmath.fdiv(ONE, mpmath.fadd(n, ONE)))
            else:
                phi_try = mpmath.fmul(phi_inf, HALF)

        elif scalar == self.symmetron:
            mu, lmbd, Mcpl = scalar_par
            numerator = mpmath.fadd(mpmath.fprod([mu, Mcpl, mu, Mcpl]), T)

            if numerator < 0:
                # No symmetry breaking, the VEV goes to zero
                phi_try = ZERO
            else:
                # The potential breaks reflection symmetry spontaneously
                phi_try = mpmath.sqrt(mpmath.fdiv(
                    numerator, mpmath.fprod([lmbd, Mcpl, Mcpl])))

                if phi_try > Mcpl:
                    # Get the order of magnitude
                    exponent = mpmath.floor(
                        mpmath.log10(mpmath.fdiv(phi_try, Mcpl)))
                    # Set phi_try to the same order
                    phi_try = mpmath.fdiv(
                        phi_try, mpmath.power(10, exponent + 0.5))

                    if phi_try > Mcpl:
                        # Compute the fractional part
                        fractional_part = mpmath.fdiv(
                            phi_try, Mcpl) - mpmath.floor(mpmath.fdiv(phi_try, Mcpl))
                        # Update phi_try with only the fractional part
                        phi_try = mpmath.fmul(fractional_part, Mcpl)

        elif scalar == self.dilaton:
            Vo, a2, fd = scalar_par
            phi_try = mpmath.fadd(mpmath.fmul(Mp, mpmath.fdiv(
                Vo, mpmath.fmul(a2, mpmath.fsub(mpmath.fmul(FOUR, Vo), T)))), fd)
            if phi_try < 0:
                phi_try = mpmath.fmul(phi_try, mpmath.fneg(2))

        else:
            raise ValueError("No scalar field was selected")

    else:
        raise ValueError("No convergence method was selected")

    # Write debug information to the file
    if file is not None:
        file.write(
            f"Initial guess for c_den={mpmath.nstr(mpmath.fdiv(mpmath.mpf(c_den), mpmath.mpf(g_cm3_to_code_units)), 3)} g/cm3: phi_try={mpmath.nstr(mpmath.fdiv(mpmath.mpf(phi_try), mpmath.mpf(Mp)), 10)} Mp (method={method})\n")

    return phi_try


def shoot_mtd(self, ode_sys, scalar, scalar_par, m, max_itr, phi_tol, phiv, dens_c, phi0, phi_inf, sigma0, rmax, dr, atm, EoS, m_rel, p_tol, crit, gamma, file, log_dir):
    """
    Shooting method to solve the TOV equations with a scalar field contribution.
    Includes fixes for numerical stability, unphysical scalar values, and improved logging.
    """
    def log_message(message, target_file):
        """
        Logs a message efficiently, avoiding redundant file opens
        """
        try:
            if isinstance(target_file, io.TextIOWrapper):  # Directly write to file object
                target_file.write(message + "\n")
                target_file.flush()
            elif isinstance(target_file, str) and log_dir is not None:  # Open file only once
                with open(os.path.join(log_dir, target_file), 'a') as f_p:
                    f_p.write(message + "\n")
        except Exception as e:
            print(f"Error writing log: {e}")

    def safe_mpf(value, name):
        """ Convert a value to mpmath.mpf safely, logging any issues. """
        try:
            return mpmath.mpf(value)
        except Exception as e:
            print(f"ERROR: Failed to convert {name}={value} to mpf: {e}")
            return ZERO  # Default to zero to avoid crashes

    def log_all(iteration, rmax, phi_var, phi_try, phi_end, rel_diff, file):

        # Ensure all values are valid before using them
        phi_var = safe_mpf(phi_var, "phi_var")
        phi_try = safe_mpf(phi_try, "phi_try")
        phi_end = safe_mpf(phi_end, "phi_end")
        rel_diff = safe_mpf(rel_diff, "rel_diff")

        log_message(
            f"{iteration} \t rmax={mpmath.nstr(mpmath.fdiv(mpmath.mpf(rmax), mpmath.mpf(km)), 5)} km \t "
            f"phi_var={mpmath.nstr(mpmath.fdiv(mpmath.mpf(phi_var), mpmath.mpf(Mp)), 50)} Mp \t "
            f"phi_try={mpmath.nstr(mpmath.fdiv(mpmath.mpf(phi_try), mpmath.mpf(Mp)), 50)} Mp \t "
            f"phi_end={mpmath.nstr(mpmath.fdiv(mpmath.mpf(phi_end), mpmath.mpf(Mp)), 50)} Mp \t "
            f"rel_diff={mpmath.nstr(rel_diff, 50)}",
            file)

    def restart_shooting(problem_type, current_rel_diff):
        nonlocal rmax, phi_try, phi_var, i, stagnation_count, slow_convergence_count, restart_count, best_attempt_id, best_attempt_count
        restart_count += 1

        if problem_type in ["initial_negative", "initial_overbackground"]:
            phi_try = mpmath.fadd(
                phi_try, phi_var) if problem_type == "initial_negative" else mpmath.fsub(phi_try, phi_var)
            # Reduce step size to prevent looping
            phi_var = mpmath.fmul(phi_var, mpmath.mpf(1e-1))

            log_message(
                f"Adjusted for {problem_type}: phi_try={mpmath.nstr(mpmath.fdiv(phi_try, Mp), 50)} Mp, phi_var={mpmath.nstr(mpmath.fdiv(phi_var, Mp), 50)} Mp.", file)

        elif problem_type in ["negative_values", "unphysical_growth", "stagnation", "slow_convergence"]:
            if best_attempt is not None:
                if best_attempt_id == id(best_attempt):
                    best_attempt_count += 1
                else:
                    best_attempt_count = 0
                    best_attempt_id = id(best_attempt)
                log_message(
                    f"Retrieving best attempt: rmax={mpmath.nstr(mpmath.fdiv(best_attempt[0], km), 5)} km, phi_var={mpmath.nstr(mpmath.fdiv(best_attempt[1], Mp), 50)} Mp, phi_try={mpmath.nstr(mpmath.fdiv(best_attempt[2], Mp), 50)} Mp, rel_diff={mpmath.nstr(best_attempt[3], 50)}.", file)
                rmax, phi_var, phi_try, rel_diff = best_attempt  # Restore best conditions

            if best_attempt_count >= 2:
                log_message(
                    "Best attempt retrieved twice, increasing adjustments.", file)
                phi_var = mpmath.fmul(phi_var, TWO)
                rmax = mpmath.fmul(rmax, mpmath.mpf(1.1))
                phi_try = mpmath.fadd(phi_try, mpmath.fmul(phi_var, mpmath.mpf(
                    2))) if rel_diff > 0 else mpmath.fsub(phi_try, mpmath.fmul(phi_var, TWO))
            else:
                phi_var = mpmath.fmul(phi_var, mpmath.mpf(1.1))
                rmax = mpmath.fmul(rmax, mpmath.mpf(1.01))
                phi_try = mpmath.fadd(
                    phi_try, phi_var) if current_rel_diff > 0 else mpmath.fsub(phi_try, phi_var)

            if best_attempt_count > 3:
                log_message(
                    "Repeated best attempt retrieval, performing soft reset.", file)
                # Move toward expected solution
                phi_try = mpmath.fmul(HALF, mpmath.fadd(phi_try, phi_inf))
                phi_var = mpmath.fmul(phi_var, TWO)

        # Make sure phi_try is positive
        phi_try = mpmath.fabs(phi_try)

        stagnation_count = 0  # Reset stagnation counter
        slow_convergence_count = 0  # Reset convergence counter
        log_message(f"Adjusted parameters: rmax={mpmath.nstr(mpmath.fdiv(mpmath.mpf(rmax), mpmath.mpf(km)), 5)} km, phi_var={mpmath.nstr(mpmath.fdiv(mpmath.mpf(phi_var), mpmath.mpf(Mp)), 50)} Mp, phi_try={mpmath.nstr(mpmath.fdiv(mpmath.mpf(phi_try), mpmath.mpf(Mp)), 50)} Mp.", file)

    i = 0  # Initialize iteration counter
    phi_try = mpmath.mpf(phi0)  # Start with initial phi value
    """
    if phi_try < phiv:
        phi_var = mpmath.fmul(mpmath.mpf(0.1), phi_try)
    else:
    """
    phi_var = mpmath.mpf(phiv)  # Start with original value
    prev_diff_sign = None  # Track difference sign for adjustments
    last_phi_end = None  # Track last phi_end for stagnation
    last_rel_diff = None  # Track last relative difference
    stagnation_count = 0  # Stagnation counter
    slow_convergence_count = 0  # Counter for slow convergence
    restart_count = 0  # Count restarts due to non-convergence
    restart_limit = int(1e2)
    best_attempt = None
    best_attempt_id = None
    best_attempt_count = 0
    best_rel_diff = mpmath.mpf('inf')

    log_message(
        f"The scalar field minimum at infinity is phi_inf={mpmath.nstr(mpmath.fdiv(mpmath.mpf(phi_inf), mpmath.mpf(Mp)), 50)} Mp\n", file)

    while i < max_itr:
        # Limit the maximum number of restarts to avoid infinite loops
        if restart_count > restart_limit:
            log_message("Restart limit reached. Exiting...", file)
            if best_attempt is not None:
                log_message(
                    f"Best attempt conditions: rmax={mpmath.nstr(mpmath.fdiv(best_attempt[0], km), 5)} km, phi_var={mpmath.nstr(mpmath.fdiv(best_attempt[1], Mp), 50)} Mp, phi_try={mpmath.nstr(mpmath.fdiv(best_attempt[2], Mp), 50)} Mp, rel_diff={mpmath.nstr(best_attempt[3], 50)}.", file)
            return None

        if phi_try < 0 and scalar != self.dilaton:
            log_message(
                f"The initial scalar field value became negative.", file)
            restart_shooting("initial_negative", ZERO)
            continue

        if phi_try > mpmath.fmul(mpmath.mpf(1e1), Mp):
            log_message(
                f"The initial scalar field went beyond the Planck mass.", file)
            restart_shooting("initial_overplanckian", ZERO)
            continue

        # Integrate the ODE system
        prf = self.scalar_ode_sol(
            ode_sys, scalar, scalar_par, dens_c, rmax, dr, phi_try, sigma0, atm, EoS, gamma)

        # Check if scalar integration was valid
        if not prf:
            log_message(
                "Unphysical growth of the scalar field during run.", file)
            restart_shooting("unphysical_growth", ZERO)
            continue

        r_prof, p_prof, m_prof, f_prof, s_prof, neg_phi, small_r, large_r = prf

        if neg_phi > 0:
            # There are four steps in the RK4 algorithm
            # Ensure steps is an integer
            steps = 4 * (len(r_prof) - 1)
            # Format percentage to 5 decimal places
            percentage = mpmath.nstr(mpmath.fmul(
                mpmath.fdiv(neg_phi, steps), mpmath.mpf(1e2)), 5)
            log_message(
                f"The scalar field has become negative {neg_phi} times out of {steps} ({percentage} %) steps between r = {mpmath.nstr(mpmath.fdiv(small_r, km), 5)} km and r = {mpmath.nstr(mpmath.fdiv(large_r, km), 5)} km, times in which we switched its sign.", file)

        # Log phi_end
        phi_end = f_prof[-1]  # Last scalar field value
        rel_diff = mpmath.fdiv(mpmath.fsub(phi_inf, phi_end), phi_inf)
        log_all(i, rmax, phi_var, phi_try, phi_end, rel_diff, file)

        if mpmath.isfinite(rel_diff) and mpmath.fabs(rel_diff) < mpmath.fabs(best_rel_diff) and neg_phi == 0:
            best_rel_diff = rel_diff
            best_attempt = (rmax, phi_var, phi_try, rel_diff)

        if log_dir is not None:
            # Unique file for each attempt
            attempt_file = f"density{float(dens_c)/float(g_cm3_to_code_units):.2e}_attempt_{i}.txt"
            with open(os.path.join(log_dir, attempt_file), 'w') as a_f:
                a_f.write("# r\tP\tm\tf\ts\n")
                for r, P, M, f, s in zip(r_prof, p_prof, m_prof, f_prof, s_prof):
                    a_f.write(
                        f"{mpmath.nstr(r, 10)}\t{mpmath.nstr(P, 10)}\t{mpmath.nstr(M, 10)}\t{mpmath.nstr(f, 10)}\t{mpmath.nstr(s, 10)}\n")

        # Detect slow convergence
        if last_rel_diff is not None and mpmath.fabs(rel_diff) >= mpmath.fabs(last_rel_diff):
            slow_convergence_count += 1
        else:
            slow_convergence_count = 0
        last_rel_diff = rel_diff

        if slow_convergence_count > 5:
            log_message(
                "Slow convergence detected: relative difference not decreasing for several iterations.", file)
            restart_shooting("slow_convergence", rel_diff)
            continue

        # Detect stagnation
        if last_phi_end is not None and mpmath.fabs(mpmath.fdiv(mpmath.fsub(phi_end, last_phi_end), phi_end)) < phi_tol:
            stagnation_count += 1
        else:
            stagnation_count = 0
        last_phi_end = phi_end

        if stagnation_count > 5:
            log_message(
                "Stagnation detected: phi_end stuck for several iterations.", file)
            restart_shooting("stagnation", rel_diff)
            continue

        # Check for values close to phi_inf
        close_indices = [idx for idx, f in enumerate(f_prof) if mpmath.fabs(
            mpmath.fdiv(mpmath.fsub(phi_inf, f), phi_inf)) < phi_tol]
        if close_indices:
            last_idx = close_indices[-1]
            r_last = r_prof[last_idx]
            f_last = f_prof[last_idx]

            # Compute stellar radius
            if crit == 0:
                ind = next(
                    (i for i in range(1, len(m_prof))
                     if m_prof[i] - m_prof[i-1] < m_rel),
                    len(m_prof) - 1
                )
            elif crit == 1:
                ind = next(
                    (i for i in range(len(p_prof))
                     if p_prof[i] < p_tol),
                    len(p_prof) - 1
                )
                if p_prof[ind] >= p_tol:
                    log_message(
                        f"Pressure never goes below the threshold. Scalar field near target: phi={mpmath.nstr(mpmath.fdiv(f_last, Mp), 50)} Mp at r={mpmath.nstr(mpmath.fdiv(r_last, km), 5)} km.", file)
                    continue
            else:
                raise ValueError("No stellar radius criterion selected.")

            if ind == 0 or ind >= len(m_prof):
                log_message(
                    "Invalid index for stellar radius criterion.", file)
                continue

            stellar_radius = r_prof[ind - 1]
            # Ensure neg_phi == 0 before logging success
            if mpmath.fmul(stellar_radius, THREE) < r_last and neg_phi == 0:
                log_message(
                    f"Success: Scalar field close at r={mpmath.nstr(mpmath.fdiv(r_last, km), 5)} km, phi={mpmath.nstr(mpmath.fdiv(f_last, Mp), 50)} Mp. Stellar radius={mpmath.nstr(mpmath.fdiv(stellar_radius, km), 5)} km.", file)
                return (r_prof[:last_idx+1], p_prof[:last_idx+1], m_prof[:last_idx+1], f_prof[:last_idx+1], s_prof[:last_idx+1])

        # Adjust phi_try
        if m == 0:
            phi_try = mpmath.fadd(phi_try, phi_var)
        elif m == 1:
            diff_sign = 1 if rel_diff >= 0 else -1
            if prev_diff_sign is not None and diff_sign != prev_diff_sign:
                # Reduce step size
                phi_var = mpmath.fmul(phi_var, mpmath.mpf(1e-1))
            if rel_diff >= 0:
                phi_try = mpmath.fadd(phi_try, phi_var)
            else:
                phi_try = mpmath.fsub(phi_try, phi_var)
            prev_diff_sign = diff_sign
        else:
            raise ValueError("No convergence method selected.")

        i += 1  # Increment iteration

    log_message(
        f"Maximum iterations reached: {max_itr} for dens_c={mpmath.nstr(mpmath.fdiv(dens_c, g_cm3_to_code_units), 2)} g/cm3", file)

    if best_attempt is not None:
        log_message(
            f"Best attempt conditions: rmax={mpmath.nstr(mpmath.fdiv(best_attempt[0], km), 5)} km, phi_var={mpmath.nstr(mpmath.fdiv(best_attempt[1], Mp), 50)} Mp, phi_try={mpmath.nstr(mpmath.fdiv(best_attempt[2], Mp), 50)} Mp, rel_diff={mpmath.nstr(best_attempt[3], 50)}.", file)

    return None