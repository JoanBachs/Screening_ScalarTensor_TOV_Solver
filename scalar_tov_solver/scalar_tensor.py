# scalar_tensor.py
from .core import *

"""
Systems of Ordinary Differential Equations
"""


def scalar_nwt_ode(self, r, P, m, f, s, scalar, scalar_par, atm_den, EoS_param):
    scalar_functions = scalar(f, scalar_par)
    if not scalar_functions:
        return scalar_functions
    A_f, V_f, dA_df, dV_df, negative_field = scalar_functions
    # Quintessence or Phantom parameter, set to 1 by default
    QP = ONE
    if scalar == self.inverse_symmetron:
        QP = mpmath.fneg(QP)

    eden = self.E_of_P_general(P, EoS_param)
    if eden < atm_den:
        dPdr = ZERO
        dmdr = ZERO

    else:
        dPdr = mpmath.fneg(mpmath.fmul(eden, mpmath.fadd(mpmath.fdiv(
            m, mpmath.fmul(r, r)), mpmath.fmul(mpmath.fdiv(dA_df, A_f), s))))
        dmdr = mpmath.fprod(
            [HALF, kappa2, r, r, mpmath.power(A_f, FOUR), eden])

    dfdr = s
    Veff = mpmath.fadd(dV_df, mpmath.fprod(
        [dA_df, mpmath.power(A_f, THREE), eden]))
    dsdr = mpmath.fadd(mpmath.fneg(mpmath.fdiv(
        mpmath.fmul(TWO, s), r)), mpmath.fmul(QP, Veff))

    return [mpmath.re(dPdr), mpmath.re(dmdr), mpmath.re(dfdr), mpmath.re(dsdr), negative_field]


def scalar_tov_ode(self, r, P, m, f, s, scalar, scalar_par, atm_den, EoS_param):
    scalar_functions = scalar(f, scalar_par)
    if not scalar_functions:
        return scalar_functions
    A_f, V_f, dA_df, dV_df, negative_field = scalar_functions

    den = mpmath.fsub(r, mpmath.fmul(TWO, m))
    if den <= ZERO:
        return False

    eden = self.E_of_P_general(P, EoS_param)
    dndr = mpmath.fadd(mpmath.fdiv(mpmath.fadd(mpmath.fdiv(m, r), mpmath.fprod([HALF, kappa2, r, r, mpmath.fsub(mpmath.fmul(
        mpmath.power(A_f, FOUR), P), V_f)])), den), mpmath.fprod([HALF, HALF, kappa2, r, s, s]))

    if eden < atm_den:
        # If we are outside the star, P and m stop evolving
        dPdr = dmdr = ZERO
    else:
        dPdr = mpmath.fneg(mpmath.fmul(mpmath.fadd(eden, P), mpmath.fadd(
            dndr, mpmath.fmul(mpmath.fdiv(dA_df, A_f), s))))
        dmdr = mpmath.fprod([HALF, kappa2, r, r, mpmath.fsum([mpmath.fmul(mpmath.power(A_f, FOUR), eden), mpmath.fprod(
            [HALF, mpmath.fsub(ONE, mpmath.fdiv(mpmath.fmul(TWO, m), r)), s, s]), V_f])])

    dfdr = s
    Veff = mpmath.fsub(dV_df, mpmath.fprod([dA_df, mpmath.power(
        A_f, THREE), mpmath.fsub(mpmath.fmul(THREE, P), eden)]))
    dsdr = mpmath.fadd(mpmath.fmul(mpmath.fsub(mpmath.fdiv(mpmath.fsum([dmdr, mpmath.fmul(THREE, mpmath.fdiv(m, r)), mpmath.fneg(
        TWO)]), den), dndr), s), mpmath.fdiv(mpmath.fmul(r, Veff), den))

    return [mpmath.re(dPdr), mpmath.re(dmdr), mpmath.re(dfdr), mpmath.re(dsdr), negative_field]



"""
Functions to Integrate the ODE Systems & Calculate Stellar Mass and Radius
"""

def scalar_ode_sol(self, ode_sys, scalar, scalar_par, c_den, rmax, dr, f, s, atm_den, EoS_param):
    h = mpmath.mpf(dr)  # Initial step size
    n_steps = int(mpmath.ceil(mpmath.fdiv(rmax, dr)))
    r = [mpmath.fmul(h, mpmath.mpf(i)) for i in range(1, n_steps + 1)]

    neg_phi = 0
    smallest_radius = mpmath.mpf('inf')  # Start with infinity for smallest
    largest_radius = ZERO  # Start with zero for largest

    eden = self.E_of_rho_general(c_den, EoS_param)
    P = self.P_of_rho_general(c_den, EoS_param)
    if scalar == self.chameleon:
        P_atm = self.P_of_rho_general(atm_den, EoS_param)

    scalar_functions = scalar(f, scalar_par)
    if not scalar_functions:
        return scalar_functions

    A_f, V_f, dA_df, dV_df, neg_phi_0 = scalar_functions
    m = mpmath.fdiv(mpmath.fprod([HALF, kappa2, mpmath.power(r[0], THREE), mpmath.fadd(
        mpmath.fmul(mpmath.power(A_f, FOUR), eden), V_f)]), THREE)

    p_R, m_R, f_R, s_R = [
        ZERO] * len(r), [ZERO] * len(r), [ZERO] * len(r), [ZERO] * len(r)
    p_R[0], m_R[0], f_R[0], s_R[0] = P, m, f, s

    def update_neg_field_tracking(neg_phi, radius, neg_flag, smallest, largest):
        if neg_flag:
            neg_phi += 1
            smallest = min(smallest, radius)
            largest = max(largest, radius)
        return neg_phi, smallest, largest

    # Track the first step outside the loop
    neg_phi, smallest_radius, largest_radius = update_neg_field_tracking(
        neg_phi, r[0], neg_phi_0, smallest_radius, largest_radius)

    for i in range(1, len(r)):
        profiles_1 = ode_sys(r[i - 1], p_R[i - 1], m_R[i - 1], f_R[i - 1],
                             s_R[i - 1], scalar, scalar_par, atm_den, EoS_param)
        if profiles_1 is False or profiles_1 is None:
            return profiles_1
        p1, m1, f1, s1, neg_phi_1 = profiles_1
        k1_p, k1_m, k1_f, k1_s = mpmath.fmul(h, p1), mpmath.fmul(
            h, m1), mpmath.fmul(h, f1), mpmath.fmul(h, s1)

        profiles_2 = ode_sys(mpmath.fadd(r[i - 1], mpmath.fmul(HALF, h)), mpmath.fadd(p_R[i - 1], mpmath.fmul(HALF, k1_p)), mpmath.fadd(m_R[i - 1], mpmath.fmul(
            HALF, k1_m)), mpmath.fadd(f_R[i - 1], mpmath.fmul(HALF, k1_f)), mpmath.fadd(s_R[i - 1], mpmath.fmul(HALF, k1_s)), scalar, scalar_par, atm_den, EoS_param)
        if profiles_2 is False or profiles_2 is None:
            return profiles_2
        p2, m2, f2, s2, neg_phi_2 = profiles_2
        k2_p, k2_m, k2_f, k2_s = mpmath.fmul(h, p2), mpmath.fmul(
            h, m2), mpmath.fmul(h, f2), mpmath.fmul(h, s2)

        profiles_3 = ode_sys(mpmath.fadd(r[i - 1], mpmath.fmul(HALF, h)), mpmath.fadd(p_R[i - 1], mpmath.fmul(HALF, k2_p)), mpmath.fadd(m_R[i - 1], mpmath.fmul(
            HALF, k2_m)), mpmath.fadd(f_R[i - 1], mpmath.fmul(HALF, k2_f)), mpmath.fadd(s_R[i - 1], mpmath.fmul(HALF, k2_s)), scalar, scalar_par, atm_den, EoS_param)
        if profiles_3 is False or profiles_3 is None:
            return profiles_3
        p3, m3, f3, s3, neg_phi_3 = profiles_3
        k3_p, k3_m, k3_f, k3_s = mpmath.fmul(h, p3), mpmath.fmul(
            h, m3), mpmath.fmul(h, f3), mpmath.fmul(h, s3)

        profiles_4 = ode_sys(mpmath.fadd(r[i - 1], h), mpmath.fadd(p_R[i - 1], k3_p), mpmath.fadd(m_R[i - 1], k3_m), mpmath.fadd(
            f_R[i - 1], k3_f), mpmath.fadd(s_R[i - 1], k3_s), scalar, scalar_par, atm_den, EoS_param)
        if profiles_4 is False or profiles_4 is None:
            return profiles_4
        p4, m4, f4, s4, neg_phi_4 = profiles_4
        k4_p, k4_m, k4_f, k4_s = mpmath.fmul(h, p4), mpmath.fmul(
            h, m4), mpmath.fmul(h, f4), mpmath.fmul(h, s4)

        p_R[i] = mpmath.fadd(p_R[i - 1], mpmath.fdiv(mpmath.fsum([k1_p,
                             mpmath.fmul(TWO, k2_p), mpmath.fmul(TWO, k3_p), k4_p]), SIX))
        m_R[i] = mpmath.fadd(m_R[i - 1], mpmath.fdiv(mpmath.fsum([k1_m,
                             mpmath.fmul(TWO, k2_m), mpmath.fmul(TWO, k3_m), k4_m]), SIX))
        f_R[i] = mpmath.fadd(f_R[i - 1], mpmath.fdiv(mpmath.fsum([k1_f,
                             mpmath.fmul(TWO, k2_f), mpmath.fmul(TWO, k3_f), k4_f]), SIX))
        s_R[i] = mpmath.fadd(s_R[i - 1], mpmath.fdiv(mpmath.fsum([k1_s,
                             mpmath.fmul(TWO, k2_s), mpmath.fmul(TWO, k3_s), k4_s]), SIX))

        # Track negative fields
        neg_phi, smallest_radius, largest_radius = update_neg_field_tracking(
            neg_phi, r[i - 1], neg_phi_1, smallest_radius, largest_radius)
        neg_phi, smallest_radius, largest_radius = update_neg_field_tracking(
            neg_phi, r[i - 1], neg_phi_2, smallest_radius, largest_radius)
        neg_phi, smallest_radius, largest_radius = update_neg_field_tracking(
            neg_phi, r[i - 1], neg_phi_3, smallest_radius, largest_radius)
        neg_phi, smallest_radius, largest_radius = update_neg_field_tracking(
            neg_phi, r[i - 1], neg_phi_4, smallest_radius, largest_radius)

        if scalar == self.chameleon and p_R[i] < P_atm:
            p_R[i] = mpmath.mpf(P_atm)

    return (r, p_R, m_R, f_R, s_R, neg_phi, smallest_radius, largest_radius)


"""
Algorithms for Saving Stellar Masses and Radii & Radial Profiles
"""

def scalar_star_MR(self, ode_sys, scalar, scalar_par, c_den, rmax, dr, atm, EoS, m_rel, p_tol, crit, s_0, f_inf, f_tol, f_var, f_ini, m_itr, mtd, log, out, attdir):
    data = []
    data.append(
        f"STTsolverRK4: atm_den(g/cm3)={float(atm)/float(g_cm3_to_code_units):.2e}\tm_rel={float(m_rel):.1e}\tp_tol(dyn/cm2)={float(p_tol)/float(dyn_cm2_to_code_units):.2e}\tcrit={crit}\n")

    def safe_div(numerator, denominator):
        if denominator == 0:
            return mpmath.fsub(numerator, denominator), "difference"
        return mpmath.fdiv(numerator, denominator), "ratio"

    for den in c_den:

        if f_ini is None:
            f_est = self.phi_ini(scalar, scalar_par, mtd,
                                 f_var, f_inf, den, EoS, log)
        else:
            f_est = mpmath.fmul(mpmath.mpf(f_ini), Mp)
        prf = self.shoot_mtd(ode_sys, scalar, scalar_par, mtd, m_itr, f_tol, f_var,
                             den, f_est, f_inf, s_0, rmax, dr, atm, EoS, m_rel, p_tol, crit, log, attdir)

        # Use the appropriate solver function based on eos_params
        if prf is not None:
            r_prf, P_prf, M_prf, f_prf, s_prf = prf
            # Target-to-calculated scalar field ratio at infinity
            print(
                f"Target-to-calculated scalar field ratio at infinity: {mpmath.nstr(mpmath.fdiv(mpmath.mpf(f_inf), mpmath.mpf(f_prf[-1])), 20)}")
            # Scalar field initial value variation
            ratio_or_diff_value, ratio_or_diff_label = safe_div(
                mpmath.mpf(f_prf[0]), mpmath.mpf(f_est))
            print(
                f"Scalar field initial value variation: phi0_ini = {mpmath.nstr(mpmath.fdiv(mpmath.mpf(f_est), mpmath.mpf(Mp)), 20)} Mp, phi0_fin = {mpmath.nstr(mpmath.fdiv(mpmath.mpf(f_prf[0]), mpmath.mpf(Mp)), 20)} Mp, {ratio_or_diff_label} = {mpmath.nstr(ratio_or_diff_value, 20)}")
            # Maximum integration distance variation
            print(
                f"Maximum integration distance variation: rmax_ini = {mpmath.nstr(mpmath.fdiv(mpmath.mpf(rmax), mpmath.mpf(km)), 5)} km, rmax_fin = {mpmath.nstr(mpmath.fdiv(mpmath.mpf(r_prf[-1]), mpmath.mpf(km)), 5)} km, ratio = {mpmath.nstr(mpmath.fdiv(mpmath.mpf(r_prf[-1]), mpmath.mpf(rmax)), 5)}")

            if crit == 0:
                ind = next(
                    (i for i in range(1, len(M_prf))
                     if M_prf[i] - M_prf[i-1] < m_rel),
                    len(M_prf) - 1
                )
            elif crit == 1:
                ind = next(
                    (i for i in range(len(P_prf))
                     if P_prf[i] < p_tol),
                    len(P_prf) - 1
                )
                if P_prf[ind] >= p_tol:
                    print("The pressure never goes below the specified threshold.")
            else:
                raise ValueError("No stellar radius criterion was selected.")
            if ind == 0 or ind >= len(M_prf):
                print("Invalid index detected for stellar radius criterion.")
                continue
            M = mpmath.fdiv(M_prf[ind - 1], mpmath.fmul(Grav_C, Msun))
            R = mpmath.fdiv(r_prf[ind - 1], km)
            data.append(
                f"{mpmath.nstr(mpmath.fdiv(mpmath.mpf(den), mpmath.mpf(g_cm3_to_code_units)), 10)}\t{mpmath.nstr(M, 10)}\t{mpmath.nstr(R, 10)}\n")

            # Create scalar field profile in Mp units
            phi_prf = [mpmath.fdiv(f, Mp) for f in f_prf]

            # Create scalar gradient * Radius profile in Mp units
            sxR_prf = [mpmath.fdiv(mpmath.fprod([s, R, km]), Mp)
                       for s in s_prf]

            # Write sample profile to file
            data.extend([f"{mpmath.nstr(r, 10)}\t{mpmath.nstr(P, 10)}\t{mpmath.nstr(M, 10)}\t{mpmath.nstr(phi, 10)}\t{mpmath.nstr(s, 10)}\t{mpmath.nstr(sxR, 10)}\n" for r,
                        P, M, phi, s, sxR in zip(r_prf, P_prf, M_prf, phi_prf, s_prf, sxR_prf)])

        else:
            # Log if shooting method failed for this central density
            log.write(
                f"Shooting method failed for c_den={mpmath.nstr(mpmath.fdiv(mpmath.mpf(den), mpmath.mpf(g_cm3_to_code_units)), 10)} g/cm3\n")

    # Write entire lists to the respective files
    out.writelines(data)
    out.flush()

    return None


"""
Scalar Fields with Screening Mechanisms
"""

def chameleon(self, field, parameters):
    """
    The exponent is kept below a certain limit
    The field is kept positive
    """
    n, a, L = parameters
    negative_field = False

    if field < 0:
        field = mpmath.fneg(field)
        negative_field = True

    exponent = mpmath.fmul(a, field)
    if exponent > mpmath.mpf('1e2'):
        return False

    A_f = mpmath.exp(exponent)
    dA_df = mpmath.fmul(a, A_f)
    d2A_df2 = mpmath.fmul(a, dA_df)

    V_f = mpmath.fdiv(mpmath.power(L, mpmath.fadd(n, FOUR)),
                      mpmath.power(field, n))
    dV_df = mpmath.fneg(mpmath.fmul(n, mpmath.fdiv(V_f, field)))
    d2V_df2 = mpmath.fneg(mpmath.fmul(
        mpmath.fadd(n, 1), mpmath.fdiv(dV_df, field)))

    return A_f, V_f, dA_df, dV_df, negative_field


def symmetron(self, field, parameters):
    """
    The field is kept positive
    """
    mu, lmbd, Mcpl = parameters
    negative_field = False

    if field < 0:
        field = mpmath.fneg(field)
        negative_field = True

    A_f = mpmath.fadd(ONE, mpmath.fprod(
        [HALF, mpmath.fdiv(field, Mcpl), mpmath.fdiv(field, Mcpl)]))
    V_f = mpmath.fadd(mpmath.fneg(mpmath.fprod([HALF, mu, mu, field, field])), mpmath.fprod(
        [mpmath.mpf('0.25'), lmbd, mpmath.power(field, FOUR)]))
    dA_df = mpmath.fdiv(field, mpmath.fmul(Mcpl, Mcpl))
    dV_df = mpmath.fadd(mpmath.fneg(mpmath.fprod([mu, mu, field])), mpmath.fmul(
        lmbd, mpmath.power(field, THREE)))

    return A_f, V_f, dA_df, dV_df, negative_field


def dilaton(self, field, parameters):
    """
    The exponent is kept below a certain limit
    """
    Vo, a2, fd = parameters
    negative_field = False

    exponent = mpmath.fdiv(mpmath.fsub(fd, field), Mp)
    if exponent > mpmath.mpf('1e2'):
        return False

    A_f = mpmath.fadd(ONE, mpmath.fprod(
        [HALF, a2, mpmath.power(exponent, TWO)]))
    V_f = mpmath.fprod([mpmath.power(A_f, FOUR), Vo, mpmath.exp(exponent)])
    dA_df = mpmath.fmul(a2, mpmath.fdiv(
        mpmath.fsub(field, fd), mpmath.fmul(Mp, Mp)))
    dV_df = mpmath.fsub(mpmath.fprod([FOUR, mpmath.power(A_f, THREE), dA_df, Vo, mpmath.exp(
        exponent)]), mpmath.fdiv(mpmath.fprod([mpmath.power(A_f, FOUR), Vo, mpmath.exp(exponent)]), Mp))

    return A_f, V_f, dA_df, dV_df, negative_field