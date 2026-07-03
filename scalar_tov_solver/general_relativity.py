# general_relativity.py
from .core import *

"""
Systems of Ordinary Differential Equations
"""


def nwt_ode(self, r, P, m, atm_den, EoS_param):
    eden = self.E_of_P_general(P, EoS_param)
    if eden < atm_den:
        # If we are outside the star, P and m stop evolving
        return [ZERO, ZERO]

    dFdr = mpmath.fdiv(m, mpmath.power(r, TWO))
    dPdr = mpmath.fneg(mpmath.fmul(eden, dFdr))
    dmdr = mpmath.fprod([HALF, kappa2, r, r, eden])

    return [mpmath.re(dPdr), mpmath.re(dmdr)]


def tov_ode(self, r, P, m, atm_den, EoS_param):
    eden = self.E_of_P_general(P, EoS_param)
    if eden < atm_den:
        # If we are outside the star, P and m stop evolving
        return [ZERO, ZERO]
    
    den = mpmath.fsub(r, mpmath.fmul(TWO, m))
    if den <= ZERO:
        return [ZERO, ZERO]

    dndr = mpmath.fdiv(mpmath.fadd(mpmath.fdiv(m, r), mpmath.fprod(
        [HALF, kappa2, r, r, P])), den)
    dPdr = mpmath.fneg(mpmath.fmul(mpmath.fadd(eden, P), dndr))
    dmdr = mpmath.fprod([HALF, kappa2, r, r, eden])

    return [mpmath.re(dPdr), mpmath.re(dmdr)]


"""
Functions to Integrate the ODE Systems & Calculate Stellar Mass and Radius
"""


def ode_sol(self, ode_sys, c_den, rmax, dr, atm_den, EoS_param, m_rel, p_tol, crit):
    h = mpmath.mpf(dr)
    n_steps = int(mpmath.ceil(mpmath.fdiv(rmax, dr)))
    r = [mpmath.fmul(h, mpmath.mpf(i)) for i in range(1, n_steps + 1)]

    eden = self.E_of_rho_general(c_den, EoS_param)
    P = self.P_of_rho_general(c_den, EoS_param)
    P_atm = self.P_of_rho_general(atm_den, EoS_param)

    p_R = [ZERO] * len(r)
    m_R = [ZERO] * len(r)

    p_R[0] = mpmath.mpf(P)
    m_R[0] = mpmath.fdiv(mpmath.fprod(
        [HALF, kappa2, mpmath.power(r[0], THREE), eden]), THREE)

    for i in range(1, len(r)):
        p1, m1 = ode_sys(r[i - 1], p_R[i - 1], m_R[i - 1], atm_den, EoS_param)
        k1_p, k1_m = mpmath.fmul(h, p1), mpmath.fmul(h, m1)

        p2, m2 = ode_sys(mpmath.fadd(r[i - 1], mpmath.fmul(HALF, h)), mpmath.fadd(p_R[i - 1], mpmath.fmul(
            HALF, k1_p)), mpmath.fadd(m_R[i - 1], mpmath.fmul(HALF, k1_m)), atm_den, EoS_param)
        k2_p, k2_m = mpmath.fmul(h, p2), mpmath.fmul(h, m2)

        p3, m3 = ode_sys(mpmath.fadd(r[i - 1], mpmath.fmul(HALF, h)), mpmath.fadd(p_R[i - 1], mpmath.fmul(
            HALF, k2_p)), mpmath.fadd(m_R[i - 1], mpmath.fmul(HALF, k2_m)), atm_den, EoS_param)
        k3_p, k3_m = mpmath.fmul(h, p3), mpmath.fmul(h, m3)

        p4, m4 = ode_sys(mpmath.fadd(r[i - 1], h), mpmath.fadd(
            p_R[i - 1], k3_p), mpmath.fadd(m_R[i - 1], k3_m), atm_den, EoS_param)
        k4_p, k4_m = mpmath.fmul(h, p4), mpmath.fmul(h, m4)

        p_R[i] = mpmath.fadd(p_R[i - 1], mpmath.fdiv(mpmath.fsum([k1_p,
                             mpmath.fmul(TWO, k2_p), mpmath.fmul(TWO, k3_p), k4_p]), SIX))
        m_R[i] = mpmath.fadd(m_R[i - 1], mpmath.fdiv(mpmath.fsum([k1_m,
                             mpmath.fmul(TWO, k2_m), mpmath.fmul(TWO, k3_m), k4_m]), SIX))

        if p_R[i] < P_atm:
            p_R[i] = mpmath.mpf(P_atm)

    if crit == 0:
        ind = next(
            (i for i in range(1, len(m_R))
             if m_R[i] - m_R[i-1] < m_rel),
            len(m_R) - 1
        )
    elif crit == 1:
        ind = next(
            (i for i in range(len(p_R))
             if p_R[i] < p_tol),
            len(p_R) - 1
        )
        if p_R[ind] >= p_tol:
            print("The pressure never goes below the specified threshold.")
    else:
        raise ValueError("No stellar radius criterion was selected.")

    M = m_R[ind - 1]
    R = r[ind - 1]

    return mpmath.fdiv(R, km), mpmath.fdiv(M, mpmath.fmul(Grav_C, Msun)), (r, p_R, m_R)



"""
Algorithms for Saving Stellar Masses and Radii & Radial Profiles
"""


def star_MR(self, ode_sys, c_den, rmax, dr, atm_den, EoS_param, m_rel, p_tol, crit, out):
    digits = 10
    data = []
    data.append(
        f"RK4-MRcurves: atm_den={mpmath.nstr(atm_den, digits)}\tm_rel={mpmath.nstr(m_rel, digits)}\tp_tol={mpmath.nstr(p_tol, digits)}\tcrit={crit}\n")

    for den in c_den:
        R, M, prf = self.ode_sol(
            ode_sys, den, rmax, dr, atm_den, EoS_param, m_rel, p_tol, crit)
        data.append(
            f"{mpmath.nstr(mpmath.fdiv(mpmath.mpf(den), mpmath.mpf(g_cm3_to_code_units)), 10)}\t{mpmath.nstr(M, digits)}\t{mpmath.nstr(R, digits)}\n")
        r_prf, P_prf, M_prf = prf
        data.extend([f"{mpmath.nstr(r, digits)}\t{mpmath.nstr(P, digits)}\t{mpmath.nstr(M, digits)}\n" for r, P, M in zip(
            r_prf, P_prf, M_prf)])

    out.writelines(data)

    return None


