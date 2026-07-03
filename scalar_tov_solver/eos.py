# eos.py
from .core import *

"""
4-Piecewise Polytropic Equation of State
"""


def PiecePolyEoS(self, P1, Gamma1, Gamma2, Gamma3, rho1, rho2, rho0, Gamma0):
    K1 = mpmath.fdiv(P1, mpmath.power(rho1, Gamma1))
    K2 = mpmath.fmul(K1, mpmath.power(rho1, mpmath.fsub(Gamma1, Gamma2)))
    K3 = mpmath.fmul(K2, mpmath.power(rho2, mpmath.fsub(Gamma2, Gamma3)))

    P0 = mpmath.fmul(K1, mpmath.power(rho0, Gamma1))
    K0 = mpmath.fdiv(P0, mpmath.power(rho0, Gamma0))
    P2 = mpmath.fmul(K3, mpmath.power(rho2, Gamma3))

    a0 = ZERO
    E0 = mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a0), rho0), mpmath.fmul(
        mpmath.fdiv(K0, mpmath.fsub(Gamma0, ONE)), mpmath.power(rho0, Gamma0)))
    a1 = mpmath.fsub(mpmath.fsub(mpmath.fdiv(E0, rho0), ONE), mpmath.fmul(mpmath.fdiv(
        K1, mpmath.fsub(Gamma1, ONE)), mpmath.power(rho0, mpmath.fsub(Gamma1, ONE))))
    E1 = mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a1), rho1), mpmath.fmul(
        mpmath.fdiv(K1, mpmath.fsub(Gamma1, ONE)), mpmath.power(rho1, Gamma1)))
    a2 = mpmath.fsub(mpmath.fsub(mpmath.fdiv(E1, rho1), ONE), mpmath.fmul(mpmath.fdiv(
        K2, mpmath.fsub(Gamma2, ONE)), mpmath.power(rho1, mpmath.fsub(Gamma2, ONE))))
    E2 = mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a2), rho2), mpmath.fmul(
        mpmath.fdiv(K2, mpmath.fsub(Gamma2, ONE)), mpmath.power(rho2, Gamma2)))
    a3 = mpmath.fsub(mpmath.fsub(mpmath.fdiv(E2, rho2), ONE), mpmath.fmul(mpmath.fdiv(
        K3, mpmath.fsub(Gamma3, ONE)), mpmath.power(rho2, mpmath.fsub(Gamma3, ONE))))

    parameters = (rho0, rho1, rho2, P0, P1, P2, E0, E1, E2, Gamma0,
                  Gamma1, Gamma2, Gamma3, K0, K1, K2, K3, a0, a1, a2, a3)

    return parameters


def P_of_rho_piecewise(self, rho, PPEoS):
    rho0, rho1, rho2, P0, P1, P2, E0, E1, E2, Gamma0, Gamma1, Gamma2, Gamma3, K0, K1, K2, K3, a0, a1, a2, a3 = PPEoS

    if rho < rho0:
        return mpmath.fmul(K0, mpmath.power(rho, Gamma0))
    if rho0 <= rho < rho1:
        return mpmath.fmul(K1, mpmath.power(rho, Gamma1))
    if rho1 <= rho < rho2:
        return mpmath.fmul(K2, mpmath.power(rho, Gamma2))

    return mpmath.fmul(K3, mpmath.power(rho, Gamma3))


def E_of_rho_piecewise(self, rho, PPEoS):
    rho0, rho1, rho2, P0, P1, P2, E0, E1, E2, Gamma0, Gamma1, Gamma2, Gamma3, K0, K1, K2, K3, a0, a1, a2, a3 = PPEoS

    if rho < rho0:
        return mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a0), rho), mpmath.fmul(mpmath.fdiv(K0, mpmath.fsub(Gamma0, ONE)), mpmath.power(rho, Gamma0)))
    if rho0 <= rho < rho1:
        return mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a1), rho), mpmath.fmul(mpmath.fdiv(K1, mpmath.fsub(Gamma1, ONE)), mpmath.power(rho, Gamma1)))
    if rho1 <= rho < rho2:
        return mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a2), rho), mpmath.fmul(mpmath.fdiv(K2, mpmath.fsub(Gamma2, ONE)), mpmath.power(rho, Gamma2)))

    return mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a3), rho), mpmath.fmul(mpmath.fdiv(K3, mpmath.fsub(Gamma3, ONE)), mpmath.power(rho, Gamma3)))


def rho_of_P_piecewise(self, P, PPEoS):
    rho0, rho1, rho2, P0, P1, P2, E0, E1, E2, Gamma0, Gamma1, Gamma2, Gamma3, K0, K1, K2, K3, a0, a1, a2, a3 = PPEoS

    if P < P0:
        return mpmath.power(mpmath.fdiv(P, K0), mpmath.fdiv(ONE, Gamma0))
    if P0 <= P < P1:
        return mpmath.power(mpmath.fdiv(P, K1), mpmath.fdiv(ONE, Gamma1))
    if P1 <= P < P2:
        return mpmath.power(mpmath.fdiv(P, K2), mpmath.fdiv(ONE, Gamma2))

    return mpmath.power(mpmath.fdiv(P, K3), mpmath.fdiv(ONE, Gamma3))


def E_of_P_piecewise(self, P, PPEoS):
    rho0, rho1, rho2, P0, P1, P2, E0, E1, E2, Gamma0, Gamma1, Gamma2, Gamma3, K0, K1, K2, K3, a0, a1, a2, a3 = PPEoS

    if P < P0:
        return mpmath.re(mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a0), mpmath.power(mpmath.fdiv(P, K0), mpmath.fdiv(ONE, Gamma0))), mpmath.fdiv(P, mpmath.fsub(Gamma0, ONE))))

    if P0 <= P < P1:
        return mpmath.re(mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a1), mpmath.power(mpmath.fdiv(P, K1), mpmath.fdiv(ONE, Gamma1))), mpmath.fdiv(P, mpmath.fsub(Gamma1, ONE))))

    if P1 <= P < P2:
        return mpmath.re(mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a2), mpmath.power(mpmath.fdiv(P, K2), mpmath.fdiv(ONE, Gamma2))), mpmath.fdiv(P, mpmath.fsub(Gamma2, ONE))))

    return mpmath.re(mpmath.fadd(mpmath.fmul(mpmath.fadd(ONE, a3), mpmath.power(mpmath.fdiv(P, K3), mpmath.fdiv(ONE, Gamma3))), mpmath.fdiv(P, mpmath.fsub(Gamma3, ONE))))


def dEdP_of_P_piecewise(self, P, PPEoS):
    rho0, rho1, rho2, P0, P1, P2, E0, E1, E2, Gamma0, Gamma1, Gamma2, Gamma3, K0, K1, K2, K3, a0, a1, a2, a3 = PPEoS

    def dEdP(P, Gamma, K, a):
        A = mpmath.fadd(ONE, a)

        # Use the existing rho_of_P function to safely compute density
        rho = self.rho_of_P_piecewise(P, PPEoS)

        # Check if rho is complex
        if hasattr(rho, 'imag') and mpmath.fabs(rho.imag) > mpmath.mpf('1e-30'):
            print(f"WARNING: rho is complex: {rho}")
            print(f"DEBUG: P={P}, Gamma={Gamma}, K={K}, a={a}")
            # Fallback to direct calculation with real part only
            rho = rho.real

        # Now use the simpler formula: dE/dP = (1+a)*rho/(Gamma*P) + 1/(Gamma-1)
        term1 = mpmath.fdiv(mpmath.fmul(A, rho), mpmath.fmul(Gamma, P))
        term2 = mpmath.fdiv(ONE, mpmath.fsub(Gamma, ONE))
        result = mpmath.fadd(term1, term2)
        return result

    if P < P0:
        return dEdP(P, Gamma0, K0, a0)
    if P0 <= P < P1:
        return dEdP(P, Gamma1, K1, a1)
    if P1 <= P < P2:
        return dEdP(P, Gamma2, K2, a2)
    return dEdP(P, Gamma3, K3, a3)

"""
Simple Polytropic Equation of State
"""

def P_of_rho_sim(self, rho, PEoS):

    Gamma, K = PEoS

    rho = mpmath.mpf(rho)
    Gamma = mpmath.mpf(Gamma)
    K = mpmath.mpf(K)

    return mpmath.fmul(K, mpmath.power(rho, Gamma))

def rho_of_P_sim(self, P, PEoS):

    Gamma, K = PEoS

    P = mpmath.mpf(P)
    Gamma = mpmath.mpf(Gamma)
    K = mpmath.mpf(K)

    exponent = mpmath.fdiv(mpmath.mpf('1'), Gamma)
    ratio = mpmath.fdiv(P, K)

    return mpmath.power(ratio, exponent)

def E_of_rho_sim(self, rho, PEoS):

    Gamma, K = PEoS

    rho = mpmath.mpf(rho)
    Gamma = mpmath.mpf(Gamma)
    K = mpmath.mpf(K)

    Gamma_minus_one = mpmath.fsub(Gamma, mpmath.mpf('1'))
    prefactor = mpmath.fdiv(K, Gamma_minus_one)

    rho_gamma = mpmath.power(rho, Gamma)

    return mpmath.fadd(rho, mpmath.fmul(prefactor, rho_gamma))

def E_of_P_sim(self, P, PEoS):

    Gamma, K = PEoS

    P = mpmath.mpf(P)
    Gamma = mpmath.mpf(Gamma)
    K = mpmath.mpf(K)

    exponent = mpmath.fdiv(mpmath.mpf('1'), Gamma)
    ratio = mpmath.fdiv(P, K)

    rho_term = mpmath.power(ratio, exponent)

    Gamma_minus_one = mpmath.fsub(Gamma, mpmath.mpf('1'))
    pressure_term = mpmath.fdiv(P, Gamma_minus_one)

    return mpmath.fadd(rho_term, pressure_term)

def dEdP_of_P_sim(self, P, PEoS):

    Gamma, K = PEoS

    P = mpmath.mpf(P)
    Gamma = mpmath.mpf(Gamma)
    K = mpmath.mpf(K)

    rho = self.rho_of_P_sim(P, PEoS)

    A = ONE
    term1 = mpmath.fdiv(mpmath.fmul(A, rho), mpmath.fmul(Gamma, P))
    term2 = mpmath.fdiv(ONE, mpmath.fsub(Gamma, ONE))

    return mpmath.fadd(term1, term2)

"""
General Functions
"""


def P_of_rho_general(self, c_den, EoS_param):
    if EoS_param is None:
        c_den = self.check_density(c_den)
        return mpmath.mpf(self.press(c_den))
    elif len(EoS_param) == 21:
        return mpmath.mpf(self.P_of_rho_piecewise(c_den, EoS_param))
    elif len(EoS_param) == 2:
        return mpmath.mpf(self.P_of_rho_sim(c_den, EoS_param))
    else:
        raise ValueError("Invalid EoS_param")


def E_of_rho_general(self, c_den, EoS_param):
    if EoS_param is None:
        c_den = self.check_density(c_den)
        return mpmath.mpf(c_den)
    elif len(EoS_param) == 21:
        return mpmath.mpf(self.E_of_rho_piecewise(c_den, EoS_param))
    elif len(EoS_param) == 2:
        return mpmath.mpf(self.E_of_rho_sim(c_den, EoS_param))
    else:
        raise ValueError("Invalid EoS_param")


def rho_of_P_general(self, P, EoS_param):
    if EoS_param is None:
        P = self.check_pressure(P)
        return mpmath.mpf(self.en_dens(P))
    elif len(EoS_param) == 21:
        return mpmath.mpf(self.rho_of_P_piecewise(P, EoS_param))
    elif len(EoS_param) == 2:
        return mpmath.mpf(self.rho_of_P_sim(P, EoS_param))
    else:
        raise ValueError("Invalid EoS_param")


def E_of_P_general(self, P, EoS_param):
    if EoS_param is None:
        P = self.check_pressure(P)
        return mpmath.mpf(self.en_dens(P))
    elif len(EoS_param) == 21:
        return mpmath.mpf(self.E_of_P_piecewise(P, EoS_param))
    elif len(EoS_param) == 2:
        return mpmath.mpf(self.E_of_P_sim(P, EoS_param))
    else:
        raise ValueError("Invalid EoS_param")


def dEdP_of_P_general(self, P, EoS_param):
    if EoS_param is None:
        raise ValueError("For a tabulated EoS we need a tabulated dEdP")
    elif len(EoS_param) == 21:
        return mpmath.mpf(self.dEdP_of_P_piecewise(P, EoS_param))
    elif len(EoS_param) == 2:
        return mpmath.mpf(self.dEdP_of_P_sim(P, EoS_param))
    else:
        raise ValueError("Invalid EoS_param")
