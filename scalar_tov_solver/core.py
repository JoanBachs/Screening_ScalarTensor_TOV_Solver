# core.py
"""
Import Libraries and Modules
"""

import io
import os
import sys
import numpy
import mpmath
import argparse
import pkg_resources
import matplotlib.pyplot as plt

"""
Precision control
"""

DEFAULT_PRECISION = 100

def set_precision(dps: int = DEFAULT_PRECISION):
    """
    Set global mpmath precision.

    Parameters
    ----------
    dps : int
        Decimal places of precision.
    """
    mpmath.mp.dps = dps

# Set default precision on import
set_precision()

"""
Numerical constants
"""

ZERO = mpmath.mpf('0')
HALF = mpmath.mpf('0.5')
ONE = mpmath.mpf('1')
TWO = mpmath.mpf('2')
THREE = mpmath.mpf('3')
FOUR = mpmath.mpf('4')
SIX = mpmath.mpf('6')
EIGHT = mpmath.mpf('8')

"""
Physical constants (code units)
"""

GeV = ONE  # We set the code units to GeV
cm = mpmath.fdiv(mpmath.mpf('5.06e13'), GeV)  # 1/GeV
gram = mpmath.fmul(mpmath.mpf('5.62e23'), GeV)  # GeV
second = mpmath.fdiv(mpmath.mpf('1.52e24'), GeV)  # 1/GeV
Grav_C = mpmath.fdiv(mpmath.mpf('6.708e-39'), mpmath.fmul(GeV, GeV))  # 1/GeV^2
clight = ONE
hbar = ONE

pi = mpmath.pi
Msun = mpmath.fmul(mpmath.mpf('1.989e33'), gram)
km = mpmath.fmul(mpmath.mpf('1e5'), cm)
Mp = mpmath.fdiv(mpmath.fmul(hbar, clight), mpmath.sqrt(
    mpmath.fprod([EIGHT, pi, Grav_C])))
kappa = mpmath.fdiv(ONE, Mp)
kappa2 = mpmath.power(kappa, TWO)

# Conversion factors
MeV_fm3_to_code_units = mpmath.fdiv(mpmath.fmul(mpmath.mpf(
    '1e-3'), GeV), mpmath.power(mpmath.fmul(mpmath.mpf('1e-13'), cm), THREE))
g_cm3_to_code_units = mpmath.fdiv(gram, mpmath.power(cm, THREE))
dyn_cm2_to_code_units = mpmath.fdiv(gram, mpmath.fprod([cm, second, second]))