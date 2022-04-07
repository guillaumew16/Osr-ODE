import argparse

import sympy as sym
from sympy.abc import F, z, s

from utils import DTA, GDA, EGM

def main():
    ODE_coeff_funs = EGM.Osr_ODE(r=3)
    sym.pprint(ODE_coeff_funs)

if __name__=="__main__":
    main()
