from numba import njit, float64, types
from numba.extending import overload
from GasNetSim.components.utils.gas_mixture.GERG2008.gerg2008 import *
from tests.global_variables import *


@njit(float64(float64))
def Tanh_numba(xx):
    return (math.exp(xx) - math.exp(-xx)) / (math.exp(xx) + math.exp(-xx))


@njit(float64(float64))
def Sinh_numba(xx):
    return (math.exp(xx) - math.exp(-xx)) / 2


@njit(float64(float64))
def Cosh_numba(xx):
    return (math.exp(xx) + math.exp(-xx)) / 2


@njit(float64(float64[:]))
def MolarMassGERG_numba(x):
    """
        Calculate the molar mass of a gas mixture using GERG-2008 reference equation.

            Inputs:
                x:  Composition (mole fraction)
                    This should be an array(of format np.array()) containing the mole fractions of each component.
                    Ensure the sum of the compositions in the x array is equal to one.
                    The order of the fluids in this array must correspond to MMiGERG.

            return:
                mm: Molar mass (g/mol)
    """
    Mm = 0.0
    for ii in range(NcGERG):
        Mm += x[ii] * MMiGERG[ii]
    return Mm


@njit(types.UniTuple(float64, 2)(float64[:]))
def PseudoCriticalPointGERG_numba(x):
    """
        Calculate a pseudo critical point as the mole fraction average of the critical temperatures and volumes.

            Inputs:
                x:   Composition (mole fraction)
                     This should be an array(of format np.array()) containing the mole fractions of each component.
                     Ensure the sum of the compositions in the x array is equal to one.
                     The order of the fluids in this array must correspond to MMiGERG.

            return:
                Tcx: Pseudo-critical temperature
                Dcx: Pseudo-critical density
    """
    Vcx = 0
    Tcx = 0
    Dcx = 0

    for ii in range(NcGERG):
        Tcx = Tcx + x[ii] * Tc[ii]
        Vcx = Vcx + x[ii] / Dc[ii]

    if Vcx > epsilon:
        Dcx = 1 / Vcx

    return Tcx, Dcx


def ReducingParametersGERG_numba(x):
    """
        Function to calculate reducing parameters in GERG equation of state.

        Inputs:
            x:   Composition (mole fraction)
                 This should be an array(of format np.array()) containing the mole fractions of each component.
                 Ensure the sum of the compositions in the x array is equal to one.
                 The order of the fluids in this array must correspond to MMiGERG.

        return:
            Tr : Reduced temperature
            Dr : Reduced density
    """
    global xold, Trold, Drold
    Tr, Dr = ReducingParametersGERG_numba_sub(x)
    for i in range(NcGERG):
        xold[i] = x[i]
    Trold = Tr
    Drold = Dr
    return Tr, Dr


@njit(types.UniTuple(float64, 2)(float64[:]))
def ReducingParametersGERG_numba_sub(x):
    """
       Sub-function to calculate reducing parameters in GERG equation of state.
       Note: Not to be used directly at any other scripts.

       Inputs:
            x:   Composition (mole fraction)
                 This should be an array(of format np.array()) containing the mole fractions of each component.
                 Ensure the sum of the compositions in the x array is equal to one.
                 The order of the fluids in this array must correspond to MMiGERG.

        return:
            Tr : Reduced temperature
            Dr : Reduced density
    """

    icheck = 0
    for i in range(NcGERG):
        if abs(x[i] - xold[i]) > 0.0000001:
            icheck = 1

    if icheck == 0:
        return Trold, Drold

    Dr = 0
    Vr = 0
    Tr = 0
    for i in range(NcGERG):
        if x[i] > epsilon:
            F = 1
            for j in range(i, NcGERG):
                if x[j] > epsilon:
                    xij = F * (x[i] * x[j]) * (x[i] + x[j])
                    Vr = Vr + xij * gvij[i][j] / (bvij[i][j] * x[i] + x[j])
                    Tr = Tr + xij * gtij[i][j] / (btij[i][j] * x[i] + x[j])
                    F = 2
    if Vr > epsilon:
        Dr = 1. / Vr

    return Tr, Dr


def CovertCompositionGERG_numba(composition):
    pass


# @staticmethod
@overload(CovertCompositionGERG_numba)
def CovertCompositionGERG_numba(composition):
    """
        Converts a dictionary representing gas compositions into a GERG composition list.
        https://numba.readthedocs.io/en/stable/reference/pysupported.html#typed-dict

        Inputs:
            composition (dict): A dictionary containing gas species and their compositions.

        Returns:
            gerg_composition: A list representing the GERG composition of gases.
    """
    gerg_composition = [0.0] * 21
    global gerg_gas_spices

    for gas_spice, composition in composition.items():
        gerg_composition[gerg_gas_spices.index(gas_spice)] = composition

    return gerg_composition
