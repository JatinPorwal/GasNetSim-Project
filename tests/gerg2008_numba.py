from numba import njit, float64, types
from numba.extending import overload
from tests.global_variables import *
# from GasNetSim.components.utils.gas_mixture.GERG2008.gerg2008 import *


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

        return:
            gerg_composition (list): A list representing the GERG composition of gases.
    """
    gerg_composition = [0.0] * 21
    global gerg_gas_spices

    for gas_spice, composition in composition.items():
        gerg_composition[gerg_gas_spices.index(gas_spice)] = composition

    return gerg_composition


def CalculateHeatingValue_numba(MolarMass, MolarDensity, comp, hhv, parameter):
    pass


@overload(CalculateHeatingValue_numba)
def CalculateHeatingValue_numba(MolarMass, MolarDensity, comp, hhv, parameter):
    """
        Calculate the heating value of a gas mixture based on its composition and other properties.

        Inputs:
            MolarMass (float64): The molar mass of the gas mixture.
            MolarDensity (float64): The molar density of the gas mixture.
            comp (dict): A dictionary representing the composition of the gas mixture.
            hhv (bool): True for Higher Heating Value (HHV) calculation, False for Lower Heating Value (LHV) calculation.
            parameter (str): Specifies the parameter for heating value calculation. Options: 'mass' or 'volume'.

        return:
            heating_value (float64): The calculated heating value based on the provided parameters.
    """
    dict_components = {'methane': {'C': 1, 'H': 4},
                       'nitrogen': {'N': 2},
                       'carbon dioxide': {'C': 1, 'O': 2},
                       'ethane': {'C': 2, 'H': 6},
                       'propane': {'C': 3, 'H': 8},
                       'isobutane': {'C': 4, 'H': 10},
                       'n-butane': {'C': 4, 'H': 10},
                       'isopentane': {'C': 5, 'H': 12},
                       'n-pentane': {'C': 5, 'H': 12},
                       'n-hexane': {'C': 6, 'H': 14},
                       'n-heptane': {'C': 7, 'H': 16},
                       'n-octane': {'C': 8, 'H': 18},
                       'n-nonane': {'C': 9, 'H': 20},
                       'n-decane': {'C': 10, 'H': 22},
                       'hydrogen': {'H': 2},
                       'oxygen': {'O': 2},
                       'carbon monoxide': {'C': 1, 'O': 1},
                       'water': {'H': 2, 'O': 1},
                       'hydrogen sulfide': {'H': 2, 'S': 1},
                       'helium': {'He': 1},
                       'argon': {'Ar': 1}}
                       # 'SO2': {'S': 1, 'O': 2}}

    # 273 K
    dict_enthalpy_mole = {'methane': -75483.51423273719,
                          'nitrogen': 0.0,
                          'carbon dioxide': -394431.82606764464,
                          'ethane': -83856.2627150042,
                          'propane': -103861.117481869,
                          'isobutane': -135360.0,
                          'n-butane': -125849.99999999999,
                          'isopentane': -178400.0,
                          'n-pentane': -173500.0,
                          'n-hexane': -198490.0,
                          'n-heptane': -223910.0,
                          'n-hctane': -249730.0,
                          'n-nonane': -274700.0,
                          'n-decane': -300900.0,
                          'hydrogen': 0.0,
                          'oxygen': -4.40676212751828,
                          'carbon monoxide': -111262.34509634285,
                          'water': -242671.7203547155,
                          'hydrogen sulfide': -20600.0,
                          'helium': 0.0,
                          'argon': 0.0,
                          'carbon': 0.0}
                          # 'H': 218000.0,
                          # 'O': 249190.0,
                          # 'SO2': -296840.0}

    # reactants
    atom_list = []

    for key, value in comp.items():
        for key1, value1 in dict_components.items():
            if key == key1:
                # njit does not support yield operations
                # atom_list.append(dict((x, y * value) for x, y in dict_components.get(key1, {}).items()))
                atoms = {}
                for x, y in dict_components.get(key1, {}).items():
                    atoms[x] = y * value
                atom_list.append(atoms)
    # njit does not support the Counter()
    # reactants_atom = Counter()
    # for d in atom_list:
    #     reactants_atom.update(d)
    reactants_atom = {}
    for d in atom_list:
        for k, v in d.items():
            reactants_atom[k] = reactants_atom.get(k, 0) + v

    # products
    n_CO2 = reactants_atom["C"]
    n_SO2 = reactants_atom["S"]
    n_H2O = reactants_atom["H"] / 2
    products_dict = {'carbon dioxide': n_CO2, 'sulfur dioxide': n_SO2, 'water': n_H2O}

    # oxygen for complete combustion
    n_O = n_CO2 * 2 + n_SO2 * 2 + n_H2O * 1  # 2 is number of O atoms in CO2 AND SO2 and 1 is number of O atoms in H2O
    n_O2 = n_O / 2
    # https://numba.discourse.group/t/deep-copy-of-typed-lists/952
    # reactants_dict = deepcopy(comp)
    # reactants_dict.update({'oxygen': n_O2})
    # Manual copy of comp dictionary
    reactants_dict = {}
    for k, v in comp.items():
        reactants_dict[k] = v
    reactants_dict['oxygen'] = n_O2

    # LHV calculation
    LHV = 0
    for key, value in dict_enthalpy_mole.items():
        for key1, value1 in reactants_dict.items():
            if key == key1:
                LHV += value * value1
        for key2, value2 in products_dict.items():
            if key == key2:
                LHV -= value * value2

    # 298 K
    hw_liq = -285825.0
    hw_gas = -241820.0

    # 273 K
    # hw_liq = -287654.96084928664
    # hw_gas = -242628.01574091613

    HHV = LHV + (hw_gas - hw_liq) * products_dict["water"]

    if parameter == 'mass':
        # returns heating value in MJ/kg
        if hhv:
            heating_value = HHV / MolarMass * 1e3
        else:
            heating_value = LHV / MolarMass * 1e3
    else:
        # returns heating value in kJ/m3
        if hhv:
            heating_value = HHV * MolarDensity
        else:
            heating_value = LHV * MolarDensity

    return heating_value
