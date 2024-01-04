from numba import njit, float64, types, int64
from numba.extending import overload
from tests.global_variables import *
#from GasNetSim.components.utils.gas_mixture.GERG2008.gerg2008 import *


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


@njit
def Alpha0GERG_numba(Temp, MolarDensity, X):
    """
            Private Sub Alpha0GERG(T, D, x, a0)

            Calculate the ideal gas Helmholtz energy and its derivatives with respect to tau and delta.
            This routine is not needed when only P (or Z) is calculated.
            Inputs:
                T: Temperature (K)
                D: Density (mol/l)
                x: Composition (mole fraction)
            return:
                a0:        a0(0) - Ideal gas Helmholtz energy (all dimensionless [i.e., divided by RT])
                           a0(1) - tau*partial(a0)/partial(tau)
                           a0(2) - tau^2*partial^2(a0)/partial(tau)^2
            """
    T = Temp
    D = MolarDensity
    x = X
    th0T = 0.0
    LogxD = 0.0
    SumHyp0 = 0.0
    SumHyp1 = 0.0
    SumHyp2 = 0.0
    hcn = 0.0
    hsn = 0.0

    a0 = [0.0] * 3
    if D > epsilon:
        LogD = math.log(D)
    else:
        LogD = math.log(epsilon)
    LogT = math.log(T)
    for i in range(NcGERG):
        if x[i] > epsilon:
            LogxD = LogD + math.log(x[i])
            SumHyp0 = 0
            SumHyp1 = 0
            SumHyp2 = 0
        for j in range(3, 7):
            if th0i[i][j] > epsilon:
                th0T = th0i[i][j] / T
                ep = math.exp(th0T)
                em = 1 / ep
                hsn = (ep - em) / 2
                hcn = (ep + em) / 2
                if j == 3 or j == 5:
                    LogHyp = math.log(abs(hsn))
                    SumHyp0 = SumHyp0 + n0i[i][j] * LogHyp
                    SumHyp1 = SumHyp1 + n0i[i][j] * th0T * hcn / hsn
                    SumHyp2 = SumHyp2 + n0i[i][j] * (th0T / hsn) * (th0T / hsn)
                else:
                    LogHyp = math.log(abs(hcn))
                    SumHyp0 = SumHyp0 - n0i[i][j] * LogHyp
                    SumHyp1 = SumHyp1 - n0i[i][j] * th0T * hsn / hcn
                    SumHyp2 = SumHyp2 + n0i[i][j] * (th0T / hcn) * (th0T / hcn)

        a0[0] += +x[i] * (LogxD + n0i[i][0] + n0i[i][1] / T - n0i[i][2] * LogT + SumHyp0)
        a0[1] += +x[i] * (n0i[i][2] + n0i[i][1] / T + SumHyp1)
        a0[2] += -x[i] * (n0i[i][2] + SumHyp2)
    return a0

def tTermsGERG_numba(lntau, x):
    """
        Private Sub tTermsGERG(lntau, x)
        Calculate temperature dependent parts of the GERG-2008 equation of state
        Inputs:
            lntau:  tau = Tr / T => lntau = math.log(tau)
            x:      Composition (mole fraction)
        return:
            null
    """
    global taup, taupijk
    taup, taupijk = tTermsGERG_numba_sub(taup, taupijk, lntau, x)


@njit
def tTermsGERG_numba_sub(taup, taupijk, lntau, x):
    """
        Calculate temperature-dependent parts of the GERG-2008 equation of state.

        Inputs:
            taup :    List containing calculated temperature-dependent values for taup.
            taupijk : List containing calculated temperature-dependent values for taupijk.
            lntau :   Natural logarithm of tau, a term used in the calculation.
            x :       Composition (mole fraction) of the components.

        returns:
            taup :    Updated taup values.
            taupijk : Updated taupijk values.
    """
    taup0 = [0] * (12)

    i = 4  # Use propane to get exponents for short form of EOS
    for k in range(int(kpol[i] + kexp[i])):  # for (int k = 1; k <= kpol[i] + kexp[i]; ++k)
        taup0[k] = math.exp(toik[i][k] * lntau)
    for i in range(NcGERG):  # for (int i = 1; i <= NcGERG; ++i)
        if x[i] > epsilon:
            if (i > 3) and (i != 14) and (i != 17) and (i != 19):
                for k in range(int(kpol[i] + kexp[i])):  # for (int k = 1; k <= kpol[i] + kexp[i]; ++k)
                    taup[i][k] = noik[i][k] * taup0[k]
            else:
                for k in range(int(kpol[i] + kexp[i])):  # for (int k = 1; k <= kpol[i] + kexp[i]; ++k)
                    taup[i][k] = noik[i][k] * math.exp(toik[i][k] * lntau)

    for i in range(NcGERG-1):  # for (int i = 1; i <= NcGERG - 1; ++i)
        if x[i] > epsilon:
            for j in range(i + 1, NcGERG):  # for (int j = i + 1; j <= NcGERG; ++j)
                if x[j] > epsilon:
                    mn = int(mNumb[i][j])
                    if mn >= 0:
                        for k in range(int(kpolij[mn])):  # for (int k = 1; k <= kpolij[mn]; ++k)
                            taupijk[mn][k] = nijk[mn][k] * math.exp(tijk[mn][k] * lntau)

    return taup, taupijk


@njit(types.UniTuple(float64, 3)(float64[:, :], float64, float64))
def PressureGERG_numba(ar, T, D):
    """
    Sub PressureGERG(T, D, x, P, Z)

    Calculate pressure as a function of temperature and density.  The derivative d(P)/d(D) is also calculated
    for use in the iterative DensityGERG subroutine (and is only returned as a common variable).

    :return:        P: Pressure (kPa)
                    Z: Compressibility factor
                    dPdDsave - d(P)/d(D) [kPa/(mol/l)] (at constant temperature)
    //          - This variable is cached in the common variables for use in the iterative density solver, but not returned as an argument.
    """
    #ar = self.AlpharGERG_numba(itau=0, idelta=0, D=D)

    Z = 1 + ar[0][1]
    P = D * RGERG * T * Z
    dPdDsave = RGERG * T * (1 + 2 * ar[0][1] + ar[0][2])
    return P, Z, dPdDsave