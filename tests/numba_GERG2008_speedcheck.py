
from GasNetSim.components.utils.gas_mixture.GERG2008.gerg2008 import *
from scipy.constants import bar
from numpy.testing import assert_almost_equal, assert_allclose
from tests.gerg2008_numba import *
import time


def speed_heating_value():
    """
        Speed check the numba version of the CalculateHeatingValue function of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = np.array([0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001])
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Measure the execution time
    start_time = time.time()
    expected_heating_value = gas_mixture.CalculateHeatingValue(comp=nist_gas_mixture, hhv=True, parameter="volume")
    end_time = time.time()
    function_time = end_time - start_time
    print(f"CalculateHeatingValue took {function_time:.6f} seconds.")
    molarmass = gas_mixture.MolarMass
    molardensity = gas_mixture.MolarDensity
    # Measure the execution time
    start_time = time.time()
    calculated_heating_value = CalculateHeatingValue_numba(MolarMass=molarmass, MolarDensity=molardensity, comp=nist_gas_mixture, hhv=True, parameter="volume")
    end_time = time.time()
    function_time = end_time - start_time
    print(f"CalculateHeatingValue_numba took {function_time:.6f} seconds.")


def speed_convert_composition_gerg():
    """
    Speed check the numba version of the ConvertCompositionGERG method of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}

    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = np.array([0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001])
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Test the ConvertCompositionGERG function
    # Measure the execution time
    start_time = time.time()
    expected_result = gas_mixture.CovertCompositionGERG(nist_gas_mixture)
    end_time = time.time()
    function_time = end_time - start_time
    print(f"CovertCompositionGERG took {function_time:.6f} seconds.")
    expected_result.pop(0)

    # Calculate the converted composition using ConvertCompositionGERG method
    # Measure the execution time
    start_time = time.time()
    converted_composition = CovertCompositionGERG_numba(nist_gas_mixture)
    end_time = time.time()
    function_time = end_time - start_time
    print(f"CovertCompositionGERG_numba took {function_time:.6f} seconds.")


def speed_molarmass_gerg():
    """
        Speed check the numba version of the MolarMassGERG method of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = np.array([0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015,
                  0.00009, 0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001])
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Calculate the expected molar mass manually based on the given mixture
    # Measure the execution time
    start_time = time.time()
    expected_molar_mass = gas_mixture.MolarMassGERG()
    end_time = time.time()
    function_time = end_time - start_time
    print(f"MolarMassGERG took {function_time:.6f} seconds.")

    # Get the calculated molar mass from the MolarMassGERG method
    # Measure the execution time
    start_time = time.time()
    calculated_molar_mass = MolarMassGERG_numba(b)
    end_time = time.time()
    function_time = end_time - start_time
    print(f"MolarMassGERG_numba took {function_time:.6f} seconds.")


def speed_pressure_gerg():
    """
        Speed check the numba version of the PressureGERG method of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = np.array([0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001])
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Define the density input for PressureGERG method
    d = 10

    # Calculate the expected pressure using an example formula or method
    # Measure the execution time
    start_time = time.time()
    expected_values = gas_mixture.PressureGERG(d)
    end_time = time.time()
    function_time = end_time - start_time
    print(f"PressureGERG took {function_time:.6f} seconds.")
    Temp = gas_mixture.T
    AR = np.array(gas_mixture.AlpharGERG(itau=0, idelta=0, D=d))
    # Call the PressureGERG method with the given diameter
    # Measure the execution time
    start_time = time.time()
    calculated_values = PressureGERG_numba(AR, Temp, d)
    end_time = time.time()
    function_time = end_time - start_time
    print(f"PressureGERG_numba took {function_time:.6f} seconds.")


def speed_density_gerg():
    """
        Speed check the numba version of the DensityGERG function of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = np.array([0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001])
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Define the density input for PressureGERG method
    d = gas_mixture.MolarDensity

    # Expected value calculated from the function call
    # Measure the execution time
    start_time = time.time()
    _, _, expected_values = gas_mixture.DensityGERG()
    end_time = time.time()
    function_time = end_time - start_time
    print(f"DensityGERG took {function_time:.6f} seconds.")

    AR = np.array(gas_mixture.AlpharGERG(itau=0, idelta=0, D=d))
    Press = gas_mixture.P
    Temp = gas_mixture.T

    # Test the DensityGERG function with iFlag=0 (default)
    # Measure the execution time
    start_time = time.time()
    _, _, calculated_values = DensityGERG_numba(AR, Press, Temp, b, iFlag=0)  # Calling the function without any argument
    end_time = time.time()
    function_time = end_time - start_time
    print(f"DensityGERG_numba took {function_time:.6f} seconds.")


def speed_alpha0_gerg():
    """
        Speed Check the numba version of the Alpha0GERG() function of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = np.array([0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001])
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Expected value calculated from the function call
    # a0(0) - Ideal gas Helmholtz energy (all dimensionless [i.e., divided by RT])
    # a0(1) - tau*partial(a0)/partial(tau)
    # a0(2) - tau^2*partial^2(a0)/partial(tau)^2
    # Measure the execution time
    start_time = time.time()
    expected_alpha0 = gas_mixture.Alpha0GERG()
    end_time = time.time()
    function_time = end_time - start_time
    print(f"Alpha0GERG took {function_time:.6f} seconds.")

    Temp = gas_mixture.T
    MolarDensity = gas_mixture.MolarDensity
    X = b

    # Call the Alpha0GERG function
    # Measure the execution time
    start_time = time.time()
    actual_alpha0 = Alpha0GERG_numba(Temp, MolarDensity, X)
    end_time = time.time()
    function_time = end_time - start_time
    print(f"Alpha0GERG_numba took {function_time:.6f} seconds.")


def speed_reducing_parameters_gerg():
    """
        Speed check the numba version of the ReducingParametersGERG() function of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = np.array([0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015,
                  0.00009, 0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001])
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Expected value calculated from the function call
    # Measure the execution time
    start_time = time.time()
    expected_reducingparametersgerg = gas_mixture.ReducingParametersGERG()
    end_time = time.time()
    function_time = end_time - start_time
    print(f"ReducingParametersGERG took {function_time:.6f} seconds.")

    # Call the ReducingParametersGERG function
    # Tr - Reducing temperature(K)
    # Dr - Reducing density(mol / l)
    # Measure the execution time
    start_time = time.time()
    actual_reducingparametersgerg = ReducingParametersGERG_numba(b)
    end_time = time.time()
    function_time = end_time - start_time
    print(f"ReducingParametersGERG_numba took {function_time:.6f} seconds.")


def speed_pseudo_critical_point_gerg():
    """
            Speed check the numba version of the PseudoCriticalPointGERG() function of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = np.array([0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015,
                  0.00009, 0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001])
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Expected value calculated from the function call
    # Measure the execution time
    start_time = time.time()
    expected_pseudocriticalpointgerg = gas_mixture.PseudoCriticalPointGERG()
    end_time = time.time()
    function_time = end_time - start_time
    print(f"PseudoCriticalPointGERG took {function_time:.6f} seconds.")

    # Call the ReducingParametersGERG function
    # Measure the execution time
    start_time = time.time()
    actual_pseudocriticalpointgerg = PseudoCriticalPointGERG_numba(b)
    end_time = time.time()
    function_time = end_time - start_time
    print(f"PseudoCriticalPointGERG_numba took {function_time:.6f} seconds.")


def speed_alphar_gerg():
    """
            Speed check the numba version of the AlpharGERG() function of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = np.array([0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001])
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Expected value calculated from the function call
    #                         ar(0,0) - Residual Helmholtz energy (dimensionless, =a/RT)
    #                         ar(0,1) -     delta*partial  (ar)/partial(delta)
    #                         ar(0,2) -   delta^2*partial^2(ar)/partial(delta)^2
    #                         ar(0,3) -   delta^3*partial^3(ar)/partial(delta)^3
    #                         ar(1,0) -       tau*partial  (ar)/partial(tau)
    #                         ar(1,1) - tau*delta*partial^2(ar)/partial(tau)/partial(delta)
    #                         ar(2,0) -     tau^2*partial^2(ar)/partial(tau)^2

    D = 15.03402741629294
    # Measure the execution time
    start_time = time.time()
    expected_alphargerg = gas_mixture.AlpharGERG(1, 0, D)
    end_time = time.time()
    function_time = end_time - start_time
    print(f"AlpharGERG took {function_time:.6f} seconds.")

    Temp = gas_mixture.T

    # Call the ReducingParametersGERG function
    # Measure the execution time
    start_time = time.time()
    actual_alphargerg = AlpharGERG_numba(Temp, b, 1, 0, D)
    end_time = time.time()
    function_time = end_time - start_time
    print(f"AlpharGERG_numba took {function_time:.6f} seconds.")


def main():
    print("Running speed tests...\n")

    # Speed test for CalculateHeatingValue function
    speed_heating_value()

    # Speed test for ConvertCompositionGERG method
    speed_convert_composition_gerg()

    # Speed test for MolarMassGERG method
    speed_molarmass_gerg()

    # Speed test for PressureGERG method
    speed_pressure_gerg()

    # Speed test for DensityGERG function
    speed_density_gerg()

    # Speed test for Alpha0GERG function
    speed_alpha0_gerg()

    # Speed test for ReducingParametersGERG function
    speed_reducing_parameters_gerg()

    # Speed test for PseudoCriticalPointGERG function
    speed_pseudo_critical_point_gerg()

    # Speed test for AlpharGERG function
    speed_alphar_gerg()

if __name__ == "__main__":
    main()