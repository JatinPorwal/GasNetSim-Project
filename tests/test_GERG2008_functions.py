
from GasNetSim.components.utils.gas_mixture.GERG2008.gerg2008 import *
from scipy.constants import bar
from numpy.testing import assert_almost_equal


# Test the tanh, sinh, and cosh functions
def test_tanh_sinh_cosh():
    """
        Test the hyperbolic tangent (tanh), hyperbolic sine (sinh), and hyperbolic cosine (cosh) functions.
    """
    assert_almost_equal(Tanh(0), 0.0)
    assert_almost_equal(Tanh(1), 0.7615941559557649)
    assert_almost_equal(Tanh(-1), -0.7615941559557649)

    assert_almost_equal(Sinh(0), 0.0)
    assert_almost_equal(Sinh(1), 1.1752011936438014)
    assert_almost_equal(Sinh(-1), -1.1752011936438014)

    assert_almost_equal(Cosh(0), 1.0)
    assert_almost_equal(Cosh(1), 1.5430806348152437)
    assert_almost_equal(Cosh(-1), 1.5430806348152437)


def test_heating_value():
    """
        Test the CalculateHeatingValue function of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = [0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001]
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Test the CalculateHeatingValue function
    expected_heating_value = 11452626.041728042
    calculated_heating_value = gas_mixture.CalculateHeatingValue(comp=nist_gas_mixture, hhv=True, parameter="volume")

    assert_almost_equal(calculated_heating_value, expected_heating_value)


def test_convert_composition_gerg():
    """
    Test the ConvertCompositionGERG method of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = [0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001]
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Test the ConvertCompositionGERG function
    expected_result = [0.0, 0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024,
                       0.00015, 9e-05, 0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001]

    # Calculate the converted composition using ConvertCompositionGERG method
    converted_composition = gas_mixture.CovertCompositionGERG(nist_gas_mixture)
    assert_almost_equal(converted_composition, expected_result)


def test_molarmass_gerg():
    """
        Test the MolarMassGERG method of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = [0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001]
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Calculate the expected molar mass manually based on the given mixture
    expected_molar_mass = 20.5427445016

    # Get the calculated molar mass from the MolarMassGERG method
    calculated_molar_mass = gas_mixture.MolarMassGERG()
    assert_almost_equal(expected_molar_mass, calculated_molar_mass)


def test_pressure_gerg():
    """
        Test the PressureGERG method of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = [0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001]
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Define the density input for PressureGERG method
    d = 1000

    # Calculate the expected pressure using an example formula or method
    expected_pressure = 2.5532549524167652e+16
    expected_compressibility_factor = 7677140991.083876
    expected_calculated_dpdd = 204183482532894.53

    # Call the PressureGERG method with the given diameter
    calculated_pressure, calculated_compressibility_factor, calculated_dpdd = gas_mixture.PressureGERG(d)
    assert_almost_equal(expected_pressure, calculated_pressure)
    assert_almost_equal(expected_compressibility_factor, calculated_compressibility_factor)
    assert_almost_equal(expected_calculated_dpdd, calculated_dpdd)


def test_density_gerg():
    """
        Test the DensityGERG function of GasMixtureGERG2008 class.
    """
    # Create the NIST gas mixture dictionary
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
         'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
         'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
         'helium', 'argon']
    b = [0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001]
    for ii in range(21):
        nist_gas_mixture[a[ii]] = b[ii]

    # Create an instance of the GasMixtureGERG2008 class with the NIST gas mixture
    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    # Expected value calculated from the function call
    expected_density = 12.79828626082062

    # Test the DensityGERG function with iFlag=0 (default)
    _, _, calculated_density = gas_mixture.DensityGERG()  # Calling the function without any argument
    assert_almost_equal(expected_density, calculated_density)
