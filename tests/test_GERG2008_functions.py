
from GasNetSim.components.utils.gas_mixture.GERG2008.gerg2008 import *
from scipy.constants import bar
from numpy.testing import assert_almost_equal
import numpy as np


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
    d = 10

    # Calculate the expected pressure using an example formula or method
    expected_pressure = 34030.02185084158
    expected_compressibility_factor = 1.0232165629652001
    expected_calculated_dpdd = 4665.32077103759

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


def test_alpha0_gerg():
    """
        Test the Alpha0GERG() function of GasMixtureGERG2008 class.
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
    # a0(0) - Ideal gas Helmholtz energy (all dimensionless [i.e., divided by RT])
    # a0(1) - tau*partial(a0)/partial(tau)
    # a0(2) - tau^2*partial^2(a0)/partial(tau)^2
    expected_alpha0 = [3.915632780302234, 0.27817599098351004, -4.4002653106914655]

    # Call the Alpha0GERG function
    actual_alpha0 = gas_mixture.Alpha0GERG()
    assert_almost_equal(actual_alpha0, expected_alpha0)


def test_reducing_parameters_gerg():
    """
        Test the ReducingParametersGERG() function of GasMixtureGERG2008 class.
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
    expected_reducingparametersgerg = [211.29730660311438, 9.389250212600038]

    # Call the ReducingParametersGERG function
    # Tr - Reducing temperature(K)
    # Dr - Reducing density(mol / l)
    actual_reducingparametersgerg = gas_mixture.ReducingParametersGERG()
    assert_almost_equal(actual_reducingparametersgerg, expected_reducingparametersgerg)


def test_pseudo_critical_point_gerg():
    """
            Test the PseudoCriticalPointGERG() function of GasMixtureGERG2008 class.
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
    expected_pseudocriticalpointgerg = [211.69335825999997, 9.378320910617676]

    # Call the ReducingParametersGERG function
    actual_pseudocriticalpointgerg = gas_mixture.PseudoCriticalPointGERG()
    assert_almost_equal(actual_pseudocriticalpointgerg, expected_pseudocriticalpointgerg)


def test_alphar_gerg():
    """
            Test the AlpharGERG() function of GasMixtureGERG2008 class.
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
    #                         ar(0,0) - Residual Helmholtz energy (dimensionless, =a/RT)
    #                         ar(0,1) -     delta*partial  (ar)/partial(delta)
    #                         ar(0,2) -   delta^2*partial^2(ar)/partial(delta)^2
    #                         ar(0,3) -   delta^3*partial^3(ar)/partial(delta)^3
    #                         ar(1,0) -       tau*partial  (ar)/partial(tau)
    #                         ar(1,1) - tau*delta*partial^2(ar)/partial(tau)/partial(delta)
    #                         ar(2,0) -     tau^2*partial^2(ar)/partial(tau)^2
    expected_alphargerg = np.array([[-0.12350275432589777, 0.02321656296520004, 0.3563380816212008, 0.2932603390400152],
                                    [-0.8745987057694501, -0.8246675326603303, 0.040189577810096334, 0],
                                    [-0.24620538247419213, 0, 0, 0],
                                    [0, 0, 0, 0]])

    # Call the ReducingParametersGERG function
    actual_alphargerg = gas_mixture.AlpharGERG(1, 0, 10)
    assert_almost_equal(actual_alphargerg, expected_alphargerg)
