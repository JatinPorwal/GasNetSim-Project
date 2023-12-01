#!/usr/bin/env python
# coding: utf-8

# **********************************************************************************************************************
# This test script serves to confirm that the Network class maintains its functionality post-modifications.
# The test script contains 3 test functions:
#   1. test_network_volume_flow_rate_balance() - to ensure the volume flow rate balance within a network.
#   2. test_network_energy_flow_balance() - to ensure the energy flow balance within a network.
#   3. test_network_composition_balance() - to ensure the accuracy of gas composition within a network simulation.
# **********************************************************************************************************************

import os
import GasNetSim as gns
from pathlib import Path
from numpy.testing import assert_almost_equal
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def find_git_root(path):
    """
        Find the root path of the Git repository starting from the given path.

        Args:
        - path: The starting directory path to search from.

        Returns:
        - The root path of the Git repository or None if not found.
    """
    # Traverse up the directory tree until finding the .git folder
    while path != '/':
        if os.path.isdir(os.path.join(path, '.git')):
            return path
        path = os.path.dirname(path)
    return None


def test_network_volume_flow_rate_balance():
    """
        Test to ensure the volume flow rate balance within a network.
        Calculates the total inflow and outflow across all nodes in the network to verify conservation
        of volume flow rates.
    """

    # Find the current absolute path
    test_directory_path = os.path.abspath(os.getcwd())
    # Find the root path of the Git repository
    root_path = find_git_root(test_directory_path)
    new_path = os.path.join(root_path, 'examples', 'Irish13')

    # Create a network instance with Irish13
    # Initialize the network with nodes and connections from a CSV file in the current directory
    # network = gns.create_network_from_csv(Path('../examples/Irish13/.'))
    network = gns.create_network_from_csv(Path(new_path))

    # Simulate the network to compute the pressures and flows
    network.simulation(tol=0.0000001)

    # Define the expected final pressure values for comparison
    # expected_final_pressure = [7000000., 7000000., 7000000., 6933669.72289829,
    #                            6835569.60774993, 6623690.36661754, 6605743.43571843,
    #                            6602395.62670504, 6600915.11259321, 6592750.20799817,
    #                            6795230.64383087, 6791480.37767578, 6753737.3583671]

    # Retrieve the final pressure values from the simulated network nodes
    # final_pressure = [node.pressure for node in network.nodes.values()]

    # Calculate total inflow and outflow over the entire Network
    # Initialize lists to store node information
    node_indices = []
    inlet_flows = []
    outlet_flows = []

    # Calculate initial flow rates from the nodes
    initial_flows = {node_index: node.volumetric_flow if node.volumetric_flow is not None else 0
                     for node_index, node in network.nodes.items()}

    # Iterate through nodes in the network to gather inlet and outlet flow information
    for node_index, node in network.nodes.items():
        inlet_flow = 0
        outlet_flow = 0

        # Iterate through connections to find flows related to the current node
        for connection in network.connections.values():
            if connection.outlet_index == node_index:
                inlet_flow += connection.flow_rate
            elif connection.inlet_index == node_index:
                outlet_flow += connection.flow_rate

        # Add initial flow to the outlet flow for each node
        outlet_flow += initial_flows.get(node_index, 0)

        # Append node information to the lists
        node_indices.append(node_index)
        inlet_flows.append(inlet_flow)
        outlet_flows.append(outlet_flow)

    # Calculate total inflow and outflow for the entire network
    total_inflow = sum(inlet_flows)
    total_outflow = sum(outlet_flows)

    # Check if the final pressure values from the simulation match the expected values
    # assert_almost_equal(final_pressure, expected_final_pressure)
    assert_almost_equal(inlet_flows, outlet_flows)
    assert_almost_equal(total_inflow, total_outflow)

    # If the assertion passes, print a message indicating that the test passed
    logger.info(f"Test passed: Results match the expected values for volume_flow_rate_balance.")


def test_network_energy_flow_balance():
    """
        Test to ensure the energy flow balance within a network.
        Calculates the total energy flowing in and out of the network to verify energy conservation.
    """
    # Find the current absolute path
    test_directory_path = os.path.abspath(os.getcwd())
    # Find the root path of the Git repository
    root_path = find_git_root(test_directory_path)
    new_path = os.path.join(root_path, 'examples', 'Irish13')

    # Create a network instance with Irish13
    # Initialize the network with nodes and connections from a CSV file in the current directory
    # network = gns.create_network_from_csv(Path('../examples/Irish13/.'))
    network = gns.create_network_from_csv(Path(new_path))

    # Simulate the network to compute the pressures and flows
    network.simulation(tol=0.0000001)

    # Calculate total energy flow going into the network
    total_energy_in = sum([node.energy_flow for node in network.nodes.values() if node.energy_flow > 0])

    # Calculate total energy flow going out of the network
    total_energy_out = sum([node.energy_flow for node in network.nodes.values() if node.energy_flow < 0])

    # Assert that the total energy going into the network equals the total energy going out (within a tolerance)
    assert_almost_equal(total_energy_in, abs(total_energy_out))

    # If the assertion passes, print a message indicating that the test passed
    logger.info(f"Test passed: Results match the expected values for energy_flow_rate_balance.")


def test_network_composition_balance():

    """
        Test to ensure the accuracy of gas composition within a network simulation.
        Calculates the total mass flow rate and component wise mass flow rate to verify the composition balance.
    """
    # Find the current absolute path
    test_directory_path = os.path.abspath(os.getcwd())
    # Find the root path of the Git repository
    root_path = find_git_root(test_directory_path)
    new_path = os.path.join(root_path, 'examples', 'Irish13')

    # Create a network instance with Irish13
    # Initialize the network with nodes and connections from a CSV file in the current directory
    # network = gns.create_network_from_csv(Path('../examples/Irish13/.'))
    network = gns.create_network_from_csv(Path(new_path))

    # Simulate the network to compute the pressures and flows
    network.simulation(tol=0.0000001)

    # Calculate the composition balance for each pipeline
    for i, pipeline in network.pipelines.items():
        # Calculate the total mass flow rate
        mass_flow_rate = pipeline.flow_rate * pipeline.gas_mixture.density

        # Calculate the mass flow rate of each component
        component_flow_rates = {component: mole_fraction * mass_flow_rate
                                for component, mole_fraction in pipeline.gas_mixture.composition.items()}

        # Calculate the total component mass flow rate
        total_component_mass_flow_rate = sum(component_flow_rates.values())

        # Check if the total component mass flow rate equals the total mass flow rate within a tolerance
        assert_almost_equal(total_component_mass_flow_rate, mass_flow_rate)

    # If the assertion passes, print a message indicating that the test passed
    logger.info(f"Test passed: Results match the expected values for composition_balance.")
