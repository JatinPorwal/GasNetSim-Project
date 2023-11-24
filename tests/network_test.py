#!/usr/bin/env python
# coding: utf-8

# **********************************************************************************************************************
# This test script serves to confirm that the Network class maintains its functionality post-modifications.
# **********************************************************************************************************************

import GasNetSim as gns
from pathlib import Path
from numpy.testing import assert_almost_equal
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def network_test():

    # Create a network instance with Irish13
    # Initialize the network with nodes and connections from a CSV file in the current directory
    network = gns.create_network_from_csv(Path('../examples/Irish13/.'))

    # Simulate the network to compute the pressures and flows
    network.simulation(tol=0.0000001)

    # Define the expected final pressure values for comparison
    expected_final_pressure = [7000000., 7000000., 7000000., 6933669.72289829,
                               6835569.60774993, 6623690.36661754, 6605743.43571843,
                               6602395.62670504, 6600915.11259321, 6592750.20799817,
                               6795230.64383087, 6791480.37767578, 6753737.3583671]

    # Retrieve the final pressure values from the simulated network nodes
    final_pressure = [node.pressure for node in network.nodes.values()]

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
    assert_almost_equal(final_pressure, expected_final_pressure)
    assert_almost_equal(inlet_flows, outlet_flows)
    assert_almost_equal(total_inflow, total_outflow)

    # If the assertion passes, print a message indicating that the test passed
    logger.info(f"Test passed: Results match the expected values.")


# Check if the script is being run directly
if __name__ == "__main__":
    network_test()
