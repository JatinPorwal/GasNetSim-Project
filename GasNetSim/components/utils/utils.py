#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2024.
#     Developed by Yifei Lu
#     Last change on 7/15/24, 2:54 PM
#     Last change by yifei
#    *****************************************************************************
import math
from pyparsing import col
from collections import OrderedDict
from scipy import sparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from timeit import default_timer as timer

try:
    import cupy as cp
    import cupy.sparse.linalg as cpsplinalg
except ImportError:
    # logging.warning(f"CuPy is not installed or not available!")
    print(f"CuPy is not installed or not available!")

from .cuda_support import create_matrix_of_zeros


def create_connection_matrix(n_nodes: int, components: dict, component_type: int,
                             use_cuda=False, sparse_matrix: bool = False):
    row_ind = list()
    col_ind = list()
    data = list()

    if not sparse_matrix:
        cnx = create_matrix_of_zeros(n_nodes, use_cuda=use_cuda, sparse_matrix=sparse_matrix)

    for comp in components.values():
        i = comp.inlet_index - 1
        j = comp.outlet_index - 1
        if sparse_matrix:
            row_ind.append(i)
            col_ind.append(j)
            data.append(component_type)
        else:
            cnx[i][j] = component_type
            cnx[j][i] = component_type

    if sparse_matrix:
        cnx = sparse.coo_matrix((data, (row_ind, col_ind)))
    return cnx


def levenberg_marquardt_damping_factor(m, s, b):
    return 10 ** (m * math.log10(s + b))


def delete_matrix_rows_and_columns(matrix, to_remove, use_cuda=False):
    new_matrix = matrix

    if use_cuda:
        new_matrix = cp.delete(new_matrix, to_remove, 0)  # delete rows
        new_matrix = cp.delete(new_matrix, to_remove, 1)  # delete columns
    else:
        new_matrix = np.delete(new_matrix, to_remove, 0)  # delete rows
        new_matrix = np.delete(new_matrix, to_remove, 1)  # delete columns

    return new_matrix


def jacobian_matrix_condition_number(matrix):
    print(f"The condition number of the matrix is {np.linalg.cond(matrix)}.")


def print_n_largest_absolute_values(n, values):
    sorted_values = sorted([abs(x) for x in values])
    print(sorted_values[-n::-1])
    return None


def calculate_nodal_inflow_states(nodes, connections, mapping_connections, flow_matrix, use_cuda=False):
    if use_cuda:
        nodal_total_inflow = cp.sum(cp.where(flow_matrix > 0, flow_matrix, 0), axis=1)
    else:
        nodal_total_inflow = np.sum(np.where(flow_matrix > 0, flow_matrix, 0), axis=1)


def gas_mixture_transportation(connection, time_step):
    composition_history = connection.composition_history
    batch_location_history = connection.batch_location_history
    length = connection.length
    velocity = connection.flow_velocity
    if velocity is None:
        velocity = 0
    if velocity >= 0:
        outlet_composition = connection.outlet.gas_mixture.composition
    else:
        outlet_composition = connection.inlet.gas_mixture.composition

    while batch_location_history[0] >= length:  # if the head of a batch reached the end of the pipeline
        outlet_composition = composition_history[0]
        composition_history, batch_location_history = composition_history[1:], batch_location_history[1:]

    batch_location_history = np.append(batch_location_history, 0)
    composition_history = np.append(composition_history, connection.inlet.gas_mixture.composition)
    batch_location_history += time_step * velocity

    # update connection composition and batch location history
    connection.composition_history = composition_history
    connection.batch_location_history = batch_location_history
    return outlet_composition, connection


def calculate_nodal_inflow_states(nodes, connections, mapping_connections, flow_matrix, composition_tracking=False,
                                  time_step=0):
    nodal_total_inflow = np.sum(np.where(flow_matrix > 0, flow_matrix, 0), axis=1)

    nodal_gas_inflow_composition = dict()
    nodal_gas_inflow_temperature = dict()

    for i_node, node in nodes.items():  # iterate over all nodes
        if use_cuda:
            inflow_from_node = cp.where(flow_matrix[i_node - 1] > 0)[0]  # find the supplying nodes
        else:
            inflow_from_node = np.where(flow_matrix[i_node-1] > 0)[0]  # find the supplying nodes

        # TODO: check inflow nodes

        if len(inflow_from_node) == 0:
            pass
        else:
            inflow_from_node += 1

        total_inflow_comp = dict()
        total_inflow = nodal_total_inflow[i_node-1]
        total_inflow_temperature_times_flow_rate = 0

        for inlet_index in inflow_from_node:
            if type(inlet_index) is not int:
                inlet_index = inlet_index.item()
            gas_composition = nodes[inlet_index].gas_mixture.composition
            connections[mapping_connections[i_node - 1][inlet_index - 1]].gas_mixture.composition = gas_composition
            connection = connections[mapping_connections[i_node - 1][inlet_index - 1]]
            if composition_tracking:
                gas_composition, connection = gas_mixture_transportation(connection, time_step=time_step)
            else:
                gas_composition = nodes[inlet_index].gas_mixture.composition
            connection.gas_mixture.composition = gas_composition
            inflow_rate = flow_matrix[i_node-1][inlet_index-1]
            inflow_temperature = connections[mapping_connections[i_node-1][inlet_index-1]].calc_pipe_outlet_temp()

            # Sum up flow rate * temperature
            total_inflow_temperature_times_flow_rate += inflow_rate * inflow_temperature
            total_inflow += inflow_rate

            # create a OrderedDict to store gas flow fractions
            gas_flow_comp = OrderedDict({gas: comp * inflow_rate for gas, comp in gas_composition.items()})
            for gas, comp in gas_flow_comp.items():
                if total_inflow_comp.get(gas) is None:
                    total_inflow_comp[gas] = comp
                else:
                    total_inflow_comp[gas] += comp

        try:
            nodal_gas_inflow_composition[i_node] = {k: v / total_inflow for k, v in total_inflow_comp.items()}
        except RuntimeWarning:
            print(total_inflow_comp)
            print(total_inflow)

        if total_inflow != .0:
            nodal_gas_inflow_temperature[i_node] = total_inflow_temperature_times_flow_rate / total_inflow
        else:
            nodal_gas_inflow_temperature[i_node] = np.nan

    return nodal_gas_inflow_composition, nodal_gas_inflow_temperature


def calculate_flow_matrix(network, pressure_bar):
    connections = network.connections
    nodes = network.nodes
    n_nodes = len(nodes)
    flow_mat = np.zeros((n_nodes, n_nodes), dtype=float)

    pressure_index = 0
    for node in nodes.values():
        if node.index not in network.non_junction_nodes:
            node.pressure = pressure_bar[pressure_index] * 1e5
            pressure_index += 1

    for connection in connections.values():
        i = connection.inlet_index - 1
        j = connection.outlet_index - 1
        connection.inlet = nodes[i+1]
        connection.outlet = nodes[j+1]

        flow_direction = connection.determine_flow_direction()

        p1 = nodes[i+1].pressure
        p2 = nodes[j+1].pressure

        slope_correction = connection.calc_pipe_slope_correction()
        temp = connection.calculate_coefficient_for_iteration()

        flow_rate = flow_direction * abs(p1 ** 2 - p2 ** 2 - slope_correction) ** (1 / 2) * temp

        flow_mat[i][j] = - flow_rate
        flow_mat[j][i] = flow_rate

    return flow_mat


def calculate_flow_vector(network, pressure_bar, target_flow):
    flow_matrix = calculate_flow_matrix(network, pressure_bar)
    n_nodes = len(network.nodes.values())
    nodal_flow = np.dot(flow_matrix, np.ones(n_nodes))
    nodal_flow = [nodal_flow[i] for i in range(len(nodal_flow)) if i + 1 not in network.non_junction_nodes]
    delta_flow = target_flow - nodal_flow

    # delta_flow = [delta_flow[i] for i in range(len(delta_flow)) if i + 1 not in network.non_junction_nodes]
    return delta_flow


def plot_network_demand_distribution(network):
    nodes = network.nodes.values()
    node_demand = [n.volumetric_flow for n in nodes if n.volumetric_flow is not None]
    sns.histplot(data=node_demand, stat="probability")
    plt.xlim((min(node_demand)-10, max(node_demand) + 10))
    plt.xlabel("Nodal volumetric flow demand [sm^3/s]")
    plt.show()
    return None


def check_square_matrix(a):
    return a.shape[0] == a.shape[1]

def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


def check_all_off_diagonal_elements(a, criterion):
    res = True

    if check_square_matrix(a):
        pass
    else:
        print("Matrix is not a square matrix!")

    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if i != j:
                if criterion == "zero":
                    res = (a[i][j] == 0)
                elif criterion == "positive":
                    res = (a[i][j] > 0)
                elif criterion == "non-negative":
                    res = (a[i][j] >= 0)
                elif criterion == "negative":
                    res = (a[i][j] < 0)
                elif criterion == "non-positive":
                    res = (a[i][j] <= 0)
                else:
                    print("Check the given criterion!")
                    return False
                if res == False:
                    return False
    return res


def check_all_diagonal_elements(a, criterion):
    res = True

    if check_square_matrix(a):
        pass
    else:
        print("Matrix is not a square matrix!")

    if criterion == "zero":
        res = (np.diagonal(a) == 0).all()
    elif criterion == "positive":
        res = (np.diagonal(a) > 0).all()
    elif criterion == "non-negative":
        res = (np.diagonal(a) >= 0).all()
    elif criterion == "negative":
        res = (np.diagonal(a) < 0).all()
    elif criterion == "non-positive":
        res = (np.diagonal(a) <= 0).all()
    else:
        print("Check the given criterion!")
        return False

    return res
