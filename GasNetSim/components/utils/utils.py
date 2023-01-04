#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 3/28/22, 12:44 AM
#     Last change by yifei
#    *****************************************************************************
import math
from pyparsing import col
from collections import OrderedDict
from scipy import sparse
import numpy as np


def create_connection_matrix(n_nodes: int, components: dict, component_type: int, sparse_matrix: bool = False):
    row_ind = list()
    col_ind = list()
    data = list()

    if not sparse_matrix:
        cnx = np.zeros((n_nodes, n_nodes))

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


def forward_substitution(L, b):
    n_row = len(b)

    z = np.zeros(n_row)

    for row in range(n_row):
        res = b[row]
        for i in range(row - 1, -1):
            res = res - L[row, i] * z[i]
        z[row] = res
    return z


def backward_substitution(U, z):
    n_row = len(z)
    x = np.zeros(n_row)

    return x


def levenberg_marquardt_damping_factor(m, s, b):
    return 10 ** (m * math.log10(s + b))


def delete_matrix_rows_and_columns(matrix, to_remove):
    new_matrix = matrix

    new_matrix = np.delete(new_matrix, to_remove, 0)  # delete rows
    new_matrix = np.delete(new_matrix, to_remove, 1)  # delete columns

    return new_matrix


def gas_mixture_transportation(connection, inlet_index, outlet_index,
                               time_step, composition, location):
    length = connection.length
    velocity = connection.flow_velocity
    outlet_composition = connection.outlet_composition

    while location[0] >= length:  # if the head of a batch reached the end of the pipeline
        outlet_composition = composition[0]
        composition, location = composition[1:], location[1:]

    location = np.append(location, 0)
    composition = np.append(composition, connection.inlet.gas_mixture.composition)
    location += time_step * velocity
    return outlet_composition, location, composition


def calculate_nodal_inflow_states(nodes, connections, mapping_connections, flow_matrix, composition_tracking=False):
    nodal_total_inflow = np.sum(np.where(flow_matrix > 0, flow_matrix, 0), axis=1)

    nodal_gas_inflow_composition = dict()
    nodal_gas_inflow_temperature = dict()

    for i_node, node in nodes.items():  # iterate over all nodes
        inflow_from_node = np.where(flow_matrix[i_node-1] > 0)[0]  # find the supplying nodes
        if len(inflow_from_node) == 0:
            pass
        else:
            inflow_from_node += 1

        total_inflow_comp = dict()
        total_inflow = nodal_total_inflow[i_node-1]
        total_inflow_temperature_times_flow_rate = 0

        for inlet_index in inflow_from_node:
            connection = connections[mapping_connections[i_node - 1][inlet_index - 1]]
            if composition_tracking:
                gas_composition, location, composition = connection.gas_mixture_transportation()
            else:
                gas_composition = nodes[inlet_index].gas_mixture.composition
            connection.gas_mixture.composition = gas_composition
            inflow_rate = flow_matrix[i_node-1][inlet_index-1]
            inflow_temperature = connections[mapping_connections[i_node-1][inlet_index-1]].calc_pipe_outlet_temp()

            # Sum up flow rate * temperature
            total_inflow_temperature_times_flow_rate += inflow_rate * inflow_temperature

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

        nodal_gas_inflow_temperature[i_node] = total_inflow_temperature_times_flow_rate / total_inflow

    return nodal_gas_inflow_composition, nodal_gas_inflow_temperature
