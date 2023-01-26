#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 7/28/22, 4:19 PM
#     Last change by yifei
#    *****************************************************************************
import pandas as pd
from pathlib import Path
import logging
import copy

from ..components.network import Network


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.WARNING)


def read_profiles(file):
    profiles = pd.read_csv(Path(file))
    logger.info(f'Reading profiles from {file}.')
    return profiles


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█'):
    """
    Call in a loop to create terminal progress bar.
    the code is mentioned in : https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    # logger.info('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix))
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end="")
    # Print New Line on Complete
    if iteration == total:
        print("\n")


def check_profiles(profiles):
    if profiles['time'].dtype == int:
        print()
    elif profiles['time'].dtype == pd.Timestamp:
        print()
        # df['time'] = df['time'].apply(lambda x: pd.Timestamp.now() + pd.Timedelta(seconds=x))


def run_snapshot(network):
    network = network.simulation()
    return network


def not_converged(time_step, ts_variables):
    logger.error(f'CalculationNotConverged at time step {time_step}.')
    if not ts_variables["continue_on_divergence"]:
        raise ts_variables['errors'][0]


def update_network_topology(network):
    full_network = copy.deepcopy(network)  # make a copy of the fully connected network
    network_nodes = full_network.nodes
    network_pipes = full_network.pipelines
    network_resistances = full_network.resistances

    # Check which nodes need to be removed
    removed_nodes = dict()
    remaining_nodes = dict()
    for i, node in list(network_nodes.items()):
        if node.flow == 0:
            removed_nodes[i] = node
        else:
            remaining_nodes[i - len(removed_nodes)] = node

    # Check which pipelines need to be removed
    removed_pipes = dict()
    remaining_pipes = dict()
    for i, pipe in list(network_resistances.items()):
        if (pipe.inlet_index in removed_nodes.keys()) or (pipe.outlet_index in removed_nodes.keys()):
            pipe.valve = 1
        else:
            pipe.valve = 0
        if pipe.valve == 1:
            removed_pipes[i] = pipe
        else:
            pipe.inlet_index = list(remaining_nodes.keys())[list(remaining_nodes.values()).index(pipe.inlet)]
            pipe.outlet_index = list(remaining_nodes.keys())[list(remaining_nodes.values()).index(pipe.outlet)]
            remaining_pipes[i - len(removed_pipes)] = pipe

    return Network(nodes=remaining_nodes, pipelines=None, resistances=remaining_pipes)


def run_time_series(network, file=None):
    # create a copy of the input network
    full_network = copy.deepcopy(network)

    # read profile
    if file is not None:
        profiles = read_profiles(file)
        time_steps = profiles.index
    else:
        time_steps = range(5)  # test with 5 fictitious time steps

    # create error log to record the time step indices where error occurs
    error_log = list()

    for t in time_steps:
        for i in full_network.nodes.keys():
            if i in full_network.reference_nodes:
                pass
            else:
                try:
                    full_network.nodes[i].demand = profiles[str(i)][t]
                    full_network.nodes[i].demand_type = 'energy'
                except KeyError:
                    pass
        simplified_network = update_network_topology(full_network)
        try:
            network = run_snapshot(simplified_network)
        except RuntimeError:
            error_log.append([simplified_network, profiles.iloc[t]])

    return network