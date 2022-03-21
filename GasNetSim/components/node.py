#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 1/17/22, 11:20 AM
#    Last change by yifei
#   *****************************************************************************
from .gas_mixture.typical_mixture_composition import *
from ..components.gas_mixture.thermo.thermo import Mixture
from ..components.gas_mixture.heating_value import *


class Node:
    """
    Class to formulate gas transmission network nodes.
    """

    def __init__(self, node_index, flow=None, pressure_pa=None, temperature=288.15, altitude=0, gas_composition=None,
                 node_type='demand', flow_type='volumetric', volumetric_flow=0, energy_flow=0):
        """
        Initial method
        :param node_index: Node index
        :param flow: Gas volumetric flow [sm3/s] or energy flow [MJ/s]
        :param pressure_pa: Gas nodal pressure [Pa]
        :param temperature: Gas nodal temperature [K]
        :param altitude: Elevation of the network node [m]
        """
        self.index = node_index
        if gas_composition is not None:
            self.gas_composition = gas_composition
        else:
            self.gas_composition = NATURAL_GAS_gri30
        self.pressure = pressure_pa
        if pressure_pa is not None:
            self.pressure_bar = pressure_pa / 101325
        if temperature is not None:
            self.temperature = temperature
        else:
            self.temperature = 288.15
        if altitude is not None:
            self.altitude = altitude
        else:
            self.altitude = 0
        if node_type is not None:
            self.node_type = node_type
        else:
            self.node_type = 'demand'
        # flow type
        if flow_type is not None:
            self.flow_type = flow_type
        else:
            self.flow_type = 'volumetric'
        if flow is None and pressure_pa is None:
            raise InitializationError("Either pressure or flow should be known.")
        try:
            self.gas_mixture = Mixture(zs=self.gas_composition,
                                       T=self.temperature_K,
                                       P=self.pressure_pa)
        except (TypeError, AttributeError):
            # If pressure or temperature is missing for some nodes
            self.gas_mixture = Mixture(zs=self.gas_composition,
                                       T=288.15,
                                       P=50 * 101325)

        self.flow = flow
        if self.flow_type == 'volumetric':
            self.volumetric_flow = flow
            try:
                self.convert_volumetric_to_energy_flow()
            except TypeError:
                self.energy_flow = None
        elif self.flow_type == 'energy':
            self.energy_flow = flow
            try:
                self.convert_energy_to_volumetric_flow()
            except TypeError:
                self.volumetric_flow = None
        else:
            raise AttributeError(f'Unknown flow type {flow_type}!')

    def update_gas_mixture(self):
        try:
            self.gas_mixture = Mixture(zs=self.get_mole_fraction(),
                                       T=self.temperature,
                                       P=self.pressure_pa)
        except (TypeError, AttributeError):
            self.gas_mixture = Mixture(zs=NATURAL_GAS_gri30,
                                       T=288.15,
                                       P=50 * 101325)

    def get_mole_fraction(self):
        """
        Get mole fraction of the gas composition at node
        :return: Gas mole fraction
        """
        mole_fraction = dict()
        for i in range(len(self.gas_mixture.components)):
            gas = self.gas_mixture.components[i]
            mole_fraction[gas] = self.gas_mixture.zs[i]
        return mole_fraction

    def convert_energy_to_volumetric_flow(self):
        """
        Convert energy flow rate (MW) into volumetric flow rate (sm^3/s)
        :return:
        """
        HHV = calc_heating_value(self.gas_mixture)
        gas_comp = self.get_mole_fraction()
        self.volumetric_flow = self.energy_flow / HHV * 1e6 / Mixture(zs=gas_comp, T=288.15, P=101325).rho

    def convert_volumetric_to_energy_flow(self):
        """
        Convert volumetric flow rate (sm^3/s) into energy flow rate (MW)
        :return:
        """
        HHV = calc_heating_value(self.gas_mixture)
        gas_comp = self.get_mole_fraction()
        self.energy_flow = self.volumetric_flow * HHV / 1e6 * Mixture(zs=gas_comp, T=288.15, P=101325).rho


if __name__ == "__main__":
    from gas_mixture.typical_mixture_composition import *
    # from ..components.gas_mixture.thermo.thermo import Mixture
    # from ..components.gas_mixture.heating_value import *
    Node(flow=None, pressure_pa=None, temperature=300)