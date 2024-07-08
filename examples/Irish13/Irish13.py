#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2024.
#     Developed by Yifei Lu
#     Last change on 4/29/24, 1:10 PM
#     Last change by yifei
#    *****************************************************************************

# import sys
# import os
# sys.path.append(os.path.dirname(__file__))

import GasNetSim as gns
from pathlib import Path
from timeit import default_timer as timer


network = gns.create_network_from_csv(Path('.'))

# start = timer()
# network.simulation(use_cuda=True, tol=0.0001)
# end = timer()
#
# print(f"Simulation time using CuPy: {end - start}")


network = gns.create_network_from_csv(Path('.'))

start = timer()
network.simulation(use_cuda=False, tol=0.0001)
end = timer()

print(f"Simulation time using NumPy: {end - start}")

