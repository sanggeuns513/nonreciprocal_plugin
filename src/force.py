# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

# Modification:
# Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

"""Non reciprocal force cacluation"""

from hoomd.nonreciprocal_plugin import _nonreciprocal_plugin
from hoomd.md import force
import hoomd

class NonReciprocalForce(force.Force):
    def __init__(self, simulation, neighbor_list, chi_par, chi_per, sigma, eps, rcut):
        super().__init__()
        self._nl = neighbor_list
        #if len(self._nl.r_cut.param_dict) == 0:
        #    raise RuntimeError("In neighborlist, you need to set a default_r_cut value.")
        self._chi_par = chi_par
        self._chi_per = chi_per
        self._sigma = sigma
        self._eps = eps
        self._rcut = rcut
        self._nl._attach(simulation)

    def _attach_hook(self):
        sim = self._simulation
        class_name = _nonreciprocal_plugin
        if isinstance(self._simulation.device, hoomd.device.CPU):
            self._cpp_obj = class_name.NonReciprocalForce(sim.state._cpp_sys_def, self._nl._cpp_obj, self._chi_par, self._chi_per, self._sigma, self._eps, self._rcut)
        else:
            self._cpp_obj = class_name.NonReciprocalForceGPU(sim.state._cpp_sys_def, self._nl._cpp_obj, self._chi_par, self._chi_per, self._sigma, self._eps, self._rcut)