# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

# Modification:
# Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

"""Non-reciprocal force plugin for HOOMD-blue."""
# This plugin implements a pairwise Lennard–Jones (LJ) force plus an optional non-reciprocal modification 
# that makes the force applied to particle i from j differ from the force applied to j from i for certain type pairs and separations. 
# The main work is done in computeForces(timestep) in nonreciprocal_plugin.cc and 
# gpu_compute_nonreciprocal_force_kernel(...) in nonreciprocal_pluginGPU.cu

from . import version
from hoomd.nonreciprocal_plugin import force
