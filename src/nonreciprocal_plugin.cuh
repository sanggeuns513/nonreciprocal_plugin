// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Modification:
// Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

#include "hip/hip_runtime.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/ParticleData.cuh"
#include <vector>

#ifndef __NONRECIPROCAL_PLUGIN_GPU_CUH__
#define __NONRECIPROCAL_PLUGIN_GPU_CUH__
namespace hoomd
{
    namespace md
    {
        namespace kernel
        {
            extern Scalar *d_chi_par_flat;
            extern Scalar *d_chi_per_flat;
            hipError_t gpu_allocate_chi(const std::vector<Scalar> chi_par_flat, const std::vector<Scalar> chi_per_flat, unsigned int ntypes);
            hipError_t gpu_deallocate_chi();
            hipError_t gpu_compute_nonreciprocal_force(Scalar4 *d_force, const Scalar4 *d_pos, const size_t *d_head_list, const unsigned int *d_n_neigh,
                                                       const unsigned int *d_nlist_array, const Scalar rcutsq, const Scalar lj1, const Scalar lj2,
                                                       const BoxDim box, const unsigned int N, unsigned int block_size, unsigned int ntypes, const Scalar drecsq);
        }
    }
}
#endif