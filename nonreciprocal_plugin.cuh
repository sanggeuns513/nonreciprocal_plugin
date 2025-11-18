#include "hip/hip_runtime.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/ParticleData.cuh"

#ifndef __NONRECIPROCAL_PLUGIN_GPU_CUH__
#define __NONRECIPROCAL_PLUGIN_GPU_CUH__
namespace hoomd
{
    namespace md
    {
        namespace kernel
        {
            hipError_t gpu_compute_nonreciprocal_force(const unsigned int group_size, 
            Scalar4* d_force, const Scalar4* d_pos, const size_t *d_head_list, const unsigned int *d_n_neigh, const unsigned int *d_nlist_array, 
                const Scalar rcutsq, const Scalar lj1, const Scalar lj2, const Scalar chi_par, const BoxDim box,
                const unsigned int N, unsigned int block_size);
        }
    }
}
#endif