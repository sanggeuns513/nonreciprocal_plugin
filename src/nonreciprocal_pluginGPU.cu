// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Modification:
// Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

#include "nonreciprocal_plugin.cuh"
#include "hoomd/RNGIdentifiers.h"
#include "hoomd/RandomNumbers.h"
#include "hoomd/TextureTools.h"

namespace hoomd
{
    namespace md
    {
        namespace kernel
        {
            Scalar *d_chi_par_flat;
            Scalar *d_chi_per_flat;
            __global__ void gpu_compute_nonreciprocal_force_kernel(Scalar4 *d_force, const Scalar4 *__restrict__ d_pos, const size_t *__restrict__ d_head_list,
                                                                   const unsigned int *__restrict__ d_n_neigh, const unsigned int *__restrict__ d_nlist_array,
                                                                   const Scalar rcutsq, const Scalar lj1, const Scalar lj2,
                                                                   Scalar *__restrict__ chi_par_flat, Scalar *__restrict__ chi_per_flat,
                                                                   const BoxDim box, const unsigned int N, const unsigned int ntypes, const Scalar drecsq)
            {
                unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
                unsigned int stride = gridDim.x * blockDim.x;
                for (unsigned int i = idx; i < N; i += stride)
                {
                    const Scalar4 pi4 = __ldg(d_pos + i);
                    const Scalar3 pi = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
                    const unsigned int type_i = __scalar_as_int(d_pos[i].w);

                    Scalar3 fi = make_scalar3(0, 0, 0);

                    const size_t myHead = d_head_list[i];
                    const unsigned int size = (unsigned int)d_n_neigh[i];
                    for (unsigned int k = 0; k < size; k++)
                    {
                        const unsigned int j = __ldg(d_nlist_array + myHead + k);
                        const Scalar4 pj4 = __ldg(d_pos + j);
                        const unsigned int type_j = __scalar_as_int(d_pos[j].w);
                        Scalar3 pj = make_scalar3(d_pos[j].x, d_pos[j].y, d_pos[j].z);
                        Scalar3 dx = pi - pj;
                        dx = box.minImage(dx);
                        Scalar rsq = dot(dx, dx);
                        if (rsq >= rcutsq)
                            continue;
                        Scalar r2inv = Scalar(1.0) / rsq;
                        Scalar r6inv = r2inv * r2inv * r2inv;
                        Scalar f_div_r = r2inv * r6inv * (Scalar(12.0) * lj1 * r6inv - Scalar(6.0) * lj2);

                        Scalar3 f_ij = f_div_r * dx;
                        if (rsq > drecsq)
                        {
                            unsigned int idx = type_i * ntypes + type_j;
                            Scalar chi_p = __ldg(chi_par_flat + idx);
                            Scalar chi_t = __ldg(chi_per_flat + idx);

                            fi.x += f_ij.x * (Scalar(1.0) + chi_p) + f_ij.y * chi_t;
                            fi.y += f_ij.y * (Scalar(1.0) + chi_p) - f_ij.x * chi_t;
                            fi.z += f_ij.z * (Scalar(1.0) + chi_p);
                        }
                        else
                        {
                            fi += f_ij;
                        }
                    }
                    atomicAdd(&d_force[i].x, fi.x);
                    atomicAdd(&d_force[i].y, fi.y);
                    atomicAdd(&d_force[i].z, fi.z);
                }
            }

            hipError_t gpu_allocate_chi(const std::vector<Scalar> chi_par_flat, const std::vector<Scalar> chi_per_flat, const unsigned int ntypes)
            {
                hipError_t err = hipMalloc(&d_chi_par_flat, sizeof(Scalar) * ntypes * ntypes);
                if (err != hipSuccess) return err;
                err = hipMalloc(&d_chi_per_flat, sizeof(Scalar) * ntypes * ntypes);
                if (err != hipSuccess)
                {
                    hipFree(d_chi_par_flat);
                    return err;
                }
                err = hipMemcpy(d_chi_par_flat, chi_par_flat.data(), sizeof(Scalar) * ntypes * ntypes, hipMemcpyHostToDevice);
                if (err != hipSuccess)
                {
                    hipFree(d_chi_par_flat);
                    hipFree(d_chi_per_flat);
                    return err;
                }
                err = hipMemcpy(d_chi_per_flat, chi_per_flat.data(), sizeof(Scalar) * ntypes * ntypes, hipMemcpyHostToDevice);
                if (err != hipSuccess)
                {
                    hipFree(d_chi_par_flat);
                    hipFree(d_chi_per_flat);
                    return err;
                }

                return hipSuccess;
            }

            hipError_t gpu_compute_nonreciprocal_force(Scalar4 *d_force, const Scalar4 *d_pos, const size_t *d_head_list, const unsigned int *d_n_neigh,
                                                       const unsigned int *d_nlist_array, const Scalar rcutsq, const Scalar lj1, const Scalar lj2,
                                                       const BoxDim box, const unsigned int N, unsigned int block_size, const unsigned int ntypes, Scalar drecsq)
            {
                // setup the grid to run the kernel
                dim3 grid(N / block_size + 1, 1, 1);
                dim3 threads(block_size, 1, 1);

                
                hipMemset(d_force, 0, sizeof(Scalar4) * N);
                // hipMemset(d_torque, 0, sizeof(Scalar4) * N);
                hipLaunchKernelGGL((gpu_compute_nonreciprocal_force_kernel),
                                   dim3(grid), dim3(threads), 0, 0,
                                   d_force, d_pos, d_head_list, d_n_neigh, d_nlist_array, rcutsq, lj1, lj2, d_chi_par_flat, d_chi_per_flat, box, N, ntypes, drecsq);

                return hipSuccess;
            }

            hipError_t gpu_deallocate_chi()
            {
                hipFree(d_chi_par_flat);
                hipFree(d_chi_per_flat);
                return hipSuccess;
            }
        }
    }
}