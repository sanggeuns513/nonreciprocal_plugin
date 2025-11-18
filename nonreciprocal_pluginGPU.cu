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
            __global__ void gpu_compute_nonreciprocal_force_kernel(const unsigned int group_size, 
                Scalar4* d_force, const Scalar4* d_pos, const size_t *d_head_list, const unsigned int *d_n_neigh, const unsigned int *d_nlist_array, 
                const Scalar rcutsq, const Scalar lj1, const Scalar lj2, const Scalar chi_par, const BoxDim box, const unsinged int N)
            {
                unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
                if (group_idx >= group_size) return;
                Scalar3 ri = make_scalar3(d_pos[group_idx].x, d_pos[group_idx].y, d_pos[group_idx].z)
                Scalar3 Fi = make_scalar3(0, 0, 0);
                unsigned int type_i = static_cast<unsigned int>(d_pos[group_idx].w);
                for (unsigned int k = 0; k < d_n_neigh[group_idx]; k++)
                {
                    unsigned int j = d_nlist_array.data[d_head_list[group_idx] + k];
                    unsigned int type_j = static_cast<unsigned int>(d_pos[j].w);
                    Scalar3 rj = make_scalar3(d_pos[j].x, d_pos[j].y, d_pos[j].z);
                    Scalar3 rij = box.minImage(rj - ri);
                    Scalar rsq = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
                    if (rsq >= rcutsq) continue;
                    Scalar r2inv = Scalar(1.0) / rsq;
                    Scalar r6inv = r2inv * r2inv * r2inv;
                    Scalar f_div_r = r2inv * r6inv * (Scalar(12.0) * lj1 * r6inv - Scalar(6.0) * lj2);
                    Scalar3 f_ij = -f_div_r * rij;
                    Scalar3 f_ji = -f_ij;
                    if (type_i == 0 && type_j == 1)
                    {
                        Fi.x += f_ij.x * (Scalar(1.0) + chi_par) + f_ij.y * chi_per;
                        Fi.y += f_ij.y * (Scalar(1.0) + chi_par) - f_ij.x * chi_per;
                        Fi.z += f_ij.z * (Scalar(1.0) + chi_par);

                        atomicAdd(&d_force[j].x, f_ji.x * (Scalar(1.0) - chi_par) - f_ji.y * chi_per);
                        atomicAdd(&d_force[j].y, f_ji.y * (Scalar(1.0) - chi_par) + f_ji.x * chi_per);
                        atomicAdd(&d_force[j].z, f_ji.z * (Scalar(1.0) - chi_par));
                    }
                    else if (type_i == 1 && type_j == 0)
                    {
                        Fi.x += f_ij.x * (1.0 - chi_par) - f_ij.y * chi_per;
                        Fi.y += f_ij.y * (1.0 - chi_par) + f_ij.x * chi_per;
                        Fi.z += f_ij.z * (1.0 - chi_par);

                        atomicAdd(&d_force[j].x, f_ji.x * (Scalar(1.0) + chi_par) + f_ji.y * chi_per);
                        atomicAdd(&d_force[j].y, f_ji.y * (Scalar(1.0) + chi_par) - f_ji.x * chi_per);
                        atomicAdd(&d_force[j].z, f_ji.z * (Scalar(1.0) + chi_par));
                    }
                    else
                    {
                        Fi += f_ij;
                        atomicAdd(&d_force[j].x, f_ji.x);
                        atomicAdd(&d_force[j].y, f_ji.y);
                        atomicAdd(&d_force[j].z, f_ji.z);
                    }
                }
                atomicAdd(&d_force[group_idx].x, Fi.x);
                atomicAdd(&d_force[group_idx].y, Fi.y);
                atomicAdd(&d_force[group_idx].z, Fi.z);

            }

            hipError_t gpu_compute_nonreciprocal_force(const unsigned int group_size, 
            Scalar4* d_force, const Scalar4* d_pos, const size_t *d_head_list, const unsigned int *d_n_neigh, const unsigned int *d_nlist_array, 
                const Scalar rcutsq, const Scalar lj1, const Scalar lj2, const Scalar chi_par, const BoxDim box,
                const unsigned int N, unsigned int block_size)
            {
                // setup the grid to run the kernel
                dim3 grid(group_size / block_size + 1, 1, 1);
                dim3 threads(block_size, 1, 1);

                hipMemset(d_force, 0, sizeof(Scalar4) * N);
                // hipMemset(d_torque, 0, sizeof(Scalar4) * N);
                hipLaunchKernelGGL((gpu_compute_nonreciprocal_force_kernel),
                                    dim3(grid), dim3(threads), 0, 0,
                                    group_size, d_force, d_pos, d_head_list, d_n_neigh, d_nlist_array, rcutsq, lj1, lj2, chi_par, box, N);
                return hipSuccess;
            }
        }
    }
}