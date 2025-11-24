// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Modification:
// Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

#include "nonreciprocal_plugin.h"

namespace hoomd
{
    namespace md
    {
        NonReciprocalForce::NonReciprocalForce(std::shared_ptr<SystemDefinition> sysdef,
                                               std::shared_ptr<NeighborList> nlist,
                                               Scalar chi_par, Scalar chi_per, Scalar sigma, Scalar eps, Scalar r_cut)
            : ForceCompute(sysdef), nlist(nlist), rcutsq(r_cut * r_cut)
        {
            // lj1 = 4.0 * epsilon * pow(sigma, 12.0)
            // lj2 = 4.0 * epsilon * pow(sigma, 6.0)
            lj1 = 4.0 * eps * sigma * sigma * sigma * sigma * sigma * sigma * sigma * sigma * sigma * sigma * sigma * sigma;
            lj2 = 4.0 * eps * sigma * sigma * sigma * sigma * sigma * sigma;
            // chi_parallel and chi_perpendicular
            unsigned int ntypes = m_pdata->getNTypes();
            chi_par_flat.assign(ntypes, 0.0);
            chi_per_flat.assign(ntypes, 0.0);
            chi_par_flat[0 * ntypes + 1] = chi_par;
            chi_par_flat[1 * ntypes + 0] = -chi_par;
            chi_per_flat[0 * ntypes + 1] = chi_per;
            chi_per_flat[1 * ntypes + 0] = -chi_per;
            // Square of d_rec, where d_rec = pow(2.0, 1.0 / 6.0) from Y-J Chiu and A. K. Omar, J. Chem. Phys. 158, 164903 (2023)
            drecsq = 1.2599210499;
        }

        NonReciprocalForce::~NonReciprocalForce()
        {
            m_exec_conf->msg->notice(5) << "Destroying NonReciprocalForce" << std::endl;
        }

        void NonReciprocalForce::computeForces(uint64_t timestep)
        {
            // Compute the neighborlist first
            nlist->compute(timestep);
            bool third_law = nlist->getStorageMode() == NeighborList::half;
            // Making another scope to avoid acquire access error
            {
                // access particle data
                ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
                ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::overwrite);
                // ArrayHandle<Scalar4> h_torque(m_torque, access_location::host, access_mode::readwrite);

                // neighbor list data
                ArrayHandle<size_t> h_head_list(nlist->getHeadList(), access_location::host, access_mode::read);
                ArrayHandle<unsigned int> h_nlist_array(nlist->getNListArray(), access_location::host, access_mode::read);
                ArrayHandle<unsigned int> h_n_neigh(nlist->getNNeighArray(), access_location::host, access_mode::read);

                unsigned int N = m_pdata->getN();
                unsigned int ntypes = m_pdata->getNTypes();
                const BoxDim box = m_pdata->getGlobalBox();

                // initialize to zero
                for (unsigned int i = 0; i < N; i++)
                {
                    h_force.data[i] = make_scalar4(0, 0, 0, 0);
                    // h_torque.data[i] = make_scalar4(0, 0, 0, 0);
                }

                // loop over all particles
                for (unsigned int i = 0; i < N; i++)
                {
                    Scalar3 pi = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);
                    Scalar3 fi = make_scalar3(0, 0, 0);
                    unsigned int type_i = __scalar_as_int(h_pos.data[i].w);
                    const size_t myHead = h_head_list.data[i];
                    const unsigned int size = (unsigned int)h_n_neigh.data[i];
                    // loop over neighbors
                    // Note that default storage mode of neighborlist is half.
                    for (unsigned int k = 0; k < size; k++)
                    {
                        unsigned int j = h_nlist_array.data[myHead + k];
                        unsigned int type_j = __scalar_as_int(h_pos.data[j].w);
                        Scalar3 pj = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);
                        Scalar3 dx = pi - pj;
                        dx = box.minImage(dx);
                        Scalar rsq = dot(dx, dx);
                        Scalar f_div_r = 0.0;
                        // Pair potential part
                        if (rsq < rcutsq)
                        {
                            // Pair interaction forces
                            Scalar r2inv = Scalar(1.0) / rsq;
                            Scalar r6inv = r2inv * r2inv * r2inv;
                            // f_div_r: magnitude of pair force divide by r
                            f_div_r = r2inv * r6inv * (Scalar(12.0) * lj1 * r6inv - Scalar(6.0) * lj2);
                        }
                        Scalar3 f_ij = f_div_r * dx;
                        Scalar3 f_ji = -f_ij;
                        // Apply non-reciprocal rule on ith particle
                        if (rsq > drecsq)
                        {
                            Scalar chi_p = chi_par_flat[ntypes * type_i + type_j];
                            Scalar chi_t = chi_per_flat[ntypes * type_i + type_j];

                            fi.x += f_ij.x * (Scalar(1.0) + chi_p) + f_ij.y * chi_t;
                            fi.y += f_ij.y * (Scalar(1.0) + chi_p) - f_ij.x * chi_t;
                            // For z component, I need to think how to modify it!
                            fi.z += f_ij.z * (Scalar(1.0) + chi_p);
                            // Force acting on jth particle
                            if (third_law)
                            {
                                h_force.data[j].x += f_ji.x * (Scalar(1.0) - chi_p) - f_ji.y * chi_t;
                                h_force.data[j].y += f_ji.y * (Scalar(1.0) - chi_p) + f_ji.x * chi_t;
                                h_force.data[j].z += f_ji.z * (Scalar(1.0) - chi_p);
                            }
                        }
                        else
                        {
                            fi += f_ij;
                            if (third_law)
                            {
                                h_force.data[j].x += f_ji.x;
                                h_force.data[j].y += f_ji.y;
                                h_force.data[j].z += f_ji.z;
                            }
                        }
                    }

                    // write force
                    h_force.data[i].x += fi.x;
                    h_force.data[i].y += fi.y;
                    h_force.data[i].z += fi.z;
                }
            } // Another scope to avoid acquire problem
        }
        namespace detail
        {
            void export_NonReciprocalForce(pybind11::module &m)
            {
                pybind11::class_<NonReciprocalForce, ForceCompute, std::shared_ptr<NonReciprocalForce>>(m, "NonReciprocalForce")
                    .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<NeighborList>, Scalar, Scalar, Scalar, Scalar, Scalar>());
            }
        }
    } // namespace md
} // namespace hoomd