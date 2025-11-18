

#include "nonreciprocal_plugin.h"

namespace hoomd
{
    namespace md
    {
        NonReciprocalForce::NonReciprocalForce(std::shared_ptr<SystemDefinition> sysdef,
                                               std::shared_ptr<NeighborList> nlist,
                                               Scalar chi_par, Scalar chi_per, Scalar sigma, Scalar eps, Scalar r_cut)
            : ForceCompute(sysdef), nlist(nlist), chi_par(chi_par), chi_per(chi_per), rcutsq(r_cut * r_cut)
        {
            // lj1 = 4.0 * epsilon * pow(sigma, 12.0)
            // lj2 = 4.0 * epsilon * pow(sigma, 6.0)
            lj1 = 4.0 * eps * sigma * sigma * sigma * sigma * sigma * sigma * sigma * sigma * sigma * sigma * sigma * sigma;
            lj2 = 4.0 * eps * sigma * sigma * sigma * sigma * sigma * sigma;
        }

        NonReciprocalForce::~NonReciprocalForce()
        {
            m_exec_conf->msg->notice(5) << "Destroying NonReciprocalForce" << std::endl;
        }

        void NonReciprocalForce::computeForces(uint64_t timestep)
        {
            // access particle data
            ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
            ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::overwrite);
            // ArrayHandle<Scalar4> h_torque(m_torque, access_location::host, access_mode::overwrite);

            // neighbor list data
            nlist->compute(timestep);
            ArrayHandle<size_t> h_head_list(nlist->getHeadList(), access_location::host, access_mode::read);
            ArrayHandle<unsigned int> h_nlist_array(nlist->getNListArray(), access_location::host, access_mode::read);
            ArrayHandle<unsigned int> h_n_neigh(nlist->getNNeighArray(), access_location::host, access_mode::read);

            unsigned int N = m_pdata->getN();
            BoxDim box = m_pdata->getGlobalBox();
            
            // initialize to zero
            for (unsigned int i = 0; i < N; i++)
            {
                h_force.data[i] = make_scalar4(0, 0, 0, 0);
                // h_torque.data[i] = make_scalar4(0, 0, 0, 0);
            }

            // loop over all particles
            for (unsigned int i = 0; i < N; i++)
            {
                Scalar3 ri = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);
                Scalar3 Fi = make_scalar3(0, 0, 0);
                unsigned int type_i = static_cast<unsigned int>(h_pos.data[i].w);
                // loop over neighbors
                for (unsigned int k = 0; k < h_n_neigh.data[i]; k++)
                {
                    unsigned int j = h_nlist_array.data[h_head_list.data[i] + k];
                    unsigned int type_j = static_cast<unsigned int>(h_pos.data[j].w);
                    Scalar3 rj = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);
                    Scalar3 rij = box.minImage(rj - ri);
                    Scalar rsq = dot(rij, rij);
                    if (rsq < rcutsq)
                    {
                        // Pair interaction forces
                        Scalar r2inv = Scalar(1.0) / rsq;
                        Scalar r6inv = r2inv * r2inv * r2inv;
                        // f_div_r: magnitude of pair force divide by r
                        Scalar f_div_r = r2inv * r6inv * (Scalar(12.0) * lj1 * r6inv - Scalar(6.0) * lj2);
                        Scalar3 f_ij = -f_div_r * rij;
                        Scalar3 f_ji = -f_ij;

                        // Apply non-reciprocal rule on ith particle
                        if (type_i == 0 && type_j == 1)
                        {
                            Fi.x += f_ij.x * (1 + chi_par) + f_ij.y * chi_per;
                            Fi.y += f_ij.y * (1 + chi_par) - f_ij.x * chi_per;
                            // For z component, I need to think how to modify it!
                            Fi.z += f_ij.z * (1 + chi_par);

                            // Force acting on jth particle
                            h_force.data[j].x += f_ji.x * (1 - chi_par) - f_ji.y * chi_per;
                            h_force.data[j].y += f_ji.y * (1 - chi_par) + f_ji.x * chi_per;
                            h_force.data[j].z += f_ji.z * (1 - chi_par);
                        }
                        else if (type_i == 1 && type_j == 0)
                        {
                            Fi.x += f_ij.x * (1 - chi_par) - f_ij.y * chi_per;
                            Fi.y += f_ij.y * (1 - chi_par) + f_ij.x * chi_per;
                            // For z component, I need to think how to modify it!
                            Fi.z += f_ij.z * (1 - chi_par);

                            // Force acting on jth particle
                            h_force.data[j].x += f_ji.x * (1 + chi_par) + f_ji.y * chi_per;
                            h_force.data[j].y += f_ji.y * (1 + chi_par) - f_ji.x * chi_per;
                            h_force.data[j].z += f_ji.z * (1 + chi_par);
                        }
                        else
                        {
                            // Apply Newton 3rd rule
                            Fi += f_ij;
                            h_force.data[j].x += f_ji.x;
                            h_force.data[j].y += f_ji.y;
                            h_force.data[j].z += f_ji.z;
                        }
                    }
                }
                // write force
                h_force.data[i].x = Fi.x;
                h_force.data[i].y = Fi.y;
                h_force.data[i].z = Fi.z;
            }
        }

        namespace detail
        {
            PYBIND11_MODULE(_nonreciprocal_plugin, m)
            {
                pybind11::class_<NonReciprocalForce, ForceCompute, std::shared_ptr<NonReciprocalForce>>(m, "NonReciprocalForce")
                    .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<NeighborList>, Scalar, Scalar, Scalar, Scalar, Scalar>());
            }
        } // namespace detail
    } // namespace md
} // namespace hoomd