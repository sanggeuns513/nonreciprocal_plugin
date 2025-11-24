// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Modification:
// Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

#include "nonreciprocal_pluginGPU.h"
#include "nonreciprocal_plugin.cuh"
#include <vector>

namespace hoomd
{
    namespace md
    {
        NonReciprocalForceGPU::NonReciprocalForceGPU(std::shared_ptr<SystemDefinition> sysdef,
                                                     std::shared_ptr<NeighborList> nlist,
                                                     Scalar chi_par, Scalar chi_per, Scalar sigma, Scalar eps, Scalar r_cut)
            : NonReciprocalForce(sysdef, nlist, chi_par, chi_per, sigma, eps, r_cut)
        {
            if (!m_exec_conf->isCUDAEnabled())
            {
                m_exec_conf->msg->error()
                    << "Creating a NonReciprocalForceGPU with no GPU in the execution configuration"
                    << std::endl;
                throw std::runtime_error("Error initializing NonReciprocalForceGPU");
            }

            // initialize autotuner
            m_tuner.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(m_exec_conf)},
                                           m_exec_conf,
                                           "nonreciprocal_plugin"));
            m_autotuners.push_back(m_tuner);

            const unsigned int ntypes = m_pdata->getNTypes();
            kernel::gpu_allocate_chi(chi_par_flat, chi_per_flat, ntypes);
        }

        NonReciprocalForceGPU::~NonReciprocalForceGPU()
        {
            m_exec_conf->msg->notice(5) << "Destroying NonReciprocalForceGPU" << std::endl;
            kernel::gpu_deallocate_chi();
        }

        void NonReciprocalForceGPU::computeForces(uint64_t timestep)
        {
            nlist->compute(timestep);
            bool third_law = nlist->getStorageMode() == NeighborList::half;
            if (third_law)
            {
                this->m_exec_conf->msg->error() << "The GPU cannot handle a half neighbor list!!";
                throw std::runtime_error("Error computing forces");
            }
            {
                // access particle data
                ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
                ArrayHandle<Scalar4> d_force(m_force, access_location::device, access_mode::readwrite);
                // ArrayHandle<Scalar4> d_torque(m_torque, access_location::device, access_mode::overwrite)

                ArrayHandle<size_t> d_head_list(nlist->getHeadList(), access_location::device, access_mode::read);
                ArrayHandle<unsigned int> d_nlist_array(nlist->getNListArray(), access_location::device, access_mode::read);
                ArrayHandle<unsigned int> d_n_neigh(nlist->getNNeighArray(), access_location::device, access_mode::read);

                const unsigned int N = m_pdata->getN();
                const BoxDim box = m_pdata->getGlobalBox();
                const unsigned int ntypes = m_pdata->getNTypes();

                // Compute the forces on the GPU
                m_tuner->begin();
                kernel::gpu_compute_nonreciprocal_force(d_force.data, d_pos.data, d_head_list.data, d_n_neigh.data, d_nlist_array.data,
                                                        rcutsq, lj1, lj2, box, N, m_tuner->getParam()[0], ntypes, drecsq);

                if (m_exec_conf->isCUDAErrorCheckingEnabled())
                    CHECK_CUDA_ERROR();

                m_tuner->end();
            }
        }
        namespace detail
        {
            void export_NonReciprocalForceGPU(pybind11::module &m)
            {
                pybind11::class_<NonReciprocalForceGPU, NonReciprocalForce, std::shared_ptr<NonReciprocalForceGPU>>(m, "NonReciprocalForceGPU")
                    .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<NeighborList>, Scalar, Scalar, Scalar, Scalar, Scalar>());
            }
        }
    }
}