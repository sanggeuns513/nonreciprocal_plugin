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
                    << "Creating a ConstantForceComputeGPU with no GPU in the execution configuration"
                    << std::endl;
                throw std::runtime_error("Error initializing ConstantForceComputeGPU");
            }

            // initialize autotuner
            m_tuner.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(m_exec_conf)},
                                           m_exec_conf,
                                           "nonreciprocal_plugin"));
            m_autotuners.push_back(m_tuner);
        }

        void NonReciprocalForceGPU::computeForces(uint64_t timestep)
        {
            // access particle data
            ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
            ArrayHandle<Scalar4> d_force(m_force, access_location::device, access_mode::overwrite);
            // ArrayHandle<Scalar4> d_torque(m_torque, access_location::device, access_mode::overwrite)

            nlist->compute(timestep);
            ArrayHandle<size_t> d_head_list(nlist->getHeadList(), access_location::device, access_mode::read);
            ArrayHandle<unsigned int> d_nlist_array(nlist->getNListArray(), access_location::device, access_mode::read);
            ArrayHandle<unsigned int> d_n_neigh(nlist->getNNeighArray(), access_location::device, access_mode::read);

            unsigned int group_size = m_group->getNumMembers();
            unsigned int N = m_pdata->getN();
            BoxDim box = m_pdata->getGlobalBox();

            // Compute the forces on the GPU
            m_tuner->begin();
            // TODO: GPU compute force under namespace kernel
            kernel::gpu_compute_nonreciprocal_force(group_size,
                                                    d_force.data, d_pos.data, d_head_list.data, d_n_neigh.data, d_nlist_array.data,
                                                    rcutsq, lj1, lj2, chi_par, box, N, m_tuner->getParam()[0]);

            if (m_exec_conf->isCUDAErrorCheckingEnabled())
                CHECK_CUDA_ERROR();

            m_tuner->end();
        }
        namespace detail
        {
            PYBIND11_MODULE(_nonreciprocal_plugin, m)
            {
                pybind11::class_<NonReciprocalForceGPU, NonReciprocalForce, std::shared_ptr<NonReciprocalForceGPU>>(m, "NonReciprocalForceGPU")
                    .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<NeighborList>, Scalar, Scalar, Scalar, Scalar, Scalar>());
            }
        }
    }
}