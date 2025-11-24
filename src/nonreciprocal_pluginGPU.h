// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Modification:
// Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

#include "nonreciprocal_plugin.h"
#include "hoomd/Autotuner.h"

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

#ifndef __NONRECIPROCAL_PLUGIN_GPU_H__
#define __NONRECIPROCAL_PLUGIN_GPU_H__

namespace hoomd
{
    namespace md
    {
        class PYBIND11_EXPORT NonReciprocalForceGPU : public NonReciprocalForce
        {
        public:
            NonReciprocalForceGPU(std::shared_ptr<SystemDefinition> sysdef,
                                  std::shared_ptr<NeighborList> nlist,
                                  Scalar chi_par, Scalar chi_per, Scalar sigma, Scalar eps, Scalar r_cut);
            
            ~NonReciprocalForceGPU();

        protected:
            std::shared_ptr<Autotuner<1>> m_tuner;
            virtual void computeForces(uint64_t timestep);
        };

        namespace detail
        {
            void export_NonReciprocalForceGPU(pybind11::module &m);
        }
    }

}
#endif