#include "nonreciprocal_plugin.h"
#include "hoomd_plugin/Autotuner.h"

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

        protected:
            std::shared_ptr<Autotuner<1>> m_tuner;
            virtual void computeForces(uint64_t timestep);
        };
    }
}
#endif