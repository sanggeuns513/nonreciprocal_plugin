// Nonreciprocal plugin for HOOMD-blue v4 and v5

#include "hoomd/HOOMDMath.h"
#include "hoomd/ForceCompute.h"
#include "hoomd/md/NeighborList.h"
#include <pybind11/pybind11.h>

#ifndef __NONRECIPROCAL_PLUGIN_H__
#define __NONRECIPROCAL_PLUGIN_H__
namespace hoomd
{
    namespace md
    {
        class PYBIND11_EXPORT NonReciprocalForce : public ForceCompute
        {
        public:
            NonReciprocalForce(std::shared_ptr<SystemDefinition> sysdef,
                               std::shared_ptr<NeighborList> nlist,
                               Scalar chi_par, Scalar chi_per, Scalar sigma, Scalar eps, Scalar r_cut);

            ~NonReciprocalForce();

        protected:
            //! compute the forces
            virtual void computeForces(uint64_t timestep);

            // Neighbor list
            std::shared_ptr<NeighborList> nlist;
            // Lennard Jones parameter: lj1 = 4.0 * epsilon * pow(sigma, 12.0), lj2 = 4.0 * epsilon * pow(sigma, 6.0)
            Scalar lj1;
            Scalar lj2;
            // Chi_parallel parameter
            Scalar chi_par;
            // Chi_perpendicular parameter
            Scalar chi_per;
            // Square of cutoff distance
            Scalar rcutsq;
        };
    } // namespace md
} // namespace hoomd
#endif