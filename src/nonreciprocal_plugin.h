// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Modification:
// Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

#include "hoomd/HOOMDMath.h"
#include "hoomd/ForceCompute.h"
#include "hoomd/md/NeighborList.h"
#include <pybind11/pybind11.h>
#include <vector>

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

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
            std::vector<Scalar> chi_par_flat;
            // Chi_perpendicular parameter
            std::vector<Scalar> chi_per_flat;
            // Square of cutoff distance
            Scalar rcutsq;
            // This is just for benchmark. Will delete this variable.
            Scalar drecsq;
        };
        namespace detail
        {
            void export_NonReciprocalForce(pybind11::module& m);
        }
    } // namespace md
} // namespace hoomd
#endif