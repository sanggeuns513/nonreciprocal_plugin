// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Modification:
// Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

#include "nonreciprocal_plugin.h"
#include "nonreciprocal_pluginGPU.h"
#include <pybind11/pybind11.h>

namespace hoomd
{
    namespace md
    {
        namespace detail
        {
            void export_NonReciprocalForce(pybind11::module &m);
            void export_NonReciprocalForceGPU(pybind11::module &m);
        }
    }
}

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
    create the md python module and define the exports here.
*/

using namespace hoomd;
using namespace hoomd::md;
using namespace hoomd::md::detail;

PYBIND11_MODULE(_nonreciprocal_plugin, m)
{
    export_NonReciprocalForce(m);
#ifdef ENABLE_HIP
    export_NonReciprocalForceGPU(m);
#endif
}