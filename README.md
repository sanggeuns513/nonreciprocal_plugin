
# **NonReciprocalForce**

NonReciprocalForce is a plug-in for HOOMD-Blue, a particle simulation toolkit, that implements specialized pair potentials for poly-disperse interacting particle system. The code is largely based on HOOMD-Blue's [class hoomd.md.force.Constant](https://hoomd-blue.readthedocs.io/en/stable/hoomd/md/force/constant.html).

The plugin is ready to use and is compatible with HOOMD-Blue v4 and v5.

## **Contents** 

Files that come with this plugin:
 - CMakeLists.txt   : main CMake configuration file for the plugin
 - README.md        : This file
 - src              : Directory containing C++ and Python source codes that interacts with HOOMD-Blue


## **Installation Instructions**

Parts of the instructions were modified from the hoomd-component-template repository provided by HOOMD-Blue. See https://hoomd-blue.readthedocs.io/en/stable/components.html for other useful information.

### Step 1: **Check Requirements**
The requirements for installing the plugin is the same as standard HOOMD, except that this plugin requires Pybind11 2.13 or higher. See [HOOMD installation page, building from source](https://hoomd-blue.readthedocs.io/en/stable/building.html) for details.

### Step 2: **Install Plugin**

The process is similar to installing HOOMD. First, git clone the project:
```bash
git clone https://github.com/mandadapu-group/nonreciprocalforce
```

Next, configure your build and install.
```bash
cd nonreciprocalforce
cmake -B build -S .
cmake --build ./build
cmake --install ./build
```

In this step, CMake will try to find the usual required packages. However, it will also try to find a HOOMD installation. Check your CMake output! 

## **Using NonreciprocalForce with HOOMD's MD**

The plugin is a complement for HOOMD's MD, which means that you need to import the plugin and hoomd.md side-by-side. Everything you do will be just like running HOOMD's MD except at the part where you're appending force. For example, in the case of two A, B Lennard-Jones particles with ratio 0.5,

```python
import hoomd
import hoomd.md as md
import hoomd.nonreciprocal_plugin as nrec
import numpy as np
from random import uniform
import gsd.hoomd
import os
import datetime

simulation = hoomd.Simulation(device = hoomd.device.CPU(notice_level = 2), seed = 0)

# Set up the parameter for molecular dynamics
delta_t = 0.001
totalsteps = 10/delta_t
snapshots = 1000
kT = 0.25
ratio = 0.5

# Set up the filename and log period
result_gsd_filename = 'MDrun.gsd'
log_filename = 'MDrun.log'
log_print_period = 1000

# Initialize configuration
# In this case, it considers 2-dimensional two-type LJ particles which are placed on square lattice.
rho = 1.00
LParticles = 32
NParticles = LParticles**2
NParticle_A = round(NParticles * ratio)
NParticle_B = round(NParticles * (1.0 - ratio))
Length = (NParticles / rho) ** 0.5

def place_square_lattice(position, typeid, diameter):
    num_A = 0
    num_B = 0
    for i in range(LParticles):
        for j in range(LParticles):
            position[i*LParticles+j,0] = Length*(i/LParticles-0.5)
            position[i*LParticles+j,1] = Length*(j/LParticles-0.5)
            position[i*LParticles+j,2] = 0
            diameter[i*LParticles+j] = 1.0
            if (uniform(0, 1) < ratio):
                if num_A < NParticle_A:
                    typeid[i*LParticles+j] = 0
                    num_A += 1
                else:
                    typeid[i*LParticles+j] = 1
                    num_B += 1
            else:
                if num_B < NParticle_B:
                    typeid[i*LParticles+j] = 1
                    num_B += 1
                else:
                    typeid[i*LParticles+j] = 0
                    num_A += 1

frame = gsd.hoomd.Frame()
frame.particles.N = NParticles
position = np.zeros((NParticles, 3))
diameter = np.zeros((NParticles,))
typeid = np.zeros((NParticles,))
place_square_lattice(position, typeid, diameter)
frame.particles.position = position
frame.particles.typeid = typeid
frame.particles.types = ['A', 'B']
frame.configuration.box = [Length, Length, 0, 0, 0, 0]
simulation.create_state_from_snapshot(frame)

# Define integrator and neighboring list
# WARNING: You need to specify explicit default_r_cut value,
#          otherwise, it will show undefined behavior!!!
integrator = md.Integrator(dt = deltat)
nl = md.nlist.Cell(buffer = 0.4, default_r_cut = 2.5)
# Add the nonreciprocal force in the integrator
# Usage: nrec.force.NonReciprocalForce(simulation, nlist, chi_par, chi_per, sigma, eps, r_cut)
nr = nrec.force.NonReciprocalForce(simulation, nl, chi_par, chi_per, 1.0, 2.0, 2.5)
integrator.forces.append(nr)

# Compute thermodynamic properties
thermodynamic_properties = md.compute.ThermodynamicQuantities(filter = hoomd.filter.All())
simulation.operations.computes.append(thermodynamic_properties)

# Set the Logger and gsd writer
samplingtime = int(totalsteps / snapshots)
if samplingtime == 0:
    samplingtime = 1
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(simulation, quantities=['timestep'])
logger.add(thermodynamic_properties, quantities=['kinetic_temperature', 'kinetic_energy'])
trigger = hoomd.trigger.Periodic(period = samplingtime, phase = 0)
gsd_writer = hoomd.write.GSD(trigger = trigger, filename = temp_gsd_filename, filter = hoomd.filter.All(), mode='wb', 
                            dynamic=['particles/position'])
simulation.operations.writers.append(gsd_writer)

# Show the status of simulation
class Status:
    def __init__(self, simulation):
        self.simulation = simulation

    @property
    def seconds_remaining(self):
        try:
            return (
                self.simulation.final_timestep - self.simulation.timestep
            ) / self.simulation.tps
        except ZeroDivisionError:
            return 0

    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining)).split('.',2)[0]

status_logger = hoomd.logging.Logger(categories=['scalar', 'string'])
status_logger.add(simulation, quantities=['timestep', 'tps'])
status = Status(simulation)
status_logger[('Status', 'etr')] = (status, 'etr', 'string')
table_status = hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=log_print_period), logger=status_logger)
simulation.operations.writers.append(table_status)

# Set up NVT thermostat
nvt = md.methods.ConstantVolume(filter = hoomd.filter.All(), thermostat=md.methods.thermostats.MTTK(kT = kT, tau = 50*deltat))
integrator.methods.append(nvt)
simulation.operations.integrator = integrator

# Run simulation
with open(log_filename, 'w') as f:
    table = hoomd.write.Table(trigger = trigger, logger = logger, output = f)
    simulation.operations.writers.append(table)
    simulation.run(totalsteps, write_at_start=True)
    gsd_writer.flush()
```

Note that the nonreciprocal force in this plugin is given by:

```math
\mathbf{F}_{\alpha \beta} = -[(1 + \chi_\parallel)\mathbf{\delta} + \chi_\perp \mathbf{\epsilon}] \nabla_{x_\alpha} u_{\alpha \beta}
```
```math
\mathbf{F}_{\beta \alpha} = -[(1 - \chi_\parallel)\mathbf{\delta} - \chi_\perp \mathbf{\epsilon}] \nabla_{x_\beta} u_{\alpha \beta},
```
where $\mathbf{\delta}$ is Kronecker delta tensor, and $\mathbf{\epsilon}$ is Levi-Civita tensor.

## **Developer Notes**

To make your own nonreciprocal force, you may need to modify `nonreciprocal_plugin.cc` and `nonreciprocal_pluginGPU.cu` file located in `src` directory.
To introduce new input to the function, you also need to modify header files, `nonreicprocal_plugin.h`, `nonreciprocal_pluginGPU.h`, and `nonreciprocal_pluginGPU.cuh`.
