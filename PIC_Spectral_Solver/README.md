# PIC Spectral Solver

The spectral solver is indifferent to the linear solver in terms of classes and functions, but the pic.py file is set to run for a Spectral Poisson solver, rather than solving the Poisson using the finite difference method. Moreover, it does not use the 'celluloid' and 'camera' modules to produce animations, rather it outputs images which are to be assembled into an animation by the user.

Currently, the settings are set to simulate a million protons interacting and colliding ellastically with the walls of a square 2D box for 200 s, and to output the Energy plot. The construction of figures for particle positions and potential plots can be initialised by setting the 'animate_particle_motion' and 'animate_potential' flags to true at the beginning of the pic.py file. This will make a plot at each full time step (0.02 s, in this simulation).
