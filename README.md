# CP(2) Worm Algorithm in C++

I wrote this code as a postdoc to perform the research carried out in this paper: https://arxiv.org/abs/1803.04767.

## Overview of the Worm Algorithm
We simulate a spin system which has two spatial dimensions and a time dimension. 

Hence a 3-dimensional lattice is used. On each lattice site there sits a CP(2) degree of freedom. The details of CP(2) - a mathematical symmetry group - are not important for those only interested in the software engineering aspects of this repository. It suffices to imagine a lattice of sites, with a (spin) value at each site that can take one of three values, let's say 0, 1 or 2.

We would like to update the current lattice configuration, a configuration with of 0s, 1s and 2s at each lattice site, so that we can generate thousands of different configurations. Once we have an ensemble of configurations we can then perform measurements one each configuration and by averaging the results we can learn about the physics of the CP(2) spin system.

To update the current configuration we use the worm algorithm:

i) First we randomly pick an initial time direction and starting site. At this point imagine the head and tail of the worm are both at the starting point.
ii) Then according to rules determined by the physics of the system, the value at the lattice site is updated and the head of the worm moves to another lattice site, this move can be:
a) spatially - if so the time direction of the worm also reverses
b) temporally
c) diagonally - spatially and temporally simultaneously.
iii) Step ii) is repeated until the head of the worm finds its tail, at which point the update is complete, and a new valid configuration has been produced.
