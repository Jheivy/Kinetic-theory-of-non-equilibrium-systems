# Kinetic Theory of Non-Equilibrium Systems  

## Overview  
This project implements a Fortran program to study non-equilibrium systems based on kinetic theory. Using Monte Carlo methods, the simulation models particle collisions and provides insights into the evolution of thermodynamic properties in non-equilibrium states.  

## Key Features  
- **H-Theorem Calculation**: Tracks the evolution of the system's entropy using Boltzmann's H-Theorem.
- ![Simulation Example](Graphics/part2/teorema_H.pdf)  
- **Velocity Distribution Analysis**: Computes histograms of particle speed magnitudes at each Monte Carlo step to observe system dynamics.  
- **Efficient Collision Modeling**: Simulates particle interactions using probabilistic methods to explore non-equilibrium behavior.  
- **Customizable Parameters**: Allows adjustments of initial velocity distributions, particle count, and simulation steps for tailored studies.  

This project contributes to understanding the microscopic foundations of thermodynamic irreversibility in non-equilibrium systems.  

Feel free to explore, fork, or contribute to the repository. Suggestions for improvements and extensions are welcome!  

```fortran
program kinetic_theory
  implicit none
  print *, "Welcome to the Kinetic Theory Simulation"
end program kinetic_theory
