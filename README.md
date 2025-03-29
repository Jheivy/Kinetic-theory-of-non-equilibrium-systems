# Kinetic Theory of Non-Equilibrium Systems  

This project implements a Fortran program to study non-equilibrium systems based on kinetic theory. The simulation employs Monte Carlo methods to model particle collisions and analyze the evolution of the system.  

## Key Features  
- **H-Theorem Calculation**: Tracks the evolution of the system's entropy using Boltzmann's H-Theorem.  
- **Velocity Distribution Analysis**: Computes histograms of particle speed magnitudes at each Monte Carlo step to observe system dynamics.  
- **Efficient Collision Modeling**: Simulates particle interactions using probabilistic methods to explore non-equilibrium behavior.  

This project contributes to understanding the microscopic foundations of thermodynamic irreversibility in non-equilibrium systems.  

Feel free to explore, fork, or contribute to the repository. Suggestions for improvements and extensions are welcome!  

```fortran
program kinetic_theory
  implicit none
  print *, "Welcome to the Kinetic Theory Simulation"
end program kinetic_theory
