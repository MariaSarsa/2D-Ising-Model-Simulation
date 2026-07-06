# 2D-Ising-Model-Simulation
A C-based Monte Carlo simulation of the 2D Ising Model on a square lattice using the **Metropolis algorithm** for state transitions and the **Parisi-Rapuano method** method for random number generation.

This codebase tracks system energy and magnetization across inverse temperature steps ($\beta$) to study phase transitions and computes statistical error bounds for higher-order thermodynamic fluctuations.

---

## Context
The Ising Model is a Physics model that explains ferromagnetism. In nature, a magnetic material is made up of tiny microscopic magnets (spins) that are arranged on a grid. In a 2D grid, each spin has two states: (+1) or (-1).

As the temperature of the system is changed, the behaviour of the material changes. At a certain temperature there is a phase transition, called the critial point, the boundary state between order and chaos.

In particular, at high temperatures (Low $\beta$), the spins are at random states, so the overall magnetism stays close to zero.
Whereas at low temperatures (High $\beta$), the spins are aligned in the same direction, making the material strongly magnetic.

---

## Features

- **Dynamic Thermalization Tracking**: Before taking real data, the code lets the virtual material settle into its new temperature.

- **Hysteresis Loop Evaluation**: Real magnets have a "memory" effect. If you heat a magnet up and then cool it down, it won't change its magnetism at the exact same temperature thresholds in both directions. This is called hysteresis.

- **Finding the Critical Transition Point ($\beta_c$)**: By observing the peaks in the Specific Heat and Magnetic Susceptibility data across different grid sizes, the point where the material flips from random chaos to an ordered magnet can be tracked.



---

## Directory Structure

The structure of the project follows as:

```text
├── .gitignore
├── README.md
├── input.txt
├── include/
│   └── Cabecera.h
└── src/
    ├── avg_deviaation.c
    ├── initialConfiguration.c
    ├── energy.c
    ├── Initialize.c
    ├── magnetization.c
    ├── Metropolis.c
    ├── probability.c
    ├── Build_Blocks.c
    └── main.c
