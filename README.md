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
```

---

## input.txt guide

The program expects 6 lines of numbers in the specific order shown below:

| Line | Value | Variable Name | Description |
| :---: | :---: | :--- | :--- |
| **1** | `0.0` | `b_initial` | **Initial Inverse Temperature ($\beta$):** The starting point of the simulation sweep. `0.0` represents infinite temperature ($T \to \infty$) where spins flip completely randomly. |
| **2** | `2.0` | `b_final` | **Final Inverse Temperature ($\beta$):** The upper boundary target where the temperature sweep stops. |
| **3** | `0.05`| `b_delta` | **Temperature Step Size:** The amount added to $\beta$ after completing each simulation run to move to the next temperature coordinate. |
| **4** | `1000000` | `Nmed` | **Number of Measurements:** The quantity of distinct data samples captured per temperature point once the system reaches equilibrium. |
| **5** | `5` | `L` | **Lattice Size ($L$):** The linear size of the square grid. A value of `5` builds a $5 \times 5$ system containing 25 total spin sites. |
| **6** | `10000` | `Nter` | **Maximum Warm-Up Steps:** The limit for thermalization passes. The code stops early if stability is detected, but will stop at this number if it takes too long. |

---

## Simulation Visualizations

### Termalization time

<table>
  <thead>
    <tr>
      <th width="10%">Condition</th>
      <th width="30%">Energy</th>
      <th width="30%">Magnetization</th>
      <th width="30%">Spin Configuration</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center">$\beta < \beta_c$<br><small>(High T)</small></td>
      <td><img  alt="Energy at low beta" src="https://github.com/user-attachments/assets/0e609740-24b9-43dc-939e-b44d7fad785a" width="100%" /></td>
      <td><img  alt="Magnetization at low beta" src="https://github.com/user-attachments/assets/6af05df5-b1e4-4c41-8a7e-828c6d3f7338" width="100%"/></td>
      <td><img  alt="Spins at low beta" src="https://github.com/user-attachments/assets/529a6755-7402-4aa7-8f99-924a178569a7" width="80%" /></td>
    </tr>
    <tr>
      <td align="center">$\beta = \beta_c$<br><small>(Critical Point)</small></td>
        <td><img  alt="Energy at critical beta" src="https://github.com/user-attachments/assets/5f83ad3d-2be9-4732-b305-bc3d16840d80" width="100%"/></td>
        <td><img  alt="Magnetization at critical beta" src="https://github.com/user-attachments/assets/6167aef8-43cf-4153-a025-080d59524af4"  width="100%"/></td>
        <td><img  alt="Spins at critical beta" src="https://github.com/user-attachments/assets/6148aba0-e30e-4704-b5b4-784021583a21" width="80%"/></td>
    </tr>
    <tr>
      <td align="center">$\beta > \beta_c$<br><small>(Low T)</small></td>
        <td><img  alt="Energy at high beta" src="https://github.com/user-attachments/assets/5a14de58-ad3c-4af5-a905-2412ef332a26" width="100%"/></td>
        <td><img  alt="Magnetization at high beta" src="https://github.com/user-attachments/assets/b9b3e550-50e2-4854-bd23-d40154f0f4d7" width="100%"/></td>
        <td><img  alt="Spins at high beta" src="https://github.com/user-attachments/assets/4086ca02-83dc-4b9d-b25e-b4798a12e981" width="80%"/></td>
    </tr>
  </tbody>
</table>


### Hysteresis Cycles

<table>
  <tbody>
    <tr>
      <td>
        <img alt="Magnetization Hysteresis L=8 & 64" src="https://github.com/user-attachments/assets/25bf94c8-2d6e-4077-9c17-dde52a5eb4f4" width="100%"/>
      </td>
      <td>
        <img alt="Magnetization Hysteresis L=16 & 32" src="https://github.com/user-attachments/assets/bc32795d-a779-4205-9a8e-05fbc534d1a9" width="100%"/>
      </td>
    </tr>
  </tbody>
</table>

### Termalized Observables

<table>
  <thead>
    <tr>
      <th width="33.3%">System Energy</th>
      <th width="33.3%">Magnetization</th>
      <th width="33.3%">Specific Heat ($C_v$)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>
        <img alt="Thermalized Energy" src="https://github.com/user-attachments/assets/2407d4ad-c3aa-455d-9336-92c1d95cfbf4" width="100%"/>
      </td>
      <td>
        <img alt="Thermalized Magnetization" src="https://github.com/user-attachments/assets/741c6ec7-6541-404d-b98a-783408546d93" width="100%"/>
      </td>
      <td>
        <img alt="Thermalized Specific Heat" src="https://github.com/user-attachments/assets/2df6c85c-be29-4d66-8f26-9e2be48b447a" width="100%"/>
      </td>
    </tr>
  </tbody>
</table>

### Critical Beta Estimation

<table>
  <thead>
    <tr>
      <th align="center">Critical Point Estimation ($\beta_c$)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center">
        <img alt="Critical Beta Estimation Plot" src="https://github.com/user-attachments/assets/0e0784c7-9723-4ff0-8cbe-348a7fab0466" width="80%" />
        <br>
      </td>
    </tr>
    <tr>
      <td align="center" bgcolor="#f6f8fa">
        <img alt="Estimated Beta Value" src="https://github.com/user-attachments/assets/ca7ccf18-2a7e-4566-8bea-19557371e11f" width="50%" />
        <br>
        <sub><b>Figure:</b> Final calculated numerical estimation for $\beta_c$.</sub>
      </td>
    </tr>
  </tbody>
</table>
