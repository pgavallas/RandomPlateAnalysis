# Stochastic analysis of random composite plate

A cantilever plate is considered, consisting of a random short-fiber composite.

<img width="500" height="250" alt="cantilever" src="https://github.com/user-attachments/assets/d33c7507-ffe8-4cfa-b83a-b9a2a3947faa" />

- Its properties are simulated using an ABD matrix consisting of 8 components:
  

```math
\begin{bmatrix}
N_x \\
N_y \\
N_{xy} \\
M_x \\
M_y \\
M_{xy}
\end{bmatrix}
=
\begin{bmatrix}
A_{11} & A_{12} & 0 & 0 & 0 & 0 \\
A_{12} & A_{22} & 0 & 0 & 0 & 0 \\
0 & 0 & A_{33} & 0 & 0 & 0 \\
0 & 0 & 0 & D_{11} & D_{12} & 0 \\
0 & 0 & 0 & D_{12} & D_{22} & 0 \\
0 & 0 & 0 & 0 & 0 & D_{33}
\end{bmatrix}
\begin{bmatrix}
\varepsilon_x^0 \\
\varepsilon_y^0 \\
\gamma_{xy}^0 \\
\kappa_x \\
\kappa_y \\
\kappa_{xy}
\end{bmatrix}
```

- The random microstructure leads to random components.

- They are described using random fields, assumed Gaussian in this case.

- The random field properties, which include the means, standard deviations, correlation lengths and cross-correlation are extracted from a real composite microstructure through homogenization and a moving window technique (See https://doi.org/10.1016/j.compstruct.2025.119432).

- This script runs Monte Carlo Simulation (MCS) of the cantilever plate under a pressure load. For each iteration, a different random field realization is generated using the above properties and assigned to the elements, and finally the displacement of all nodes is computed. The statistics of the maximum displacement can then be investigated for a sufficient number of MCS.

<h3 align="center">Generated random field realization at one Monte Carlo iteration</h3>

<p align="center">
  <img src="https://github.com/user-attachments/assets/31e68f97-16d5-4edc-9217-2d4aa9e430ea" width="700"/>
</p>

>[!NOTE]
>random field generation is performed using https://github.com/pgavallas/MultivariateRF
