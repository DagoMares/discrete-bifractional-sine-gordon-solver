# Numerical Scheme for the 2D Dissipative Bifractional Sine-Gordon Equation
<p align="center">
  <img src="https://github.com/DagoMares/discrete-bifractional-sine-gordon-solver/blob/main/Soliton.gif" alt="Simulation of Soliton Wave evolution" />
</p>

## 📌 Project Overview

This repository hosts a high-performance **MATLAB** implementation of a second-order finite difference scheme designed to solve the **2D Dissipative Sine-Gordon Equation** featuring Riesz space-fractional operators of orders $\alpha, \beta \in (1, 2]$.

The primary objective of this project was to develop a robust numerical model capable of preserving the inherent physical structure of the system—specifically, the conservation of energy in the undamped case and the correct dissipation rate when damping is present. This structure-preserving property is critical for long-term simulation stability and physical accuracy.

### 📄 Associated Publication
The theoretical analysis and numerical results derived from this code have been published in the journal **Fractal and Fractional**.

---

## ⚛️ Mathematical Model

The project investigates the following nonlinear Partial Differential Equation (PDE):

$$
\frac{\partial^{2}u}{\partial t^{2}} + \gamma\frac{\partial u}{\partial t} - \lambda\Delta^{\alpha,\beta}u = -\phi(x,y)\sin(u) + F(x,y,t) - G^{\prime}(u)
$$

Where:
-   **$u(x,y,t)$**: The scalar field representing the wave solution.
-   **$\gamma > 0$**: Damping coefficient (dissipation).
-   **$\Delta^{\alpha,\beta}$**: The **Riesz Fractional Laplacian Operator** defined as:
    $$
    \Delta^{\alpha,\beta}u = \frac{\partial^{\alpha}u}{\partial|x|^{\alpha}} + \frac{\partial^{\beta}u}{\partial|y|^{\beta}}
    $$
    This fractional operator models nonlocal interactions in the medium, providing a generalization of standard diffusion.

---

## 🛠️ Numerical Methodology & Resources

This project leverages advanced numerical methods implemented in **MATLAB (R2024)**. The core algorithm preserves the discrete energy of the system, ensuring stability.

### Key Algorithms
*   **Temporal Discretization**: A two-step **Crank-Nicolson** method is used for time integration, ensuring unconditional stability and second-order accuracy.
*   **Spatial Approximation**: Riesz fractional derivatives are approximated using **fractional centered differences**, maintaining second-order spatial precision ($\mathcal{O}(h^2)$).
*   **Nonlinear Solver**: The arising nonlinear algebraic system at each time step is solved using a **Fixed-Point Iteration** scheme (Picard iteration) for efficiency and robustness.

### Theoretical Guarantees
1.  **Structure Preservation**: The scheme rigorously conserves discrete energy when $\gamma=0$ and correctly dissipates it when $\gamma > 0$, mirroring the continuous theorem (Theorem 2).
2.  **Consistency**: The method is formally proved to have a truncation error of $\mathcal{O}(h^2 + \tau^2)$.
3.  **Stability & Convergence**: The scheme is unconditionally stable and convergent in the discrete $L^2$-norm.

### Project Resources
The implementation relies on standard MATLAB libraries and custom optimization techniques:
-   **`simulation.m`**: The main driver script. It handles grid generation, matrix assembly, the time-stepping loop, and on-the-fly visualization.
-   **Matrix Operations**: Efficient sparse matrix handling for the discrete Laplacian operators.
-   **Special Functions**: Utilization of the `gamma` function to accurately compute fractional weights for the difference operators.
-   **Visualization**: Built-in `VideoWriter` for generating high-quality animations (MPEG-4) and `surf` for 3D surface rendering of the wave propagation.

---

## 🧪 Simulation Results

Extensive simulations validate the theoretical properties of the scheme:
*   **Accuracy Verification**: Comparison against analytical solutions yields $L^2$ errors consistently below $1.6 \times 10^{-2}$.
*   **Energy Dynamics**: Numerical experiments confirm the exact preservation of the discrete Hamiltonian in conservative regimes.
*   **Fractional Dynamics**: Varying $\alpha$ and $\beta$ reveals that **wave  amplitude significantly increases** as the fractional orders deviate from the classical integer case ($\alpha=\beta=2$), demonstrating the impact of nonlocality.

---

## 💻 Usage

To replicate the results or run your own simulations:

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/DagoMares/discrete-bifractional-sine-gordon-solver.git
    ```
2.  **Open MATLAB**: Navigate to the cloned directory.
3.  **Run the Simulation**:
    Execute the `simulation.m` script.
    ```matlab
    simulation
    ```
    *   Parameters such as `alpha`, `beta`, grid size, and time steps can be modified directly in the script header.

---

## 👥 Authors & Contact

**Dagoberto Mares-Rincón**, Siegfried Macías, Jorge E. Macías-Díaz, José A. Guerrero-Díaz-de-León, and Tassos Bountis.

For inquiries or collaborations, please feel free to reach out:

[![Gmail Badge](https://img.shields.io/badge/-dagobertomares0@gmail.com-c14438?style=flat&logo=Gmail&logoColor=white&link=mailto:dagobertomares0@gmail.com)](mailto:dagobertomares0@gmail.com)
[![Linkedin Badge](https://img.shields.io/badge/-dagoberto--mares-0072b1?style=flat&logo=Linkedin&logoColor=white&link=https://www.linkedin.com/in/dagoberto-mares/)](https://www.linkedin.com/in/dagoberto-mares/)
