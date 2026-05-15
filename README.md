# Analytical and Numerical Methods for Differential Equations (SF1682)

This repository contains my solutions for the course **SF1682 Analytical and Numerical Methods for Differential Equations** at KTH Royal Institute of Technology, taken during the Autumn Semester 2024 (HT24).

The project involves modeling the kinematics of a differential-steered robot using MATLAB to simulate trajectories through various numerical schemes.

---

## 🛠️ Course Overview & Technologies
* **Language:** MATLAB
* **Concepts:** ODEs, Numerical Integration (Euler, RK4), Convergence Analysis, and Newton-Raphson.

---

## 📂 Project Structure: Project 1A 

### [Task 1-3] Kinematics & Numerical Schemes
* **Focus:** Differential Steering & Integration
* Derived analytical models for robot motion and implemented **Forward Euler** and **Runge-Kutta 4** solvers. Conducted a convergence study to verify the global truncation error and accuracy order ($p=1$ and $p=4$) for each method.

### [Task 4-5] Adaptive Solvers
* **Focus:** Error Control
* Simulated complex, non-circular trajectories using MATLAB’s adaptive `ode45` solver, comparing numerical approximations against analytical benchmarks under strict tolerance levels.

### [Task 6] Root-finding Optimization
* **Focus:** Newton-Raphson Method
* Designed a solver to determine specific acceleration parameters required to hit target coordinates. Verified the algorithm's performance through quadratic convergence analysis.

## 📂 Project Structure: Project 1B (Suspension Systems)

### [Task 1-2] Quarter-Car Modeling & Simulation
* **Focus:** Mass-Spring-Damper Systems
* Modeled a "quarter-car" suspension system as a system of second-order ODEs. Implemented `quartercar.m` and used MATLAB’s `ode45` to simulate the vehicle's vertical response when driving over a road bump, analyzing chassis and wheel displacement ($z_1, z_2$).

### [Task 3] Optimization & Newton-Raphson
* **Focus:** Multivariable Root-finding
* Optimized spring constants $k_1$ and $k_2$ to balance passenger comfort and road handling. Developed a multivariable **Newton-Raphson** solver using a custom Jacobian matrix (`Jacobian_transfer_functions.m`) to find the optimal parameters where the transfer functions meet specified safety and comfort criteria.

### [Task 4] Stability & Stiffness Analysis
* **Focus:** Eigenvalues & Explicit Methods
* Analyzed the system's stability by calculating the eigenvalues of the state-space matrix $H$. Determined the maximum stable time step ($\Delta t_{max}$) for the **Forward Euler** method and demonstrated how "stiff" differential equations cause numerical instability when the step size exceeds theoretical limits.

### [Task 5] Implicit Methods for Stiff Problems
* **Focus:** Trapezoidal Rule
* Implemented the **Implicit Trapezoidal Method** to handle stiff ODEs. Verified that implicit solvers maintain stability even with significantly larger time steps compared to explicit methods, and conducted a convergence study to confirm the method's second-order accuracy.

## 📂 Project Structure: Project 2A (Heat Equation & Finite Differences)

### [Task 1] Numerical Differentiation
* **Focus:** Taylor Expansion & Error Analysis
* Evaluated central, forward, and backward difference schemes for first-order derivatives. Performed a convergence study using log-log plots to identify the transition point where truncation error is superseded by floating-point round-off error.

### [Task 2] Boundary Value Problems (BVP)
* **Focus:** Steady-state Heat & Boundary Conditions
* Solved 1D BVPs using the **Finite Difference Method (FDM)**. Constructed system matrices to handle **Dirichlet** and **Neumann** boundary conditions, verifying second-order accuracy ($O(h^2)$) through grid refinement studies.

### [Task 3] Diffusion & Explicit Methods
* **Focus:** Heat Equation & Stability
* Modeled time-dependent heat distribution using the **FTCS** (Forward Time Central Space) scheme. Applied **Von Neumann stability analysis** to determine the critical time step $\Delta t_{max}$ and demonstrated how exceeding this limit causes numerical divergence.

### [Task 4] Implicit Methods & Crank-Nicolson
* **Focus:** Unconditional Stability
* Implemented the **Crank-Nicolson** method for solving non-homogeneous heat equations with time-varying source terms. Confirmed the method's second-order convergence in both space and time and its ability to maintain stability with larger time steps compared to explicit solvers.

## 📂 Project Structure: Project 2B (Spectral Methods & Fourier Transforms)

### [Task 1] Spectral Differentiation
* **Focus:** Fast Fourier Transform (FFT)
* Implemented numerical differentiation in the frequency domain using `fft` and `ifftshift`. Demonstrated the "spectral accuracy" of the method, showing that for smooth, periodic functions, the approximation error decreases exponentially with $N$ until reaching machine precision.

### [Task 2] Diffusion & Resolution Analysis
* **Focus:** High-resolution Reference Modeling
* Solved the 1D diffusion equation using a spectral approach to analyze how spatial resolution affects the numerical solution. Compared low-resolution results against a high-resolution reference ($N=2^{11}$) to determine the $L_2$-error convergence for non-trivial initial conditions.

### [Task 3] Time-stepping with Runge-Kutta
* **Focus:** Hybrid Spectral-RK4 Solver
* Developed a solver for the non-homogeneous heat equation by combining spectral spatial discretization with a fourth-order **Runge-Kutta (RK4)** method for time integration. Visualized the results using 3D waterfall plots to analyze the temporal evolution of the heat distribution.

### [Task 4] Stability & Grid Refinement
* **Focus:** Numerical Stability Analysis
* Investigated the stability limits of spectral methods in time-dependent problems. Performed grid refinement studies to observe how the choice of time step $\Delta t$ and the number of modes $N$ interact to maintain a stable and physically accurate solution.

