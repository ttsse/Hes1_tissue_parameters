# Parameter Sweep for Hes1 Tissue Model
We try to find parameters which result in stable oscillations for the given PDE system

$$
\begin{align}
    \frac{\partial d}{\partial t} &= D_d \frac{\partial^2 d}{\partial x^2} + \alpha_d n - \mu_d d,\\
    \frac{\partial m}{\partial t} &= \frac{\alpha_m}{1+p^h} - \mu_m m,\\
    \frac{\partial p}{\partial t} &= \alpha_p m - \mu_p p,\\
    \frac{\partial n}{\partial t} &= \frac{\alpha_n (1+d)}{1+p^{\gamma}} - \mu_n n, \\
\end{align}
$$

describing how the protein Hes1 is involved in signalling between cells during embryonal development using Matlab's PDE Toolbox. In this case the involved molecules are

$$
\begin{align}
    &\text{d(t,x): Dll1 protein}, \\
    &\text{m(t,x): Hes1 mRNA}, \\
    &\text{p(t,x): Hes1 protein}, \\
    &\text{n(t,x): Ngn2 protein}.
\end{align}
$$

Since data concerning the required parameter values are difficult to determine, we perform a parameter sweep over several parameters in the system to determine their influence on the system's behaviour and, hopefully, find more fitting parameter values than the initially chosen ones.


## Computational Environment
The code is written in the MATLAB computational environment. The used release is MATLAB release R2021b Update 6 (version number: 9.11). Apart from that, the code (including the tests) uses the following toolboxes:

| Toolbox                                   | Version |
|-------------------------------------------|---------|
| Signal Processing Toolbox                 | '8.7'   |
| Parallel Computing Toolbox                | '7.5'   |
| MATLAB Parallel Server                    | '7.5'   |
| Polyspace Bug Finder                      | '3.5'   |
| Antenna Toolbox                           | '5.1'   |
| Curve Fitting Toolbox                     | '3.6'   |
| Deep Learning HDL Toolbox                 | '1.2'   |
| Fixed-Point Designer                      | '7.3'   |
| System Identification Toolbox             | '9.15'  |
| MATLAB Coder                              | '5.3'   |
| Optimization Toolbox                      | '9.2'   |
| Simulink                                  | '10.4'  |
| Statistics and Machine Learning Toolbox   | '12.2'  |
| Computer Vision Toolbox                   | '10.1'  |
