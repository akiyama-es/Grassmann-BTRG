# Grassmann-BTRG

A sample implementation of the bond-weighted tensor renormalization group (BTRG) for lattice fermions.

## Model

The two-dimensional single-flavor Gross-Neveu-Wilson model at finite density is provided as a specific example. The model is defined by the following action,
$$S=-\frac{1}{2}\sum_{n\in\Lambda_{2}}\sum_{\nu=1,2}\left[{\mathrm{e}}^{\mu\delta_{\nu,2}} \bar{\psi}(n)(r1-\gamma_{\nu})\psi(n+\hat{\nu})+{\mathrm{e}}^{-\mu\delta_{\nu,2}} \bar{\psi}(n+\hat{\nu})(r1+\gamma_{\nu})\psi(n)\right]+(M+2r)\sum_{n}\bar{\psi}(n)\psi(n)-\frac{g^2}{2}\sum_{n}\left[(\bar{\psi}(n)\psi(n))^{2}+(\bar{\psi}(n){\mathrm{i}}\gamma_{5}\psi(n))^{2}\right]$$
with $\psi(n)$ and $\bar{\psi}(n)$ are the two-component Grassmann fields. $r$ is the Wilson parameter, which is set to be $r=1$ in the code. The model is defined on a square lattice, where we have the periodic boundary condition for the direction of $\nu=1$ and the anti-periodic one for $\nu=2$. The two-dimensional $\gamma$-matrices are given by the Pauli matrices via $\gamma_{1}=\sigma_{x}$, $\gamma_{2}=\sigma_{y}$ and $\gamma_{5}=\sigma_{z}$. 
**$M$, $\mu$, $g^{2}$ denote the mass, chemical potential, and four-point coupling, which are the parameters you can set freely.**

## Algorithmic parameters

The BTRG approximately calculates the path integral,
$$Z=\int[{\mathrm{d}}\psi{\mathrm{d}}\bar{\psi}]{\mathrm e}^{-S}$$,
which is represented as the Grassmann tensor network,
$$Z={\mathrm{gTr}}\left[\prod_{n}\mathcal{T}_{n}\right]$$.
**You can freely choose the bond dimension $D$, iteration number $N$, which is interpreted as the system volume $V$ via $V=2^{N}$, and the hyperparameter $k$.** With $k=0$, the BTRG is reduced to be the Levin-Nave TRG. The optimal choice of $k$ on a square lattice is $k=-0.5$. Please see the references shown below.

## Note

You can specify the parameters $M$, $\mu$, $g^{2}$, $D$, $N$, and $k$ from the command line by executing the file named "GBTRG". If you set $g^{2}=0$, then the code can calculate the relative error of ${\mathrm{ln}}Z$, based on the diagonalization of the Dirac matrix in the momentum space. However, the diagonalization may require much longer time compared with the BTRG calculation, so you can omit it from the command line if you want.

## Reference

- [D. Adachi, T. Okubo, and S. Todo, Bond-weighted Tensor Renormalization Group, PRB105(2022)L060402](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.L060402).
- [S. Akiyama, Bond-weighting method for the Grassmann tensor renormalization group, JHEP11(2022)030](https://link.springer.com/article/10.1007/JHEP11(2022)030).
- [S. Akiyama and D. Kadoh, More about the Grassmann tensor renormalization group, JHEP10(2021)188](https://link.springer.com/article/10.1007/JHEP10(2021)188).
