# Grassmann-BTRG

A sample implementation of the bond-weighted tensor renormalization group (BTRG) for lattice fermions.

## Model

The two-dimensional single-flavor Gross-Neveu-Wilson model at finite density is provided as a specific example.
$$S=-\frac{1}{2}\sum_{n\in\Lambda_{2}}\sum_{\nu=1,2}\left[{\mathrm{e}}^{\mu\delta_{\nu,2}}\bar{\psi}(n)(r1-\gamma_{\nu})\psi(n+\hat{\nu})+{\mathrm{e}}^{-\mu\delta_{\nu,2}}\bar{\psi}(n+\hat{\nu})(r1+\gamma_{\nu})\psi(n)\right]+(M+2r)\sum_{n}\bar{\psi}(n)\psi(n)-\frac{g^2}{2}\sum_{n}\left[(\bar{\psi}(n)\psi(n))^{2}-(\bar{\psi}(n){\mathrm{i}}\gamma_{5}\psi(n))^{2}\right]$$
with $\psi(n)$ and $\bar{\psi}(n)$ are the two-component Grassmann fields. $r$ is the Wilson parameter, which is set to be $r=1$ in the code. 
**$M$, $\mu$, $g^{2}$ denote the mass, chemical potential, and four-point coupling, which are the parameters you can set freely.**

## Algorithmic parameters

BTRG approximately calculates a tensor contraction. 
**You can freely choose the bond dimension $D$, iteration number $n$, which is interpreted as the system volume $V$ via $V=2^{n}$, and the hyperparameter $k$.** With $k=0$, the BTRG is reduced to be the Levin-Nave TRG. The optimal choice of $k$ on a square lattice is $k=-0.5$. Please see the references shown below.

## Note

You can specify the parameters (both for the model and the BTRG) from the command line by executing the file named "GBTRG".

## Reference

- [D. Adachi, T. Okubo, and S. Todo, Bond-weighted Tensor Renormalization Group, PRB105(2022)L060402](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.L060402).
- [S. Akiyama, Bond-weighting method for the Grassmann tensor renormalization group, JHEP11(2022)030](https://link.springer.com/article/10.1007/JHEP11(2022)030).
