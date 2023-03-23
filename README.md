# FidelityEchoes GUI

Matlab App interactive GUI to compute Fidelity and Loschmidt Echoes for quadratic bosonic Hamiltonians.

Given: 
  - a Quadratic Bosonic Hamiltonian

  $$ H=(\alpha + \beta)(a^\dagger a + \frac{1}{2}) + \frac{(\beta - \alpha)}{2}(a^{\dagger 2} + a^2 ) $$

defined by the bosonic creation (anhilation) operators $a^\dagger$ ($a$) and parameters $\alpha$ and $\beta$.

  - a perturbated Hamiltonian $H_2 = H + \epsilon P$, with $P$ also in the quadratic form
  
  $$ P=(a + b)(a^\dagger a + \frac{1}{2}) + \frac{(b - a)}{2}(a^{\dagger 2} + a^2 ) $$
  
  and a,b defined as $a=\cos(\theta)$ and $b=\sin(\theta)$.
