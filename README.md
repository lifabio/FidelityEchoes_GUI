# FidelityEchoes GUI

Matlab App interactive GUI to compute Fidelity and Loschmidt Echoes for quadratic bosonic Hamiltonians.

Given: 
  - a Quadratic Bosonic Hamiltonian

  $$ H=(\alpha + \beta)(a^\dagger a + \frac{1}{2}) + \frac{(\beta - \alpha)}{2}(a^{\dagger 2} + a^2 ) $$

defined by the bosonic creation (anhilation) operators $a^\dagger$ ($a$) and parameters $\alpha$ and $\beta$.

  - a perturbated Hamiltonian $H_2 = H + \epsilon P$, with $P$ also in the quadratic form
  
  $$ P=(c + s)(a^\dagger a + \frac{1}{2}) + \frac{(s - c)}{2}(a^{\dagger 2} + a^2 ) $$
  
  and perturbation parameters $c$, $s$ defined as $c=\cos(\theta)$ and $s=\sin(\theta)$.
  
  The algorithm takes as input parameters: $\alpha$, $\beta$, $\epsilon$ and $\theta$ and computes the quantum Fidelity $F$ and Loschmidt Echo $M$ as
  
  $$ F= \frac{||\langle\Psi_2(t)|\Psi_1(t)\rangle||^2}{\langle\Psi_1(t)|\Psi_1(t)\rangle\langle\Psi_2(t)|\Psi_2(t)\rangle} $$
  
  $$ M= $$
