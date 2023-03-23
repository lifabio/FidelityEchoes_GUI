# FidelityEchoes GUI

Matlab App interactive GUI to compute Fidelity and Loschmidt Echoes for quadratic bosonic Hamiltonians.

Given: 
  - a Quadratic Bosonic Hamiltonian

  $$ H=(\alpha + \beta)(a^\dagger a + \frac{1}{2}) + \frac{(\beta - \alpha)}{2}(a^{\dagger 2} + a^2 ) $$

defined by the bosonic creation (anhilation) operators $a^\dagger$ ($a$) and parameters $\alpha$ and $\beta$.

  - a perturbated Hamiltonian $H_2 = H + \epsilon P$, with $P$ also in the quadratic form
  
  $$ P=(c + s)(a^\dagger a + \frac{1}{2}) + \frac{(s - c)}{2}(a^{\dagger 2} + a^2 ) $$
  
  and perturbation parameters $c$, $s$ defined as $c=\cos(\Theta)$ and $s=\sin(\Theta)$.
  
  The algorithm takes as input parameters: $\alpha$, $\beta$, $\epsilon$, $\Theta$ and initial state $|\Psi_0\rangle$ and computes the quantum Fidelity $F(t)$ and Loschmidt Echo $M(t)$ as a function of time as
  
  $$ F(t)= \frac{||\langle\Psi_2(t)|\Psi_1(t)\rangle||^2}{\langle\Psi_1(t)|\Psi_1(t)\rangle\langle\Psi_2(t)|\Psi_2(t)\rangle} $$
  
  $$ M(t)= \frac{||\langle\Psi_0|\Psi_f(t)\rangle||^2}{\langle\Psi_0|\Psi_0\rangle\langle\Psi_f(t)|\Psi_f(t)\rangle} $$

where $|\Psi_1(t)\rangle = e^{-i H t}|\Psi_0\rangle$ and $|\Psi_2(t)\rangle = e^{-i H_2 t}|\Psi_0\rangle$ are the forward propagation of the initial state $|\Psi_0\rangle$ according to $H$ and $H_2$ respectivelly. $|\Psi_f(t)\rangle = e^{+i H_2 t}e^{-i H t}|\Psi_0\rangle$ is the forward and then backward propagation of $|\Psi_0\rangle$ according to $H$ and $H_2$.

The initial state $|\Psi_0\rangle$ can be chosen from:
  - Equal superposition of Fock states $|\Psi_0\rangle = \frac{|1\rangle + |2\rangle + ... + |N\rangle}{\sqrt{N}}$
  - Coherent State $|\Psi_0\rangle = |z\rangle$ with $z=|\alpha_0|e^{i\sigma}$
  - Vacuum State $|\Psi_0\rangle = |0\rangle$
  - Single Fock state $|\Psi_0\rangle = |n_0\rangle$

Figure  shows the visual appearence of the GUI.
  


