

# 1. Two-spin Hilbert space

Take two electron spins, one in each dot. The Hilbert space is

$$
\mathcal H = \mathbb C^2 \otimes \mathbb C^2
$$

with basis

$$
{|{\uparrow\uparrow}\rangle,\ |{\uparrow\downarrow}\rangle,\ |{\downarrow\uparrow}\rangle,\ |{\downarrow\downarrow}\rangle}.
$$

Let $\sigma_\alpha^{(1)} = \sigma_\alpha \otimes I$ and $\sigma_\alpha^{(2)} = I \otimes \sigma_\alpha$, for $\alpha=x,y,z$.

Also define spin operators

$$
S_\alpha^{(i)} = \frac{1}{2}\sigma_\alpha^{(i)}.
$$

---

# 2. Concrete Hamiltonian

A useful model is

$$
H(t)=H_Z + H_{\mathrm{OH}} + H_{\mathrm{ex}}(t) + H_{\mathrm{SOC}}(t).
$$

I’ll define each term.

## 2.1 Static Zeeman term

Assume an external magnetic field mainly along (z):

$$
H_Z = \frac{\omega_1}{2}\sigma_z^{(1)} + \frac{\omega_2}{2}\sigma_z^{(2)}.
$$

Equivalently,

$$
H_Z = \frac{\omega_{\mathrm{avg}}}{2}\left(\sigma_z^{(1)}+\sigma_z^{(2)}\right)+\frac{\Delta\omega}{2}\left(\sigma_z^{(1)}-\sigma_z^{(2)}\right),
$$

where

$$
\omega_{\mathrm{avg}}=\frac{\omega_1+\omega_2}{2},\qquad
\Delta\omega=\frac{\omega_1-\omega_2}{2}.
$$

The average term is often just a rotating-frame offset. The difference term is what matters most for controllability.


## 2.2 Overhauser field

Model the nuclear environment as effective local magnetic fields (\mathbf b_1,\mathbf b_2):

$$
H_{\mathrm{OH}} =
\frac{1}{2}\mathbf b_1\cdot \boldsymbol{\sigma}^{(1)}
+
\frac{1}{2}\mathbf b_2\cdot \boldsymbol{\sigma}^{(2)}.
$$

Expanding:

$$
H_{\mathrm{OH}} =
\frac{1}{2}\sum_{\alpha=x,y,z}
\left(
b_{1\alpha}\sigma_\alpha^{(1)} + b_{2\alpha}\sigma_\alpha^{(2)}
\right).
$$

Often the dominant piece is the (z)-gradient:

$$
H_{\mathrm{OH}}^{(z)} =
\frac{b_{1z}}{2}\sigma_z^{(1)}+\frac{b_{2z}}{2}\sigma_z^{(2)}
\frac{b_{\mathrm{avg}}}{2}\left(\sigma_z^{(1)}+\sigma_z^{(2)}\right)
+
\frac{\Delta b_z}{2}\left(\sigma_z^{(1)}-\sigma_z^{(2)}\right).
$$

This $\Delta b_z$ is the usual singlet-triplet dephasing source.


## 2.3 Exchange control

Barrier or detuning control gives an exchange interaction

$$
H_{\mathrm{ex}}(t)=J(t),\mathbf S_1\cdot\mathbf S_2
$$

or

$$
H_{\mathrm{ex}}(t)=\frac{J(t)}{4}
\left(
\sigma_x^{(1)}\sigma_x^{(2)}
+\sigma_y^{(1)}\sigma_y^{(2)}
+\sigma_z^{(1)}\sigma_z^{(2)}
\right).
$$

This is the standard isotropic Heisenberg coupling.

Usually $J(t)$ is voltage-controlled:

$$
J(t)=J(V(t)).
$$

---

## 2.4 Rashba-induced spin-orbit term

A practical effective form is an antisymmetric exchange, i.e. Dzyaloshinskii–Moriya-like coupling:

$$
H_{\mathrm{SOC}}(t)= \mathbf D(t)\cdot(\mathbf S_1\times \mathbf S_2).
$$

In Pauli form:

$$
H_{\mathrm{SOC}}(t)=
\frac{1}{4}\mathbf D(t)\cdot
\begin{pmatrix}
\sigma_y^{(1)}\sigma_z^{(2)}-\sigma_z^{(1)}\sigma_y^{(2)}\\
\sigma_z^{(1)}\sigma_x^{(2)}-\sigma_x^{(1)}\sigma_z^{(2)}\\
\sigma_x^{(1)}\sigma_y^{(2)}-\sigma_y^{(1)}\sigma_x^{(2)}
\end{pmatrix}.
$$

A common simplification is that $\mathbf D(t)$ points along one axis, say (z):

$$
H_{\mathrm{SOC}}(t)=
\frac{D_z(t)}{4}
\left(
\sigma_x^{(1)}\sigma_y^{(2)}-\sigma_y^{(1)}\sigma_x^{(2)}
\right).
$$

This is often enough to study controllability.

You can also add symmetric anisotropic exchange:

$$
H_{\mathrm{aniso}}(t)=
\sum_{\alpha,\beta=x,y,z}\Gamma_{\alpha\beta}(t)S_\alpha^{(1)}S_\beta^{(2)},
$$

but for most first studies the DM term is the key SOC contribution.

# 3. Full concrete Hamiltonian

A clean practical model is:

$$
\boxed{
H(t)=
\frac{\omega_1+b_{1z}}{2}\sigma_z^{(1)}
+
\frac{\omega_2+b_{2z}}{2}\sigma_z^{(2)}
+
\frac{J(t)}{4}
\left(
\sigma_x^{(1)}\sigma_x^{(2)}+
\sigma_y^{(1)}\sigma_y^{(2)}+
\sigma_z^{(1)}\sigma_z^{(2)}
\right)
+
\frac{D_z(t)}{4}
\left(
\sigma_x^{(1)}\sigma_y^{(2)}-
\sigma_y^{(1)}\sigma_x^{(2)}
\right)
}
$$

with optional transverse Overhauser terms added if needed.

That is a very usable starting Hamiltonian.

# 4. Lindblad noise

Now add open-system dynamics:

$$
\dot\rho = -i[H(t),\rho] + \sum_k \mathcal D[L_k]\rho
$$

with

$$
\mathcal D[L]\rho = L\rho L^\dagger - \frac{1}{2}{L^\dagger L,\rho}.
$$

## 4.1 Spin relaxation

$$
L_1=\sqrt{\gamma_{1}^{(1)}},\sigma_-^{(1)},\qquad
L_2=\sqrt{\gamma_{1}^{(2)}},\sigma_-^{(2)}.
$$

These describe $T_1$-type decay.


## 4.2 Pure dephasing

$$
L_3=\sqrt{\gamma_{\phi}^{(1)}},\sigma_z^{(1)},\qquad
L_4=\sqrt{\gamma_{\phi}^{(2)}},\sigma_z^{(2)}.
$$

Depending on convention, some people use $\sqrt{\gamma_\phi/2}\sigma_z$. Just keep the convention consistent.


## 4.3 Correlated dephasing

If both spins see a common noise source:

$$
L_5=\sqrt{\gamma_c}\left(\sigma_z^{(1)}+\sigma_z^{(2)}\right).
$$

This is useful when modeling partially common Overhauser or charge-noise effects.


## 4.4 Optional charge-noise-induced exchange fluctuations

Instead of only Lindblad dephasing, you may include quasi-static or stochastic noise in (J):

$$
J(t)\to J(t)+\delta J(t).
$$

This is often more realistic than forcing all orbital noise into Lindblad form.


# 5. How the controls enter

Suppose the experimentally accessible control is a voltage (V(t)). Then typically

$$
J(t)=J(V(t)),
\qquad
D_z(t)=D_z(V(t)).
$$

A common phenomenological model is

$$
J(V)=J_0 e^{\beta (V-V_0)},
\qquad
D_z(V)=D_0 e^{\beta_D (V-V_0)}
$$

or sometimes

$$
D_z(V)=\eta_{\mathrm{SO}},J(V).
$$

That last relation is a very convenient simplified SOC model.

Then the Hamiltonian is entirely controlled by a single scalar voltage waveform $V(t)$.


# 6. Reachable controls: what can you do?

Now the key part.

## Case A: only exchange control (J(t))

Take

$$
H(t)=
\frac{\omega_1}{2}\sigma_z^{(1)}+\frac{\omega_2}{2}\sigma_z^{(2)}
+
\frac{J(t)}{4}\sigma_1\cdot \sigma_2.
$$

### If $\omega_1=\omega_2$ and there is no gradient

Then the only nontrivial control term is Heisenberg exchange.

Reachable unitary family is essentially

$$
U(T)=\exp\left[-i \theta, \mathbf S_1\cdot \mathbf S_2 \right]
$$

up to global phase, where

$$
\theta = \int_0^T J(t),dt.
$$

This gives:

* SWAP
* $\sqrt{\mathrm{SWAP}}$
* partial SWAP family

So:

* **entangling gates are reachable**
* **full SU(4) is not reachable**


### If $\omega_1\neq\omega_2$ or there is a static gradient

Then you have drift

$$
H_{\mathrm{grad}}=\frac{\Delta}{2}(\sigma_z^{(1)}-\sigma_z^{(2)})
$$

plus control $J(t)\sigma_1\cdot\sigma_2$.

This is much richer. In practice:

* exchange mixes $|\uparrow\downarrow\rangle$ and $|\downarrow\uparrow\rangle$
* gradient gives relative phase accumulation

That is enough to generate nontrivial entangling operations, and in many cases gives effective controllability on the relevant subspace.

Still, if the gradient is not itself controllable, you do not generally get arbitrary local single-qubit rotations. But you can often realize useful entangling gates much more flexibly than with exchange alone.


## Case B: exchange plus Rashba SOC, both voltage controlled

Now use

$$
H(t)=
H_{\mathrm{drift}}
+
u(t) H_J
+
v(t) H_D
$$

with

$$
H_J=\frac{1}{4}\sigma_1\cdot\sigma_2,
\qquad
H_D=\frac{1}{4}\left(\sigma_x^{(1)}\sigma_y^{(2)}-\sigma_y^{(1)}\sigma_x^{(2)}\right),
$$

and $u(t)=J(t)$, $v(t)=D_z(t)$.

This is much better for controllability, because the DM term does not commute with isotropic exchange in the same way. It breaks the symmetry and introduces effectively spin-dependent structure.

### Reachability improves to:

* arbitrary entangling gates become easier
* basis rotations within the two-spin manifold become possible
* with a Zeeman gradient present, full logical controllability is often achievable on the computational space

This is where all-electric control starts looking much closer to genuine spin control.


## Case C: exchange only, but encoded qubits

If you insist on no spin-dependent control, another route is encoding logical qubits into more than two spins.

With only exchange interactions between multiple spins, universal computation becomes possible. But with **just two physical spins**, exchange alone is not universal.


# 7. Lie-algebra intuition for reachability

A simple way to think about reachability is to ask what algebra is generated by your drift and control Hamiltonians.

## Exchange only

Generated algebra is too small:
$$
\mathrm{Lie}\{i \sigma_1\cdot\sigma_2\},
$$
just gives a one-parameter family.

## Exchange + gradient

Generated algebra:
$$
\mathrm{Lie}\left\{
i(\sigma_z^{(1)}-\sigma_z^{(2)}),
i\sigma_1\cdot\sigma_2
\right\}
$$
is much larger and gives useful entangling control.

## Exchange + DM + gradient

Generated algebra:
$$
\mathrm{Lie}\left\{
i(\sigma_z^{(1)}-\sigma_z^{(2)}),
i\sigma_1\cdot\sigma_2,
i(\sigma_x^{(1)}\sigma_y^{(2)}-\sigma_y^{(1)}\sigma_x^{(2)})
\right\}
$$
is larger still, and in practical terms this is the regime where electrically driven universal-like control becomes realistic.


# 8. A concrete control problem

Suppose you control only one voltage (V(t)), and it sets

$$
J(t)=J_0 e^{\beta V(t)},\qquad
D_z(t)=\eta_{\mathrm{SO}} J(t).
$$

Then your dynamics is

$$
\dot\rho = -i[H(V(t)),\rho] + \sum_k \mathcal D[L_k]\rho.
$$

A standard gate optimization is:

$$
\max_{V(t)} F_{\mathrm{avg}}\big(\mathcal E_T,\mathcal U_{\mathrm{target}}\big)
$$

for a target such as:

* $\sqrt{\mathrm{SWAP}}$
* iSWAP-like gate
* CZ-equivalent gate up to local phases

subject to:

* $V_{\min}\le V(t)\le V_{\max}$
* bandwidth constraints
* robustness to $\delta b_i$, $\delta J$, $\delta D_z$

This is a very standard open-system quantum control setup.


# 9. What is realistically reachable?

Here’s the practical summary.

## With only $J(t)$, no gradients, no SOC

Reachable:

* partial SWAP family
* $\sqrt{\mathrm{SWAP}}$
* Bell-state generation from suitable inputs

Not reachable:

* arbitrary two-qubit gate
* CNOT/CZ directly as a native control result


## With $J(t)$ and static Zeeman/Overhauser gradient

Reachable:

* richer entangling gates
* useful two-qubit logical control in relevant subspaces
* noise-adapted pulse shaping

Still limited if you cannot separately control local terms.


## With $J(t)$ plus Rashba-induced $D_z(t)$

Reachable:

* substantially richer electrically driven control
* effective spin-dependent manipulation using only voltages
* higher-quality entangling gates
* better chances for robust all-electric control

This is the most interesting two-spin all-electrical control regime.


# 10. Minimal model I recommend you start with

If you want one concrete model to work with numerically, use this:

$$
\boxed{
H(t)=
\frac{\Delta_z + \delta b_z}{2}(\sigma_z^{(1)}-\sigma_z^{(2)})
+
\frac{J(V(t))}{4}
\left(
\sigma_x^{(1)}\sigma_x^{(2)}+
\sigma_y^{(1)}\sigma_y^{(2)}+
\sigma_z^{(1)}\sigma_z^{(2)}
\right)
+
\frac{D(V(t))}{4}
\left(
\sigma_x^{(1)}\sigma_y^{(2)}-\sigma_y^{(1)}\sigma_x^{(2)}
\right)
}
$$

with Lindblad operators

$$
L_1=\sqrt{\gamma_1^{(1)}},\sigma_-^{(1)},\quad
L_2=\sqrt{\gamma_1^{(2)}},\sigma_-^{(2)},\quad
L_3=\sqrt{\gamma_\phi^{(1)}},\sigma_z^{(1)},\quad
L_4=\sqrt{\gamma_\phi^{(2)}},\sigma_z^{(2)}.
$$

This is compact, realistic, and already rich enough to study:

* exchange gates
* Overhauser robustness
* SOC-assisted all-electric control
* Lindblad-limited fidelity

