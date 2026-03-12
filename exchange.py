import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# 2) Atomic units:
# -------------------------
Eh = 27211.  # Eh = 27211 meV
Th = 2.42e-5 # Th = 2.42 x 10^-5 ps

# ----------------------------
# Parameters
# ----------------------------
J = .5/Eh       # (meV)
dt = .1/Th      # (ps)
n_steps = 200

# ----------------------------
# Hamiltonian S1.S2
# basis: |00>, |01>, |10>, |11>
# ----------------------------
H = np.array([
    [1,  0,  0, 0],
    [0, -1,  2, 0],
    [0,  2, -1, 0],
    [0,  0,  0, 1]
], dtype=complex)*J/4

# ----------------------------
# Initial state: |01>
# ----------------------------
psi0 = np.array([0, 1, 0, 0], dtype=complex)

# ----------------------------
# Initialization of the 3-point method
# we need psi_0 and psi_1
# psi_1 is obtained from one Euler step
# ----------------------------
psi_prev = psi0.copy()
psi_curr = psi0 - 1j * dt * (H @ psi0)

# auxiliary normalization
psi_curr /= np.linalg.norm(psi_curr)

# ----------------------------
# Simulation - iteration
# ----------------------------
times = [0.0, dt]
states = [psi_prev.copy(), psi_curr.copy()]

# ----------------------------
# #-point iteration
# psi_{n+1} = psi_{n-1} - 2i dt/hbar H psi_n
# ----------------------------
for n in range(1, n_steps):
    psi_next = psi_prev - 2j * dt * (H @ psi_curr)

    # normalization (H is hermitian but its discretization might be not!)
    psi_next /= np.linalg.norm(psi_next)

    psi_prev = psi_curr
    psi_curr = psi_next

    times.append((n + 1) * dt)
    states.append(psi_curr.copy())
times = np.array(times)
# ----------------------------
# Populations
# ----------------------------
states = np.array(states)
pop_00 = np.abs(states[:, 0])**2
pop_01 = np.abs(states[:, 1])**2
pop_10 = np.abs(states[:, 2])**2
pop_11 = np.abs(states[:, 3])**2

# ----------------------------
# Plot
# ----------------------------
fig, ax = plt.subplots()
ax.plot(times*Th, pop_00, label='|00>')
ax.plot(times*Th, pop_01, label='|01>')
ax.plot(times*Th, pop_10, label='|10>')
ax.plot(times*Th, pop_11, label='|11>')
ax.set_xlabel("t (ps)")
ax.set_ylabel("Populations")
ax.legend()
ax.grid(True)
fig.savefig('swap.png')   # save the figure to file
plt.close(fig) 
