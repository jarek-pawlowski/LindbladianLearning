import numpy as np
from scipy.linalg import expm

# ---------- basic operators ----------
I2 = np.eye(2, dtype=complex)

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
sm = np.array([[0, 0], [1, 0]], dtype=complex)
sx1, sy1, sz1 = np.kron(sx, I2), np.kron(sy, I2), np.kron(sz, I2)
sx2, sy2, sz2 = np.kron(I2, sx), np.kron(I2, sy), np.kron(I2, sz)
sm1, sm2 = np.kron(sm, I2), np.kron(I2, sm)
I4 = np.eye(4, dtype=complex)

# ---------- vectorization ----------
def vec(rho):
    return rho.reshape(-1, order="F")

def mat(v):
    return v.reshape((4, 4), order="F")

# ---------- voltage-to-coupling models ----------
def J_of_V(V, J0=1.0, beta=1.0, V0=0.0):
    return J0 * np.exp(beta * (V - V0))

def D_of_V(V, eta_so=0.05, J0=1.0, beta=1.0, V0=0.0):
    return eta_so * J_of_V(V, J0, beta, V0)

# ---------- Hamiltonian ----------
def H_two_spin(
    V,
    delta_bz=0.0,
    J0=1.0,
    beta=1.0,
    V0=0.0,
    eta_so=0.05,
):
    J = J_of_V(V, J0, beta, V0)
    D = D_of_V(V, eta_so, J0, beta, V0)

    H_grad = 0.5 * delta_bz * (sz1 - sz2)

    H_ex = (J / 4.0) * (
        sx1 @ sx2 + sy1 @ sy2 + sz1 @ sz2
    )

    H_soc = (D / 4.0) * (
        sx1 @ sy2 - sy1 @ sx2
    )

    return H_grad + H_ex + H_soc

# ---------- Lindblad superoperator ----------
def lindblad_superoperator(L):
    LdL = L.conj().T @ L
    return (
        np.kron(L.conj(), L)
        - 0.5 * np.kron(I4, LdL)
        - 0.5 * np.kron(LdL.T, I4)
    )

def liouvillian(H, collapse_ops):
    # -i[H,rho]
    L_unitary = -1j * (
        np.kron(I4, H) - np.kron(H.T, I4)
    )

    L_diss = np.zeros((16, 16), dtype=complex)
    for c in collapse_ops:
        L_diss += lindblad_superoperator(c)

    return L_unitary + L_diss

'''
# Alternatively we could avoid vectorization, then the master equation would be:
# But this is less stable for numerical integration (needs RK4 or similar).
def drho_dt(rho, H, Ls):
    comm = -1j * (H @ rho - rho @ H)

    dissipator = 0
    for L in Ls:
        dissipator += (
            L @ rho @ L.conj().T
            - 0.5 * (L.conj().T @ L @ rho)
            - 0.5 * (rho @ L.conj().T @ L)
        )
    return comm + dissipator
'''

# ---------- evolution ----------
def evolve_lindblad(
    rho0,
    dt, Nt,
    V_of_t,
    params=None,
    rates=None,
):
    if params is None:
        params = {}
    if rates is None:
        rates = {}

    gamma1_1 = rates.get("gamma1_1", 0.0)
    gamma1_2 = rates.get("gamma1_2", 0.0)
    gamma_phi_1 = rates.get("gamma_phi_1", 0.0)
    gamma_phi_2 = rates.get("gamma_phi_2", 0.0)
    gamma_c = rates.get("gamma_c", 0.0)  # correlated dephasing 

    collapse_ops = [
        np.sqrt(gamma1_1) * sm1,
        np.sqrt(gamma1_2) * sm2,
        np.sqrt(gamma_phi_1) * sz1,
        np.sqrt(gamma_phi_2) * sz2,
        np.sqrt(gamma_c) * (sz1 + sz2),
    ]

    collapse_ops = [c for c in collapse_ops if np.linalg.norm(c) > 0]

    rhos = [rho0.copy()]
    rho_vec_0 = vec(rho0)
    
    # method starting:
    H = H_two_spin(V_of_t(0.), **params)
    L = liouvillian(H, collapse_ops)
    rho_vec_1 = expm(L * dt) @ rho_vec_0
    rho_vec_1 /= np.linalg.norm(rho_vec_1)
    rho = mat(rho_vec_1)
    # numerical cleanup
    rho = 0.5 * (rho + rho.conj().T)
    rho = rho / np.trace(rho)
    rhos.append(rho)
    
    for n in range(Nt-1):
        H = H_two_spin(V_of_t((n+1)*dt), **params)
        L = liouvillian(H, collapse_ops)
        rho_vec_2 = rho_vec_0 + 2*dt*(L @ rho_vec_1)
        rho_vec_2 /= np.linalg.norm(rho_vec_2)
        
        rho_vec_0 = rho_vec_1
        rho_vec_1 = rho_vec_2
        
        rho = mat(rho_vec_1)
        # numerical cleanup
        rho = 0.5 * (rho + rho.conj().T)
        rho = rho / np.trace(rho)
        rhos.append(rho)

    return np.array(rhos)

# Example evolution
# initial state |up down><up down|
up = np.array([1, 0], dtype=complex)
down = np.array([0, 1], dtype=complex)
psi0 = np.kron(up, down)
rho0 = np.outer(psi0, psi0.conj())

# custom voltage pulse
def V_of_t(t):
    # smooth pulse example
    T = 10.0
    return 0.3*np.sin(np.pi * t / T)**2 if 0 <= t <= T else 0.0

dt = 0.02
Nt = 1000

# simple exchange
params = dict(
    delta_bz=0.0,
    J0=1.0,
    beta=3.0,
    V0=0.1,
    eta_so=0.0,
)

rates = dict(
    gamma1_1=0.1,
    gamma1_2=0.1,
    gamma_phi_1=0.05,
    gamma_phi_2=0.05,
    gamma_c=0.000,
)

'''
params = dict(
    delta_bz=0.1,
    J0=1.0,
    beta=2.0,
    V0=0.0,
    eta_so=0.05,
)

rates = dict(
    gamma1_1=0.001,
    gamma1_2=0.001,
    gamma_phi_1=0.005,
    gamma_phi_2=0.005,
    gamma_c=0.000,
)
'''

rhos = evolve_lindblad(rho0, dt, Nt, V_of_t, params, rates)
rho_final = rhos[-1]

print("Final density matrix:")
print(np.round(rho_final, 4))

print("Trace:", np.trace(rho_final))
print("Purity:", np.real(np.trace(rho_final @ rho_final)))

spins = np.array([[np.trace(rho@sz1).real/2., np.trace(rho@sz2).real/2.] for rho in rhos])
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
times = np.linspace(0, dt*Nt, Nt+1)
ax.plot(times, spins[:, 0], label='<sz1>')
ax.plot(times, spins[:, 1], label='<sz2>')
ax.plot(times, [V_of_t(t) for t in times], label='V(t)')
ax.set_xlabel("t (ps)")
ax.set_ylabel("spins")
ax.legend()
ax.grid(True)
fig.savefig('spins.png')
plt.close(fig) 