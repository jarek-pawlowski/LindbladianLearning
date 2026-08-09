import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import matplotlib.pyplot as plt

# Atomic units
Eh = 27211.  # Eh = 27211 meV
Th = 2.42e-5 # Th = 2.42 x 10^-5 ps

# Dwie funkcje pożyczone od Ziemka

def plot_signals(signals, t_axis = None, labels=None, points=None, pointlabels=None):
    fig, ax = plt.subplots()
    for s, sig in enumerate(signals):
        if t_axis is None:
            t_axis = np.linspace(0.,1., num=len(sig), endpoint=True)
        ax.plot(t_axis, sig, label=labels[s])
    if points is not None:
        for p, pts in enumerate(points):
            ax.scatter(pts[0], pts[1], marker='x', label=pointlabels[p])
    if labels is not None or pointlabels is not None: ax.legend()
    plt.show()

def FixRho(rho):
    # Handle both 2D (N,N) and 3D batched (B,N,N) inputs
    squeeze = rho.ndim == 2
    if squeeze:
        rho = rho.unsqueeze(0)

    # 1. Force Hermitian (more aggressively)
    rho_dag = rho.conj().transpose(1, 2)
    rho = (rho + rho_dag) / 2
    # Add tiny regularization to diagonal to help convergence
    eps = 1e-10
    eye = torch.eye(rho.shape[-1], dtype=rho.dtype, device=rho.device)
    rho = rho + eps * eye.unsqueeze(0)

    # 2. Batched eigen-decomposition
    eigvals, eigvecs = torch.linalg.eigh(rho)

    # 3. Clip negative eigenvalues
    eigvals = torch.clamp(eigvals, min=0)

    # 4. Reconstruct:  V @ diag(λ) @ V†
    rho_fixed = eigvecs @ torch.diag_embed(eigvals.to(eigvecs.dtype)) @ eigvecs.conj().transpose(1, 2)
    
    # 5. Normalize trace
    trace = rho_fixed.diagonal(dim1=1, dim2=2).sum(dim=1)
    rho_fixed = rho_fixed / trace[:, None, None]

    return rho_fixed.squeeze(0) if squeeze else rho_fixed


"""LINDBLAD SOLVER FOR ARBITRARY NUMBER OF QUBITS"""
class Lindblad_solver:
    """
    Solves the Lindblad master equation for the density matrix ρ:

        dρ/dt = -i[H, ρ] + Σ_k γ_k D[L_k](ρ)

    where D[L](ρ) = L ρ L† - ½ {L†L, ρ}  (dissipator)
    """

    def __init__(self, H, L_ops, gammas, dt):
        """
        H      : Hamiltonian (NxN complex array)
        L_ops  : list of jump operators [L_k], each (NxN)
        gammas : list of decay rates [g_k], one per jump operator
        dt     : time step
        """
        self.H      = torch.tensor(np.array(H, dtype=np.complex64))  # complex64, not complex128
        self.L_ops  = [torch.tensor(np.array(L, dtype=np.complex64)) for L in L_ops]
        self.gammas = list(gammas)
        self.dt     = dt
        self.LdagL  = [L.conj().mT @ L for L in self.L_ops]

    def lindblad_rhs(self, rho):
        H = self.H.to(device=rho.device, dtype=rho.dtype)
        drho = -1j * (H @ rho - rho @ H)

        for gamma, Lk, LdLk in zip(self.gammas, self.L_ops, self.LdagL):
            Lk = Lk.to(device=rho.device, dtype=rho.dtype)
            LdLk = LdLk.to(device=rho.device, dtype=rho.dtype)
            drho += gamma * (Lk @ rho @ Lk.conj().mT - 0.5 * (LdLk @ rho + rho @ LdLk))
        return drho

    def _renormalize(self, rho):
        """Enforce Tr(ρ) = 1 after each step to suppress floating-point drift."""
        return rho / np.trace(rho)

    # RK4 step

    def step(self, rho):
        """
        Advance ρ by one time step using classical RK4.
        Returns rho_{n+1}.
        """
        L  = self.lindblad_rhs

        k1 = L(rho)
        k2 = L(rho + 0.5 * self.dt * k1)
        k3 = L(rho + 0.5 * self.dt * k2)
        k4 = L(rho +       self.dt * k3)

        rho_next = rho + (self.dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        return FixRho(rho_next)

    # Full trajectory

    def run_simulation(self, rho0, n_steps=1000):
        """
        Propagate from rho0 for n_steps time steps.
        Returns array of shape (n_steps+1, N, N).
        """
        rho = rho0.clone()
        rhos = [rho.clone()]
        times = [0]

        for i in range(n_steps):
            rho = self.step(rho)
            rhos.append(rho.clone())
            times.append((i + 1) * dt)
            #print(np.amin(np.linalg.eigvals(rho.cpu().numpy())))

        return [np.array(times), np.array(rhos)]

    def populations(rhos):
        """Extract diagonal populations from a trajectory array."""
        return np.real(np.einsum('tii->ti', rhos))   # shape (T, N)


# Plotting
def plot_density_matrix_elements(times, rhos, filename=None, show=True):
    """Plot every element of a density-matrix trajectory.

    Diagonal elements are plotted using only their real parts. For every
    off-diagonal element, both the real and imaginary parts are shown.

    Parameters
    ----------
    times : array-like, shape (T,)
        Time axis.
    rhos : array-like or torch.Tensor
        Density-matrix trajectory with shape (T, N, N) or (T, 1, N, N).
    filename : str or pathlib.Path, optional
        Output path used to save the figure.
    show : bool, default=True
        Whether to display the figure with ``plt.show()``.
    """
    if isinstance(rhos, torch.Tensor):
        rhos_np = rhos.detach().cpu().numpy()
    else:
        rhos_np = np.asarray(rhos)

    # Remove a singleton batch dimension: (T, 1, N, N) -> (T, N, N).
    if rhos_np.ndim == 4 and rhos_np.shape[1] == 1:
        rhos_np = rhos_np[:, 0]

    if rhos_np.ndim != 3 or rhos_np.shape[1] != rhos_np.shape[2]:
        raise ValueError(
            "rhos must have shape (T, N, N) or (T, 1, N, N); "
            f"received {rhos_np.shape}."
        )

    times = np.asarray(times)
    if len(times) != len(rhos_np):
        raise ValueError(
            f"times and rhos must have the same trajectory length; "
            f"received {len(times)} and {len(rhos_np)}."
        )

    dim = rhos_np.shape[-1]
    basis_width = max(1, int(np.ceil(np.log2(dim))))
    basis_labels = [format(i, f"0{basis_width}b") for i in range(dim)]

    fig, axes = plt.subplots(
        dim,
        dim,
        figsize=(4.0 * dim, 3.0 * dim),
        sharex=True,
        constrained_layout=True,
    )
    axes = np.atleast_2d(axes)

    for i in range(dim):
        for j in range(dim):
            ax = axes[i, j]
            rho_ij = rhos_np[:, i, j]
            element_label = rf"$\rho_{{{basis_labels[i]},{basis_labels[j]}}}$"

            ax.plot(times, rho_ij.real, label="Re")
            if i != j:
                ax.plot(times, rho_ij.imag, linestyle="--", label="Im")

            ax.set_title(element_label)
            ax.grid(True, alpha=0.3)

            if i == dim - 1:
                ax.set_xlabel("time (ps)")
            if j == 0:
                ax.set_ylabel("value")
            ax.legend(loc="best", fontsize="small")

    fig.suptitle(r"Two-qubit density matrix $\rho(t)$", fontsize=16)
    
    if filename is not None:
        fig.savefig(filename, dpi=200, bbox_inches="tight")
    if show:
        plt.show()

    return fig, axes


"""NUMERICAL SOLUTIONS"""

# 2 qubit - dissipation

J = .5/Eh       # (meV)
dt = .1/Th      # (ps)
#dt = .02/Th      # (ps)
n_steps = 200

H = np.array([
    [1,  0,  0, 0],
    [0, -1,  2, 0],
    [0,  2, -1, 0],
    [0,  0,  0, 1]
], dtype=complex)*J/4
I2 = np.eye(2, dtype=complex)
L1 = np.array([[0, 1],
               [0, 0]], dtype=complex)
L2 = np.array([[1,  0],
               [0, -1]], dtype=complex)
L1 = np.kron(I2, L1)
L2 = np.kron(L2, I2)
L = [L1, L2]

gamma_relax = 0.05 * J   # amplitude damping  (T1-type, relaxation to ground state)
gamma_deph  = 0.02 * J   # dephasing          (T2-type, pure dephasing)
gammas = [gamma_relax, gamma_deph]

psi0 = np.array([1, 0, 0, 1], dtype=complex)/np.sqrt(2.)
psi0 /= np.linalg.norm(psi0)
psi0 = torch.tensor(psi0, dtype=torch.complex64)
rho0 = torch.outer(psi0, psi0.conj()).unsqueeze(0)

solver_3 = Lindblad_solver(H, L_ops=L, gammas=gammas, dt=dt)
[times_3, rhos_3]  = solver_3.run_simulation(rho0, n_steps=n_steps)
times_3 *= Th  # convert to ps

# Plotting
plot_density_matrix_elements(
    times_3,
    rhos_3,
    filename="lindbladians.png",
)