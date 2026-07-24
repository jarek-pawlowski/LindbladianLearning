"""Reconstruct the two-photon polarization density matrix with NumPy/SciPy.

The implementation uses the augmented A-S-b formalism

    x_dot   = A(t) x
    eta_dot = (I_4 kron A(t)) eta + S_stack x
    m_dot   = (I_4 kron B_read) eta

where x = vec(rho) uses column-major (Fortran) vectorization.  The four
auxiliary vectors eta correspond to first-emission source labels
(j,l) = HH, HV, VH, VV.  For each source, four readouts (k,m) are accumulated,
so m has 16 complex components.

The physical parameters and atomic-unit conversions match cascade.py.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

ComplexArray = NDArray[np.complex128]
RealArray = NDArray[np.float64]


@dataclass(frozen=True)
class Parameters:
    # Atomic-unit conversions copied from cascade(3).py
    Eh_meV: float = 27211.0       # 1 Hartree in meV
    Th_ps: float = 2.42e-5        # 1 atomic time unit in ps
    hbar: float = 1.0

    # Physical parameters
    delta_meV: float = 0.001
    EB_meV: float = 4.0
    gammaX_per_ps: float = 0.005
    gammaB_per_ps: float = 0.010

    # Pulse parameters
    tL_ps: float = 100.0
    FWHM_ps: float = 20.0
    alphaH: float = 1.0
    alphaV: float = 0.0

    # Integration window and output sampling
    post_pulse_lifetimes: float = 8.0
    n_output: int = 1200
    rtol: float = 2.0e-8
    atol: float = 2.0e-10

    @property
    def delta(self) -> float:
        return self.delta_meV / self.Eh_meV

    @property
    def EB(self) -> float:
        return self.EB_meV / self.Eh_meV

    @property
    def gammaX(self) -> float:
        # A rate in ps^{-1} is multiplied by Th_ps to obtain atomic units.
        return self.gammaX_per_ps * self.Th_ps

    @property
    def gammaB(self) -> float:
        return self.gammaB_per_ps * self.Th_ps

    @property
    def tL(self) -> float:
        return self.tL_ps / self.Th_ps

    @property
    def FWHM(self) -> float:
        return self.FWHM_ps / self.Th_ps

    @property
    def Delta_XL(self) -> float:
        return self.EB / 2.0

    @property
    def theta(self) -> float:
        # Same parametrization as cascade(3).py
        denominator = np.sqrt(2.0 * np.pi * np.log(2.0))
        return np.pi * np.sqrt(self.EB * self.FWHM / denominator)

    @property
    def t_final(self) -> float:
        # The slowest decay controls convergence.  This corrects the original
        # use of max(gammaX,gammaB), which can truncate the exciton tail.
        slowest_rate = min(self.gammaX, self.gammaB)
        pulse_end = self.tL + 4.0 * self.FWHM
        return pulse_end + self.post_pulse_lifetimes / slowest_rate


def ket(index: int, dimension: int = 4) -> ComplexArray:
    v = np.zeros(dimension, dtype=complex)
    v[index] = 1.0
    return v


def projector(v: ComplexArray) -> ComplexArray:
    return np.outer(v, v.conj())


def vec(matrix: ComplexArray) -> ComplexArray:
    """Column-major vectorization compatible with vec(P rho Q)."""
    return np.asarray(matrix, dtype=complex).reshape(-1, order="F")


def unvec(vector: ComplexArray, dimension: int = 4) -> ComplexArray:
    return np.asarray(vector, dtype=complex).reshape((dimension, dimension), order="F")


def omega(t: float | RealArray, p: Parameters) -> float | RealArray:
    arg = 4.0 * np.log(2.0) * ((t - p.tL) / p.FWHM) ** 2
    prefactor = np.sqrt(4.0 * np.log(2.0) / np.pi) * p.theta / p.FWHM
    return prefactor * np.exp(-arg)


def build_operators(p: Parameters) -> dict[str, ComplexArray | list[ComplexArray]]:
    G, XH, XV, B = (ket(i) for i in range(4))

    PG = projector(G)
    PXH = projector(XH)
    PXV = projector(XV)
    PB = projector(B)

    sigma_H = np.outer(G, XH.conj()) + np.outer(XH, B.conj())
    sigma_V = np.outer(G, XV.conj()) + np.outer(XV, B.conj())
    sigma_L = p.alphaH * sigma_H + p.alphaV * sigma_V

    # Correct rotating-frame Hamiltonian.  At exact two-photon resonance,
    # 2 Delta_XL - EB = 0, so the biexciton diagonal entry is zero.
    H0 = (
        (p.Delta_XL + p.delta / 2.0) * PXH
        + (p.Delta_XL - p.delta / 2.0) * PXV
        + (2.0 * p.Delta_XL - p.EB) * PB
    )
    H_drive = -(p.hbar / 2.0) * (sigma_L + sigma_L.conj().T)

    collapse_ops = [
        np.sqrt(p.gammaX) * np.outer(G, XH.conj()),
        np.sqrt(p.gammaX) * np.outer(G, XV.conj()),
        np.sqrt(p.gammaB / 2.0) * np.outer(XH, B.conj()),
        np.sqrt(p.gammaB / 2.0) * np.outer(XV, B.conj()),
    ]

    return {
        "PG": PG,
        "H0": H0,
        "H_drive": H_drive,
        "sigma_H": sigma_H,
        "sigma_V": sigma_V,
        "collapse_ops": collapse_ops,
    }


def liouvillian(H: ComplexArray, collapse_ops: list[ComplexArray], hbar: float = 1.0) -> ComplexArray:
    """Return A such that d vec(rho)/dt = A vec(rho), using column stacking."""
    d = H.shape[0]
    identity = np.eye(d, dtype=complex)

    # vec(H rho) = (I kron H) vec(rho)
    # vec(rho H) = (H^T kron I) vec(rho)
    A = -(1j / hbar) * (np.kron(identity, H) - np.kron(H.T, identity))

    for L in collapse_ops:
        LdL = L.conj().T @ L
        A += (
            np.kron(L.conj(), L)
            - 0.5 * np.kron(identity, LdL)
            - 0.5 * np.kron(LdL.T, identity)
        )
    return A


def build_ASb(
    p: Parameters,
) -> tuple[
    ComplexArray,
    ComplexArray,
    Callable[[float], ComplexArray],
    ComplexArray,
    ComplexArray,
    list[tuple[str, str]],
]:
    """Construct A0, A_drive, A(t), stacked S, and readout matrix B."""
    ops = build_operators(p)
    H0 = np.asarray(ops["H0"])
    H_drive = np.asarray(ops["H_drive"])
    collapse_ops = list(ops["collapse_ops"])

    A0 = liouvillian(H0, collapse_ops, p.hbar)
    # The pulse Liouvillian is purely Hamiltonian; no dissipator is duplicated.
    A_drive = liouvillian(H_drive, [], p.hbar)

    def A_of_t(t: float) -> ComplexArray:
        return A0 + omega(t, p) * A_drive

    sigma = {
        "H": np.asarray(ops["sigma_H"]),
        "V": np.asarray(ops["sigma_V"]),
    }
    pair_basis = [("H", "H"), ("H", "V"), ("V", "H"), ("V", "V")]

    # Source label a=(j,l): vec(sigma_l rho sigma_j^dagger)
    # = (sigma_j^* kron sigma_l) vec(rho).
    source_blocks = []
    for j, ell in pair_basis:
        source_blocks.append(np.kron(sigma[j].conj(), sigma[ell]))
    S_stack = np.vstack(source_blocks)  # shape (64,16)

    # Readout label r=(k,m): Tr(sigma_k^dagger sigma_m eta) = b_km^T vec(eta).
    readout_rows = []
    for k, m in pair_basis:
        observable = sigma[k].conj().T @ sigma[m]
        readout_rows.append(vec(observable.T))
    B_read = np.vstack(readout_rows)  # shape (4,16)

    return A0, A_drive, A_of_t, S_stack, B_read, pair_basis


def augmented_rhs(
    t: float,
    z: ComplexArray,
    A_of_t: Callable[[float], ComplexArray],
    S_stack: ComplexArray,
    B_read: ComplexArray,
) -> ComplexArray:
    """RHS for z=(x, eta_HH, eta_HV, eta_VH, eta_VV, m_16)."""
    x = z[:16]
    eta = z[16:80].reshape((4, 16))

    A = A_of_t(t)
    dx = A @ x

    # Each row eta[a] is stored as an ordinary 16-vector.  Right multiplication
    # by A.T applies A to every row vector: (A @ eta[a]).
    source = (S_stack @ x).reshape((4, 16))
    deta = eta @ A.T + source

    # dm[a,r] = b_r^T eta_a.  Shape is (4 sources, 4 readouts).
    dm_matrix = eta @ B_read.T

    return np.concatenate((dx, deta.reshape(-1), dm_matrix.reshape(-1)))


def raw_moments_to_matrix(
    m_vector: ComplexArray,
    pair_basis: list[tuple[str, str]],
) -> ComplexArray:
    """Assemble R[row=(l,m), col=(j,k)] from m[source=(j,l), readout=(k,m)]."""
    m_by_source_readout = m_vector.reshape((4, 4))
    pair_to_index = {pair: i for i, pair in enumerate(pair_basis)}
    raw = np.zeros((4, 4), dtype=complex)

    for source_index, (j, ell) in enumerate(pair_basis):
        for readout_index, (k, m) in enumerate(pair_basis):
            row = pair_to_index[(ell, m)]
            col = pair_to_index[(j, k)]
            raw[row, col] = m_by_source_readout[source_index, readout_index]
    return raw


def make_physical_density_matrix(raw: ComplexArray, clip_negative: bool = True) -> ComplexArray:
    """Hermitize, normalize, and optionally remove tiny numerical negativity."""
    rho = 0.5 * (raw + raw.conj().T)
    trace = np.trace(rho)
    if abs(trace) < 1.0e-14:
        raise RuntimeError("The integrated coincidence matrix has nearly zero trace.")
    rho /= trace

    if clip_negative:
        eigenvalues, eigenvectors = np.linalg.eigh(rho)
        if np.min(eigenvalues) < -1.0e-7:
            print(
                "Warning: appreciable negative eigenvalue before projection:",
                float(np.min(eigenvalues)),
            )
        eigenvalues = np.clip(eigenvalues.real, 0.0, None)
        if eigenvalues.sum() == 0.0:
            raise RuntimeError("Positivity projection removed the entire matrix.")
        rho = (eigenvectors * eigenvalues) @ eigenvectors.conj().T
        rho /= np.trace(rho)

    return rho


def concurrence(rho: ComplexArray) -> float:
    """Wootters concurrence for the basis |HH>, |HV>, |VH>, |VV>."""
    sigma_y = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
    spin_flip = np.kron(sigma_y, sigma_y)
    R = rho @ spin_flip @ rho.conj() @ spin_flip
    eigenvalues = np.linalg.eigvals(R)
    roots = np.sqrt(np.clip(np.real_if_close(eigenvalues).real, 0.0, None))
    roots = np.sort(roots)[::-1]
    return float(max(0.0, roots[0] - roots[1] - roots[2] - roots[3]))


def reconstruct_photonic_density_matrix(
    p: Parameters = Parameters(),
) -> tuple[ComplexArray, ComplexArray, float, solve_ivp]:
    """Run one augmented forward solve and reconstruct rho_2p."""
    ops = build_operators(p)
    _, _, A_of_t, S_stack, B_read, pair_basis = build_ASb(p)

    x0 = vec(np.asarray(ops["PG"]))
    z0 = np.zeros(96, dtype=complex)
    z0[:16] = x0

    t_eval = np.linspace(0.0, p.t_final, p.n_output)
    solution = solve_ivp(
        fun=lambda t, z: augmented_rhs(t, z, A_of_t, S_stack, B_read),
        t_span=(0.0, p.t_final),
        y0=z0,
        t_eval=t_eval,
        method="DOP853",
        rtol=p.rtol,
        atol=p.atol,
    )
    if not solution.success:
        raise RuntimeError(f"ODE solver failed: {solution.message}")

    m_final = solution.y[80:96, -1]
    raw = raw_moments_to_matrix(m_final, pair_basis)
    rho_2p = make_physical_density_matrix(raw)
    C = concurrence(rho_2p)
    return rho_2p, raw, C, solution


def print_matrix(name: str, matrix: ComplexArray) -> None:
    print(f"{name} =")
    with np.printoptions(precision=7, suppress=True, linewidth=160):
        print(matrix)


def plot_electronic_density_matrix(
    solution: solve_ivp,
    p: Parameters,
    filename: str = "biexciton.png",
    show: bool = False,
) -> None:
    """Plot the pulse and selected electronic density-matrix elements.

    The plot follows cascade.py and displays Omega(t), rho_GG(t),
    rho_BB(t), and |rho_BG(t)| as functions of time.
    """
    time_ps = solution.t * p.Th_ps

    # The first 16 augmented variables are x(t) = vec(rho(t)).
    x_history = solution.y[:16, :]
    rho_history = np.stack(
        [unvec(x_history[:, i]) for i in range(x_history.shape[1])],
        axis=0,
    )

    rho_GG = rho_history[:, 0, 0].real
    rho_BB = rho_history[:, 3, 3].real
    rho_BG_abs = np.abs(rho_history[:, 3, 0])

    # Convert the pulse from inverse atomic time to ps^{-1}.
    pulse_per_ps = omega(solution.t, p) / p.Th_ps

    fig, ax = plt.subplots()
    ax.plot(time_ps, pulse_per_ps, label=r"$\Omega(t)$")
    ax.plot(time_ps, rho_GG, label=r"$\rho_{GG}(t)$")
    ax.plot(time_ps, rho_BB, label=r"$\rho_{BB}(t)$")
    ax.plot(time_ps, rho_BG_abs, label=r"$|\rho_{BG}(t)|$")

    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Population / pulse amplitude")
    ax.legend()
    fig.tight_layout()
    fig.savefig(filename, dpi=200)

    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    p = Parameters()
    rho_2p, raw, C, solution = reconstruct_photonic_density_matrix(p)
    plot_electronic_density_matrix(solution, p, filename="biexciton.png")

    print(f"Theta = {p.theta:.10g}")
    print(f"t_final = {p.t_final * p.Th_ps:.3f} ps")
    print(f"ODE evaluations = {solution.nfev}")
    print_matrix("Raw integrated pair matrix", raw)
    print_matrix("rho_2p in [HH, HV, VH, VV]", rho_2p)
    print("trace(rho_2p) =", np.trace(rho_2p))
    print("eigenvalues(rho_2p) =", np.linalg.eigvalsh(rho_2p))
    print("Concurrence =", C)


if __name__ == "__main__":
    main()
