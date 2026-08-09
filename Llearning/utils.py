import numpy as np
import torch
import matplotlib.pyplot as plt


# Plotting
def plot_density_matrix_elements_t(times, rhos_t, filename=None, show=True):
    """Plot every element of a density-matrix trajectory.

    Diagonal elements are plotted using only their real parts. For every
    off-diagonal element, both the real and imaginary parts are shown.

    """

    dim = rhos_t[0].shape[-1]
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
            for k, rho_t in enumerate(rhos_t):
                rho_ij = rho_t[:, i, j]
                element_label = rf"$\rho_{{{basis_labels[i]},{basis_labels[j]}}}$"
                if k == 0:
                    ax.plot(times, rho_ij.real, label="Re")
                else:
                    ax.plot(times, rho_ij.real)
                if i != j:
                    if k == 0:
                        ax.plot(times, rho_ij.imag, linestyle="--", label="Im")
                    else:
                        ax.plot(times, rho_ij.imag, linestyle="--")
                        
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