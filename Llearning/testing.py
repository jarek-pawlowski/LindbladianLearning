import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
from torch.utils.data import DataLoader

from lindblad_solver import Lindblad_solver


import dataset as ds
import training as tr
from model import PINN


def _as_path(path: Path | str) -> Path:
    return Path(path).expanduser() if not isinstance(path, Path) else path


def load_test_subset(path: Path | str, device: torch.device = "cpu"):
    path = _as_path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing test subset file: {path}")
    return torch.load(path, map_location=device, weights_only=False)


def build_model(checkpoint_path: Path | str, feature_dim: int, device: torch.device = "cpu"):
    checkpoint_path = _as_path(checkpoint_path)
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Missing checkpoint: {checkpoint_path}")

    model = PINN(
        feature_dim=feature_dim,
        H=tr.H,
        L_ops=tr.L_ops,
        dt=tr.dt,
        hidden_dims=(3208, 802, 128, 64, 32),
    )
    state_dict = torch.load(checkpoint_path, map_location=device, weights_only=True)
    model.load_state_dict(state_dict)
    model.to(device)
    model.eval()
    return model


def plot_predictions(true_taus, pred_taus, save_path: Path | str | None = None):
    true_taus = true_taus.detach().cpu().numpy()
    pred_taus = pred_taus.detach().cpu().numpy()

    if save_path is not None:
        save_path = _as_path(save_path)

    n_params = true_taus.shape[1]
    fig, axes = plt.subplots(1, n_params, figsize=(4.2 * n_params, 4.2))
    if n_params == 1:
        axes = [axes]

    for idx, ax in enumerate(axes):
        param_true = true_taus[:, idx]
        param_pred = pred_taus[:, idx]

        min_val = min(param_true.min(), param_pred.min())
        max_val = max(param_true.max(), param_pred.max())
        padding = 0.05 * (max_val - min_val + 1e-12)
        axis_min = min_val - padding
        axis_max = max_val + padding

        ax.scatter(param_true, param_pred, s=20, color="tab:orange", alpha=0.8)
        ax.plot([axis_min, axis_max], [axis_min, axis_max],
                color="tab:blue", linestyle="--", linewidth=1.0, label="y = x")
        ax.set_xlim(axis_min, axis_max)
        ax.set_ylim(axis_min, axis_max)

        ax.set_xlabel("ground truth")
        ax.set_ylabel("prediction")
        ax.set_title(ds.TAU_KEYS[idx])
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best")
        ax.set_aspect('equal', adjustable='box')

    fig.suptitle("Predicted vs ground-truth parameters")
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    if save_path is not None:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=200)
        print(f"Saved plot to {save_path}")

    plt.close(fig)


def plot_rho_trajectory_comparison(model, test_ds, n_samples: int = 3, save_dir: Path | str = './results'):
    save_dir = _as_path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    dataset = test_ds.dataset
    indices = np.asarray(test_ds.indices)
    if len(indices) == 0:
        raise ValueError("The test subset is empty.")

    rng = np.random.default_rng(0)
    sample_indices = rng.choice(indices, size=min(n_samples, len(indices)), replace=False)

    model_device = next(model.parameters()).device

    saved_paths = []
    for sample_pos, sample_idx in enumerate(sample_indices):
        sample_idx = int(sample_idx)
        raw = dataset._load_raw(sample_idx)
        rhos_true = raw["rhos"]
        times = raw["times"]

        features, _, _, _ = dataset[sample_idx]
        features = features.unsqueeze(0).to(model_device)  # add batch dimension

        with torch.no_grad():
            _, taus_pred, _ = model(features)

        taus_pred = taus_pred.squeeze(0).cpu()
        gammas = (1.0 / taus_pred) * 1e12

        solver = Lindblad_solver(tr.H, tr.L_ops, gammas.tolist(), tr.dt)
        rho = torch.tensor(rhos_true[0].astype(np.complex64), dtype=torch.complex64, device=model_device)
        pred_rhos = [rho.detach().cpu().numpy()]

        for _ in range(len(rhos_true) - 1):
            rho = solver.step(rho)
            pred_rhos.append(rho.detach().cpu().numpy())

        pred_rhos = np.stack(pred_rhos, axis=0)

        fig, axes = plt.subplots(4, 4, figsize=(14, 10), sharex=True, constrained_layout=True)
        axes = np.atleast_2d(axes)

        for i in range(4):
            for j in range(4):
                ax = axes[i, j]
                gt_real = rhos_true[:, i, j].real
                pred_real = pred_rhos[:, i, j].real
                ax.plot(times, gt_real, label="ground truth (Re)", color="tab:blue", linewidth=1.3)
                ax.plot(times, pred_real, label="prediction (Re)", color="tab:orange", linewidth=1.2, alpha=0.9)

                if i != j:
                    gt_imag = rhos_true[:, i, j].imag
                    pred_imag = pred_rhos[:, i, j].imag
                    ax.plot(times, gt_imag, linestyle="--", color="tab:blue", alpha=0.6, label="ground truth (Im)")
                    ax.plot(times, pred_imag, linestyle="--", color="tab:orange", alpha=0.6, label="prediction (Im)")

                ax.set_title(rf"$\rho_{{{i}{j}}}$")
                ax.grid(True, alpha=0.3)
                if i == 3:
                    ax.set_xlabel("time")
                if j == 0:
                    ax.set_ylabel("value")
                if i == 0 and j == 0:
                    ax.legend(loc="best", fontsize="small")

        fig.suptitle(f"Density-matrix trajectory comparison (sample {sample_pos}, original idx {sample_idx})")
        out_path = save_dir / f"rho_comparison_sample_{sample_pos}.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        saved_paths.append(out_path)
        print(f"Saved rho comparison plot to {out_path}")

    return saved_paths


def evaluate_model(model, test_ds, batch_size: int = 64, num_samples: int = 3):
    loader = DataLoader(test_ds, batch_size=batch_size, shuffle=False, num_workers=0)

    all_true_taus = []
    all_pred_taus = []

    with torch.no_grad():
        for features, labels_log, _, _ in loader:
            features = features.to(next(model.parameters()).device)
            labels_log = labels_log.to(next(model.parameters()).device)
            _, taus_pred, _ = model(features)
            true_taus = torch.exp(labels_log)
            all_true_taus.append(true_taus.cpu())
            all_pred_taus.append(taus_pred.cpu())

    true_taus = torch.cat(all_true_taus, dim=0)
    pred_taus = torch.cat(all_pred_taus, dim=0)

    mae = torch.mean(torch.abs(pred_taus - true_taus)).item()
    mse = torch.mean((pred_taus - true_taus) ** 2).item()

    print("Overall test metrics:")
    print(f"  MAE: {mae:.6f}")
    print(f"  MSE: {mse:.6f}")

    print("\nPer-parameter metrics:")
    for idx, name in enumerate(ds.TAU_KEYS):
        param_true = true_taus[:, idx]
        param_pred = pred_taus[:, idx]
        mae_i = torch.mean(torch.abs(param_pred - param_true)).item()
        mse_i = torch.mean((param_pred - param_true) ** 2).item()
        print(f"  {name:>12s}: MAE={mae_i:.6f}, MSE={mse_i:.6f}")

    print("\nExample predictions:")
    for i in range(min(num_samples, len(true_taus))):
        print(f"  sample {i}")
        print(f"    true: {true_taus[i].tolist()}")
        print(f"    pred: {pred_taus[i].tolist()}")

    return {
        "true_taus": true_taus,
        "pred_taus": pred_taus,
        "mae": mae,
        "mse": mse,
    }
