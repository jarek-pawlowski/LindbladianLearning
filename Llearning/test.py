from pathlib import Path

import torch
import testing as ts


DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

results_dir = Path("./results")
results_dir.mkdir(parents=True, exist_ok=True)

checkpoint_path = "./best_pinn_dataset1.pt"
test_subset_path = "./test_ds.pt"

test_ds = ts.load_test_subset(test_subset_path, device=DEVICE)
feature_dim = test_ds.dataset.feature_dim
model = ts.build_model(checkpoint_path, feature_dim, device=DEVICE)

results = ts.evaluate_model(model, test_ds)
ts.plot_predictions(results["true_taus"], results["pred_taus"], save_path=results_dir / "parameter_predictions.png")
ts.plot_rho_trajectory_comparison(model, test_ds, n_samples=3, save_dir=results_dir)
