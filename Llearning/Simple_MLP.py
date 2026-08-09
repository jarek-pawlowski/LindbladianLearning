'''
Simple MPL
Autor: Maria Hadam
28.06.2026
Edited: Ziemowit Olinkiewicz
20.07.2026
'''

import json
import torch
import sys
sys.path.append("./source")

import model as mlp
import numpy as np
import lindblad_solver as solver
import dataset as ds
import training as t
import plotting as plt
import utils

# Load the dataset
root = "data/dataset_1_rel20_decoh40_5000t"     # name of the dataset directory
# root = "data/dataset_2_rel40_decoh20_5000t"

dataset = ds.TrajectoryDataset(root)
x, y, rho_t, rho_tp1 = dataset[0]                # load the first trajectory
print(f"number of trajectories: {len(dataset)}")
print(f"feature shape (input): {x.shape}  (n_steps+1={dataset._n_steps_plus_1}, 32 per step)")
print(f"label shape: {y.shape}  -> {ds.TAU_KEYS}")
print(f"label of the first trajectory (log-tau): {y}")
print(f"rho_t shape: {rho_t.shape}, rho_tp1 shape: {rho_tp1.shape}  (random consecutive pair)")

# Preprocessing: prepare the dataset for training
train_ds, val_ds, test_ds = ds.split_dataset(dataset)   # split into train/val/test
print(f"split sizes: train={len(train_ds)} val={len(val_ds)} test={len(test_ds)}")

mean, std, labels = ds.compute_label_stats(dataset, train_ds.indices)
print(f"train label mean: {mean}")
print(f"train label std:  {std}")

# Visualize training flows
rhos = []
for i in range(len(dataset)):
    rhos.append(dataset._load_raw(i)["rhos"])
times = dataset._load_raw(0)["times"]
utils.plot_density_matrix_elements_t(times, rhos[:100], filename="./results/train_density_matrices_flow.png")

# Test the MLP model
input_dim = dataset.feature_dim
output_dim = 5
model = mlp.MLP(input_dim, output_dim, hidden_dims=[256, 128, 64, 32], activation='relu', use_batchnorm=True, dropout_rate=0.2)
print(model)

dummy_x = torch.randn(8, input_dim)  # batch of 8
out = model(dummy_x)
print("output shape:", out.shape)  # expect (8, 5)
assert out.shape == (8, output_dim)
print("OK")

# Perform training
# pinn, dataset, (train_ds, val_ds, test_ds), history = t.train_pinn(
#     root=root,
#     n_epochs=100,
#     batch_size=64,
#     lr=1e-3,
#     warmup_epochs=10,
#     min_lambda_phys=0,#1e4,
#     max_lambda_phys=0,#1e8,
#     lambda_data=1.0,
#     hidden_dims=(3208, 802, 128, 64, 32),
#     num_workers=4,
#     checkpoint_path="best_pinn_dataset1.pt",
#     seed=0,
# )
pinn, dataset, (train_ds, val_ds, test_ds), history = t.train_pinn(
    root=root,
    n_epochs=1000,
    batch_size=64,
    lr=1e-5,
    warmup_epochs=10,
    min_lambda_phys=1.e6,#1e4,
    max_lambda_phys=1.e6,#1e8,
    lambda_data=.0,
    hidden_dims=(3208, 802, 128, 64, 32),
    num_workers=4,
    checkpoint_path="best_pinn_dataset1.pt",
    seed=0,
)

print("Training complete. Best checkpoint loaded.")

# save the test split so it can be reused later
test_subset_path = "test_ds.pt"
torch.save(test_ds, test_subset_path)
print(f"Saved test subset to {test_subset_path}")

# save history so you can replot later without retraining
import json
with open("history_dataset1.json", "w") as f:
    json.dump(history, f)

plt.plot_history(history, save_path="./results/pinn_training_dataset3.png")

