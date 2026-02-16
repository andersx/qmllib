import numpy as np
from pathlib import Path
from qmllib.utils.xyz_format import read_xyz
from qmllib.representations import generate_fchl19

ASSETS = Path("tests/assets")
REP_PARAMS = {
    "elements": [1, 6, 7, 8],
    "nRs2": 24,
    "nRs3": 20,
    "nFourier": 1,
    "eta2": 0.32,
    "eta3": 2.7,
    "zeta": np.pi,
    "rcut": 8.0,
    "acut": 8.0,
    "two_body_decay": 1.8,
    "three_body_decay": 0.57,
    "three_body_weight": 13.4,
}

coordinates, nuclear_charges = read_xyz(ASSETS / "qm7/0101.xyz")

print("Generating with gradients...")
(repa, anal_grad) = generate_fchl19(
    nuclear_charges, coordinates, gradients=True, **REP_PARAMS
)

print("Generating without gradients...")
repb = generate_fchl19(nuclear_charges, coordinates, gradients=False, **REP_PARAMS)

print(f"\nrepa shape: {repa.shape}")
print(f"repb shape: {repb.shape}")

# Find differences
diff = np.abs(repa - repb)
max_diff = np.max(diff)
rel_diff = diff / (np.abs(repb) + 1e-10)
max_rel_diff = np.max(rel_diff)

print(f"\nMax absolute difference: {max_diff}")
print(f"Max relative difference: {max_rel_diff}")

# Find where differences are largest
large_diff_mask = diff > 1e-6
if np.any(large_diff_mask):
    print(f"\nNumber of large differences (>1e-6): {np.sum(large_diff_mask)}")
    indices = np.where(large_diff_mask)
    for i in range(min(10, len(indices[0]))):
        row, col = indices[0][i], indices[1][i]
        print(
            f"  [{row}, {col}]: repa={repa[row, col]:.10e}, repb={repb[row, col]:.10e}, diff={diff[row, col]:.10e}"
        )

# Check column-wise differences
print("\nColumn-wise max differences:")
col_diffs = np.max(diff, axis=0)
problem_cols = np.where(col_diffs > 1e-6)[0]
print(f"Columns with large differences: {problem_cols[:20]}")

# Print some values to understand the pattern
print("\nFirst row comparison (first 10 elements):")
print("repa:", repa[0, :10])
print("repb:", repb[0, :10])
print("diff:", diff[0, :10])
