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

n_elements = len(REP_PARAMS["elements"])
nRs2 = REP_PARAMS["nRs2"]
nRs3 = REP_PARAMS["nRs3"]
nFourier = REP_PARAMS["nFourier"]

# Calculate sizes
twobody_size = n_elements * nRs2  # 4 * 24 = 96
threebody_size_per_pair = nRs3 * 2 * nFourier  # 20 * 2 * 1 = 40
n_element_pairs = n_elements * (n_elements + 1)  # 4 * 5 = 20
threebody_size = n_element_pairs * threebody_size_per_pair  # 20 * 40 = 800
total_size = twobody_size + threebody_size  # 96 + 800 = 896

print(f"Two-body size: {twobody_size}")
print(f"Three-body size: {threebody_size}")
print(f"Total expected size: {total_size}")

print("\nGenerating with gradients...")
(repa, anal_grad) = generate_fchl19(
    nuclear_charges, coordinates, gradients=True, **REP_PARAMS
)

print("Generating without gradients...")
repb = generate_fchl19(nuclear_charges, coordinates, gradients=False, **REP_PARAMS)

print(f"\nActual shape: {repa.shape}")

# Check two-body terms separately
twobody_a = repa[:, :twobody_size]
twobody_b = repb[:, :twobody_size]

twobody_diff = np.abs(twobody_a - twobody_b)
print(f"\nTwo-body max diff: {np.max(twobody_diff)}")
print(f"Two-body differences > 1e-6: {np.sum(twobody_diff > 1e-6)}")

# Check three-body terms separately
threebody_a = repa[:, twobody_size:]
threebody_b = repb[:, twobody_size:]

threebody_diff = np.abs(threebody_a - threebody_b)
print(f"\nThree-body max diff: {np.max(threebody_diff)}")
print(f"Three-body differences > 1e-6: {np.sum(threebody_diff > 1e-6)}")

# Show first few two-body term comparisons
print("\nFirst atom, first few two-body terms:")
for i in range(min(10, twobody_size)):
    if twobody_diff[0, i] > 1e-10:
        print(
            f"  [{i}]: a={twobody_a[0, i]:.10e}, b={twobody_b[0, i]:.10e}, diff={twobody_diff[0, i]:.10e}"
        )
