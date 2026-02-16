import numpy as np
from pathlib import Path
from qmllib.utils.xyz_format import read_xyz

ASSETS = Path("tests/assets")
coordinates, nuclear_charges = read_xyz(ASSETS / "qm7/0101.xyz")

print(f"Coordinates shape: {coordinates.shape}")
print(f"Coordinates flags:\n{coordinates.flags}")
print(f"First atom coordinates: {coordinates[0]}")
print(f"\nNuclear charges shape: {nuclear_charges.shape}")
print(f"Nuclear charges: {nuclear_charges}")
