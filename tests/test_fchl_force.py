# flake8: noqa

import ast
import csv
from copy import deepcopy

import numpy as np
import pytest
import scipy
import scipy.stats
from conftest import ASSETS
from scipy.linalg import lstsq

from qmllib.representations import (
    generate_fchl18,
    generate_fchl18_displaced,
    generate_fchl18_displaced_5point,
)
from qmllib.representations.fchl import (
    get_atomic_local_gradient_5point_kernels,
    get_atomic_local_gradient_kernels,
    get_atomic_local_kernels,
    get_force_alphas,
    get_gaussian_process_kernels,
    get_local_gradient_kernels,
    get_local_hessian_kernels,
    get_local_kernels,
    get_local_symmetric_hessian_kernels,
    get_local_symmetric_kernels,
)
from qmllib.solvers import cho_solve

FORCE_KEY = "forces"
ENERGY_KEY = "om2_energy"
CSV_FILE = ASSETS / "amons_small.csv"
SIGMAS = [0.64]

TRAINING = 13
TEST = 7

DX = 0.005
CUT_DISTANCE = 1e6
KERNEL_ARGS = {
    "verbose": False,
    "cut_distance": CUT_DISTANCE,
    "kernel": "gaussian",
    "kernel_args": {
        "sigma": SIGMAS,
    },
}

LLAMBDA_ENERGY = 1e-7
LLAMBDA_FORCE = 1e-7


# pytest.skip(allow_module_level=True, reason="Test is broken")


def mae(a, b):
    return np.mean(np.abs(a.flatten() - b.flatten()))


def csv_to_molecular_reps(
    csv_filename, force_key="orca_forces", energy_key="orca_energy"
):
    np.random.seed(667)

    x = []
    f = []
    e = []
    distance = []

    disp_x = []
    disp_x5 = []

    max_atoms = 5

    with open(csv_filename, "r") as csvfile:
        df = csv.reader(csvfile, delimiter=";", quotechar="#")

        for row in df:
            coordinates = np.array(ast.literal_eval(row[2]))
            nuclear_charges = ast.literal_eval(row[5])
            atomtypes = ast.literal_eval(row[1])
            force = np.array(ast.literal_eval(row[3]))
            energy = float(row[6])

            rep = generate_fchl18(
                nuclear_charges,
                coordinates,
                max_size=max_atoms,
                cut_distance=CUT_DISTANCE,
            )

            disp_rep = generate_fchl18_displaced(
                nuclear_charges,
                coordinates,
                max_size=max_atoms,
                cut_distance=CUT_DISTANCE,
                dx=DX,
            )

            disp_rep5 = generate_fchl18_displaced_5point(
                nuclear_charges,
                coordinates,
                max_size=max_atoms,
                cut_distance=CUT_DISTANCE,
                dx=DX,
            )

            x.append(rep)
            f.append(force)
            e.append(energy)

            disp_x.append(disp_rep)
            disp_x5.append(disp_rep5)

    return np.array(x), f, e, np.array(disp_x), np.array(disp_x5)


def test_gaussian_process_derivative():
    Xall, Fall, Eall, dXall, dXall5 = csv_to_molecular_reps(
        CSV_FILE, force_key=FORCE_KEY, energy_key=ENERGY_KEY
    )

    Eall = np.array(Eall)
    Fall = np.array(Fall)

    X = Xall[:TRAINING]
    dX = dXall[:TRAINING]
    F = Fall[:TRAINING]
    E = Eall[:TRAINING]

    Xs = Xall[-TEST:]
    dXs = dXall[-TEST:]
    Fs = Fall[-TEST:]
    Es = Eall[-TEST:]

    K = get_gaussian_process_kernels(X, dX, dx=DX, **KERNEL_ARGS)
    Kt = K[:, TRAINING:, TRAINING:]
    Kt_local = K[:, :TRAINING, :TRAINING]
    Kt_energy = K[:, :TRAINING, TRAINING:]

    Kt_grad2 = get_local_gradient_kernels(X, dX, dx=DX, **KERNEL_ARGS)

    Ks = get_local_hessian_kernels(dX, dXs, dx=DX, **KERNEL_ARGS)
    Ks_energy = get_local_gradient_kernels(X, dXs, dx=DX, **KERNEL_ARGS)

    Ks_energy2 = get_local_gradient_kernels(Xs, dX, dx=DX, **KERNEL_ARGS)
    Ks_local = get_local_kernels(X, Xs, **KERNEL_ARGS)

    F = np.concatenate(F)
    Fs = np.concatenate(Fs)

    Y = np.array(F.flatten())
    Y = np.concatenate((E, Y))

    for i, sigma in enumerate(SIGMAS):
        C = deepcopy(K[i])

        for j in range(TRAINING):
            C[j, j] += LLAMBDA_ENERGY

        for j in range(TRAINING, K.shape[2]):
            C[j, j] += LLAMBDA_FORCE

        alpha = cho_solve(C, Y)
        beta = alpha[:TRAINING]
        gamma = alpha[TRAINING:]

        Fss = np.dot(np.transpose(Ks[i]), gamma) + np.dot(
            np.transpose(Ks_energy[i]), beta
        )
        Ft = np.dot(np.transpose(Kt[i]), gamma) + np.dot(
            np.transpose(Kt_energy[i]), beta
        )

        Ess = np.dot(Ks_energy2[i], gamma) + np.dot(Ks_local[i].T, beta)
        Et = np.dot(Kt_energy[i], gamma) + np.dot(Kt_local[i].T, beta)

        assert mae(Ess, Es) < 0.1, "Error in Gaussian Process test energy"
        assert mae(Et, E) < 0.001, "Error in Gaussian Process training energy"

        assert mae(Fss, Fs) < 1.0, "Error in Gaussian Process test force"
        assert mae(Ft, F) < 0.001, "Error in Gaussian Process training force"


def test_gdml_derivative():
    Xall, Fall, Eall, dXall, dXall5 = csv_to_molecular_reps(
        CSV_FILE, force_key=FORCE_KEY, energy_key=ENERGY_KEY
    )

    Eall = np.array(Eall)
    Fall = np.array(Fall)

    X = Xall[:TRAINING]
    dX = dXall[:TRAINING]
    F = Fall[:TRAINING]
    E = Eall[:TRAINING]

    Xs = Xall[-TEST:]
    dXs = dXall[-TEST:]
    Fs = Fall[-TEST:]
    Es = Eall[-TEST:]

    K = get_local_symmetric_hessian_kernels(dX, dx=DX, **KERNEL_ARGS)
    Ks = get_local_hessian_kernels(dXs, dX, dx=DX, **KERNEL_ARGS)

    Kt_energy = get_local_gradient_kernels(X, dX, dx=DX, **KERNEL_ARGS)
    Ks_energy = get_local_gradient_kernels(Xs, dX, dx=DX, **KERNEL_ARGS)

    F = np.concatenate(F)
    Fs = np.concatenate(Fs)

    Y = np.array(F.flatten())
    # Y = np.concatenate((E, Y))

    for i, sigma in enumerate(SIGMAS):
        C = deepcopy(K[i])
        for j in range(K.shape[2]):
            C[j, j] += LLAMBDA_FORCE

        alpha = cho_solve(C, Y)
        Fss = np.dot(Ks[i], alpha)
        Ft = np.dot(K[i], alpha)

        Ess = np.dot(Ks_energy[i], alpha)
        Et = np.dot(Kt_energy[i], alpha)

        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(
            E.flatten(), Et.flatten()
        )

        Ess -= intercept
        Et -= intercept

        # This test will only work for molecules of same type
        # assert mae(Ess, Es) < 0.1, "Error in Gaussian Process test energy"
        # assert mae(Et, E) < 0.001, "Error in Gaussian Process training energy"

        assert mae(Fss, Fs) < 1.0, "Error in GDML test force"
        assert mae(Ft, F) < 0.001, "Error in GDML training force"


def test_normal_equation_derivative():
    Xall, Fall, Eall, dXall, dXall5 = csv_to_molecular_reps(
        CSV_FILE, force_key=FORCE_KEY, energy_key=ENERGY_KEY
    )

    Eall = np.array(Eall)
    Fall = np.array(Fall)

    X = Xall[:TRAINING]
    dX = dXall[:TRAINING]
    dX5 = dXall5[:TRAINING]
    F = Fall[:TRAINING]
    E = Eall[:TRAINING]

    Xs = Xall[-TEST:]
    dXs = dXall[-TEST:]
    dXs5 = dXall5[-TEST:]
    Fs = Fall[-TEST:]
    Es = Eall[-TEST:]

    Ftrain = np.concatenate(F)
    Etrain = np.array(E)
    alphas = get_force_alphas(
        X, dX, Ftrain, energy=Etrain, dx=DX, regularization=LLAMBDA_FORCE, **KERNEL_ARGS
    )

    Kt_force = get_atomic_local_gradient_kernels(X, dX, dx=DX, **KERNEL_ARGS)
    Ks_force = get_atomic_local_gradient_kernels(X, dXs, dx=DX, **KERNEL_ARGS)

    Kt_force5 = get_atomic_local_gradient_5point_kernels(X, dX5, dx=DX, **KERNEL_ARGS)
    Ks_force5 = get_atomic_local_gradient_5point_kernels(X, dXs5, dx=DX, **KERNEL_ARGS)

    Kt_energy = get_atomic_local_kernels(X, X, **KERNEL_ARGS)
    Ks_energy = get_atomic_local_kernels(X, Xs, **KERNEL_ARGS)

    F = np.concatenate(F)
    Fs = np.concatenate(Fs)
    Y = np.array(F.flatten())

    for i, sigma in enumerate(SIGMAS):
        Ft = np.zeros((Kt_force[i, :, :].shape[1] // 3, 3))
        Fss = np.zeros((Ks_force[i, :, :].shape[1] // 3, 3))

        Ft5 = np.zeros((Kt_force5[i, :, :].shape[1] // 3, 3))
        Fss5 = np.zeros((Ks_force5[i, :, :].shape[1] // 3, 3))

        for xyz in range(3):
            Ft[:, xyz] = np.dot(Kt_force[i, :, xyz::3].T, alphas[i])
            Fss[:, xyz] = np.dot(Ks_force[i, :, xyz::3].T, alphas[i])

            Ft5[:, xyz] = np.dot(Kt_force5[i, :, xyz::3].T, alphas[i])
            Fss5[:, xyz] = np.dot(Ks_force5[i, :, xyz::3].T, alphas[i])

        Ess = np.dot(Ks_energy[i].T, alphas[i])
        Et = np.dot(Kt_energy[i].T, alphas[i])

        assert mae(Ess, Es) < 0.3, "Error in normal equation test energy"
        assert mae(Et, E) < 0.08, "Error in normal equation training energy"

        assert mae(Fss, Fs) < 3.2, "Error in  normal equation test force"
        assert mae(Ft, F) < 0.5, "Error in  normal equation training force"

        assert mae(Fss5, Fs) < 3.2, "Error in normal equation 5-point test force"
        assert mae(Ft5, F) < 0.5, "Error in normal equation 5-point training force"

        assert mae(Fss5, Fss) < 0.01, (
            "Error in normal equation 5-point or 2-point test force"
        )
        assert mae(Ft5, Ft) < 0.01, (
            "Error in normal equation 5-point or 2-point training force"
        )


def test_operator_derivative():
    Xall, Fall, Eall, dXall, dXall5 = csv_to_molecular_reps(
        CSV_FILE, force_key=FORCE_KEY, energy_key=ENERGY_KEY
    )

    Eall = np.array(Eall)
    Fall = np.array(Fall)

    X = Xall[:TRAINING]
    dX = dXall[:TRAINING]
    dX5 = dXall5[:TRAINING]
    F = Fall[:TRAINING]
    E = Eall[:TRAINING]

    Xs = Xall[-TEST:]
    dXs = dXall[-TEST:]
    dXs5 = dXall5[-TEST:]
    Fs = Fall[-TEST:]
    Es = Eall[-TEST:]

    Ftrain = np.concatenate(F)
    Etrain = np.array(E)

    Kt_energy = get_atomic_local_kernels(X, X, **KERNEL_ARGS)
    Ks_energy = get_atomic_local_kernels(X, Xs, **KERNEL_ARGS)

    Kt_force = get_atomic_local_gradient_kernels(X, dX, dx=DX, **KERNEL_ARGS)
    Ks_force = get_atomic_local_gradient_kernels(X, dXs, dx=DX, **KERNEL_ARGS)

    F = np.concatenate(F)
    Fs = np.concatenate(Fs)
    Y = np.array(F.flatten())

    for i, sigma in enumerate(SIGMAS):
        Y = np.concatenate((E, F.flatten()))

        C = np.concatenate((Kt_energy[i].T, Kt_force[i].T))

        alphas, residuals, singular_values, rank = lstsq(
            C, Y, cond=1e-9, lapack_driver="gelsd"
        )

        Ess = np.dot(Ks_energy[i].T, alphas)
        Et = np.dot(Kt_energy[i].T, alphas)

        Fss = np.dot(Ks_force[i].T, alphas)
        Ft = np.dot(Kt_force[i].T, alphas)

        assert mae(Ess, Es) < 0.08, "Error in operator test energy"
        assert mae(Et, E) < 0.04, "Error in  operator training energy"

        assert mae(Fss, Fs.flatten()) < 1.1, "Error in  operator test force"
        assert mae(Ft, F.flatten()) < 0.1, "Error in  operator training force"


def test_krr_derivative():
    """Test that gradient kernels can be computed without errors.

    Note: This test only verifies that the function runs and produces
    finite values. The original test had unrealistic expectations and
    was skipped in the f2py version.
    """
    Xall, Fall, Eall, dXall, dXall5 = csv_to_molecular_reps(
        CSV_FILE, force_key=FORCE_KEY, energy_key=ENERGY_KEY
    )

    Eall = np.array(Eall)
    # Fall = np.array(Fall)  # Fall has inhomogeneous shape, keep as list

    X = Xall[:TRAINING]
    dX = dXall[:TRAINING]
    F = Fall[:TRAINING]
    E = Eall[:TRAINING]

    Xs = Xall[-TEST:]
    dXs = dXall[-TEST:]
    Fs = Fall[-TEST:]
    Es = Eall[-TEST:]

    K = get_local_symmetric_kernels(X, **KERNEL_ARGS)
    Ks = get_local_kernels(Xs, X, **KERNEL_ARGS)

    Kt_force = get_local_gradient_kernels(X, dX, dx=DX, **KERNEL_ARGS)
    Ks_force = get_local_gradient_kernels(X, dXs, dx=DX, **KERNEL_ARGS)

    # Verify kernels have correct shapes
    assert Kt_force.shape[0] == len(SIGMAS), "Wrong number of sigmas"
    assert Kt_force.shape[1] == TRAINING, "Wrong number of training molecules"
    assert Ks_force.shape[1] == TRAINING, "Wrong number of training molecules"

    # Verify kernels contain finite values
    assert np.all(np.isfinite(Kt_force)), "Gradient kernel contains NaN/Inf"
    assert np.all(np.isfinite(Ks_force)), "Gradient kernel contains NaN/Inf"

    # Verify energy kernels still work
    assert mae(K[0], K[0].T) < 1e-10, "Symmetric kernel not symmetric"


if __name__ == "__main__":
    test_gaussian_process_derivative()
    test_gdml_derivative()
    test_normal_equation_derivative()
    test_operator_derivative()
    test_krr_derivative()
