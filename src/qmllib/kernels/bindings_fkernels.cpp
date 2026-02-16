#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>

namespace py = pybind11;

// Declare C ABI Fortran functions
extern "C" {
    void fkpca(const double* k, int n, int centering, double* kpca);
    void fwasserstein_kernel(const double* a, int rep_size, int na,
                            const double* b, int nb,
                            double* k, double sigma, int p, int q);
}

// Wrapper for fkpca
py::array_t<double> kpca_wrapper(
    py::array_t<double, py::array::f_style | py::array::forcecast> k,
    int n,
    bool centering
) {
    auto bufK = k.request();
    
    if (bufK.ndim != 2) {
        throw std::runtime_error("K must be a 2D array");
    }
    
    int size = static_cast<int>(bufK.shape[0]);
    
    if (bufK.shape[0] != bufK.shape[1]) {
        throw std::runtime_error("K must be a square matrix");
    }
    
    if (size != n) {
        throw std::runtime_error("K dimensions must match n parameter");
    }
    
    // Create Fortran-style (column-major) output array
    std::vector<ssize_t> shape = {n, n};
    std::vector<ssize_t> strides = {sizeof(double), sizeof(double) * n};
    auto kpca = py::array_t<double>(shape, strides);
    auto bufKPCA = kpca.request();
    
    // Call Fortran function (0=false, 1=true for centering)
    fkpca(
        static_cast<const double*>(bufK.ptr),
        n,
        centering ? 1 : 0,
        static_cast<double*>(bufKPCA.ptr)
    );
    
    return kpca;
}

// Wrapper for fwasserstein_kernel
py::array_t<double> wasserstein_kernel_wrapper(
    py::array_t<double, py::array::f_style | py::array::forcecast> a,
    int na,
    py::array_t<double, py::array::f_style | py::array::forcecast> b,
    int nb,
    double sigma,
    int p,
    int q
) {
    auto bufA = a.request();
    auto bufB = b.request();
    
    if (bufA.ndim != 2 || bufB.ndim != 2) {
        throw std::runtime_error("A and B must be 2D arrays");
    }
    
    int rep_size = static_cast<int>(bufA.shape[0]);
    
    if (bufA.shape[0] != bufB.shape[0]) {
        throw std::runtime_error("A and B must have same representation size");
    }
    
    if (bufA.shape[1] != na) {
        throw std::runtime_error("A second dimension must match na");
    }
    
    if (bufB.shape[1] != nb) {
        throw std::runtime_error("B second dimension must match nb");
    }
    
    // Create Fortran-style (column-major) output array
    std::vector<ssize_t> shape = {na, nb};
    std::vector<ssize_t> strides = {sizeof(double), sizeof(double) * na};
    auto k = py::array_t<double>(shape, strides);
    auto bufK = k.request();
    
    // Initialize to zero
    std::memset(bufK.ptr, 0, na * nb * sizeof(double));
    
    fwasserstein_kernel(
        static_cast<const double*>(bufA.ptr),
        rep_size, na,
        static_cast<const double*>(bufB.ptr),
        nb,
        static_cast<double*>(bufK.ptr),
        sigma, p, q
    );
    
    return k;
}

PYBIND11_MODULE(_fkernels, m) {
    m.doc() = "QMLlib kernel functions (KPCA and Wasserstein)";

    m.def("fkpca", &kpca_wrapper,
        py::arg("k"), py::arg("n"), py::arg("centering"),
        "Kernel PCA decomposition");

    m.def("fwasserstein_kernel", &wasserstein_kernel_wrapper,
        py::arg("a"), py::arg("na"), py::arg("b"), py::arg("nb"),
        py::arg("sigma"), py::arg("p"), py::arg("q"),
        "Wasserstein kernel computation");
}
