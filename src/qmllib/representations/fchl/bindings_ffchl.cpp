#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// Fortran function declarations for scalar kernels
extern "C" {
    void fget_kernels_fchl_wrapper(
        const double* x1, const double* x2, const int* verbose,
        const int* n1, const int* n2, 
        const int* nneigh1, const int* nneigh2,
        const int* nm1, const int* nm2, const int* nsigmas,
        const double* t_width, const double* d_width,
        const double* cut_start, const double* cut_distance,
        const int* order, const double* pd,
        const double* distance_scale, const double* angular_scale,
        const int* alchemy, const double* two_body_power, const double* three_body_power,
        const int* kernel_idx, const double* parameters,
        double* kernels,
        // Dimension parameters
        const int* dim_x1_1, const int* dim_x1_2, const int* dim_x1_3, const int* dim_x1_4,
        const int* dim_x2_1, const int* dim_x2_2, const int* dim_x2_3, const int* dim_x2_4,
        const int* dim_n1, const int* dim_n2,
        const int* dim_nneigh1_1, const int* dim_nneigh1_2,
        const int* dim_nneigh2_1, const int* dim_nneigh2_2,
        const int* dim_pd_1, const int* dim_pd_2,
        const int* dim_parameters_1, const int* dim_parameters_2
    );
}

// Wrapper for fget_kernels_fchl
py::array_t<double> get_kernels_fchl_wrapper(
    py::array_t<double, py::array::f_style | py::array::forcecast> x1_in,
    py::array_t<double, py::array::f_style | py::array::forcecast> x2_in,
    bool verbose,
    py::array_t<int, py::array::f_style | py::array::forcecast> n1_in,
    py::array_t<int, py::array::f_style | py::array::forcecast> n2_in,
    py::array_t<int, py::array::f_style | py::array::forcecast> nneigh1_in,
    py::array_t<int, py::array::f_style | py::array::forcecast> nneigh2_in,
    int nm1, int nm2, int nsigmas,
    double t_width, double d_width,
    double cut_start, double cut_distance,
    int order,
    py::array_t<double, py::array::f_style | py::array::forcecast> pd_in,
    double distance_scale, double angular_scale,
    bool alchemy,
    double two_body_power, double three_body_power,
    int kernel_idx,
    py::array_t<double, py::array::f_style | py::array::forcecast> parameters_in
) {
    // Ensure converted arrays stay alive
    auto x1 = py::array_t<double, py::array::f_style | py::array::forcecast>(x1_in);
    auto x2 = py::array_t<double, py::array::f_style | py::array::forcecast>(x2_in);
    auto n1 = py::array_t<int, py::array::f_style | py::array::forcecast>(n1_in);
    auto n2 = py::array_t<int, py::array::f_style | py::array::forcecast>(n2_in);
    auto nneigh1 = py::array_t<int, py::array::f_style | py::array::forcecast>(nneigh1_in);
    auto nneigh2 = py::array_t<int, py::array::f_style | py::array::forcecast>(nneigh2_in);
    auto pd = py::array_t<double, py::array::f_style | py::array::forcecast>(pd_in);
    auto parameters = py::array_t<double, py::array::f_style | py::array::forcecast>(parameters_in);
    
    auto buf_x1 = x1.request();
    auto buf_x2 = x2.request();
    auto buf_n1 = n1.request();
    auto buf_n2 = n2.request();
    auto buf_nneigh1 = nneigh1.request();
    auto buf_nneigh2 = nneigh2.request();
    auto buf_pd = pd.request();
    auto buf_parameters = parameters.request();
    
    // Get dimensions
    int dim_x1_1 = buf_x1.shape[0];
    int dim_x1_2 = buf_x1.shape[1];
    int dim_x1_3 = buf_x1.shape[2];
    int dim_x1_4 = buf_x1.shape[3];
    
    int dim_x2_1 = buf_x2.shape[0];
    int dim_x2_2 = buf_x2.shape[1];
    int dim_x2_3 = buf_x2.shape[2];
    int dim_x2_4 = buf_x2.shape[3];
    
    int dim_n1 = buf_n1.shape[0];
    int dim_n2 = buf_n2.shape[0];
    
    int dim_nneigh1_1 = buf_nneigh1.shape[0];
    int dim_nneigh1_2 = buf_nneigh1.shape[1];
    int dim_nneigh2_1 = buf_nneigh2.shape[0];
    int dim_nneigh2_2 = buf_nneigh2.shape[1];
    
    int dim_pd_1 = buf_pd.shape[0];
    int dim_pd_2 = buf_pd.shape[1];
    
    int dim_parameters_1 = buf_parameters.shape[0];
    int dim_parameters_2 = buf_parameters.shape[1];
    
    // Create output array - Fortran column-major
    std::vector<ssize_t> shape = {nsigmas, nm1, nm2};
    std::vector<ssize_t> strides = {sizeof(double), sizeof(double) * nsigmas, sizeof(double) * nsigmas * nm1};
    auto kernels = py::array_t<double>(shape, strides);
    auto buf_kernels = kernels.request();
    
    // Convert booleans to int for Fortran
    int verbose_int = verbose ? 1 : 0;
    int alchemy_int = alchemy ? 1 : 0;
    
    fget_kernels_fchl_wrapper(
        static_cast<const double*>(buf_x1.ptr),
        static_cast<const double*>(buf_x2.ptr),
        &verbose_int,
        static_cast<const int*>(buf_n1.ptr),
        static_cast<const int*>(buf_n2.ptr),
        static_cast<const int*>(buf_nneigh1.ptr),
        static_cast<const int*>(buf_nneigh2.ptr),
        &nm1, &nm2, &nsigmas,
        &t_width, &d_width,
        &cut_start, &cut_distance,
        &order,
        static_cast<const double*>(buf_pd.ptr),
        &distance_scale, &angular_scale,
        &alchemy_int,
        &two_body_power, &three_body_power,
        &kernel_idx,
        static_cast<const double*>(buf_parameters.ptr),
        static_cast<double*>(buf_kernels.ptr),
        // Dimensions
        &dim_x1_1, &dim_x1_2, &dim_x1_3, &dim_x1_4,
        &dim_x2_1, &dim_x2_2, &dim_x2_3, &dim_x2_4,
        &dim_n1, &dim_n2,
        &dim_nneigh1_1, &dim_nneigh1_2,
        &dim_nneigh2_1, &dim_nneigh2_2,
        &dim_pd_1, &dim_pd_2,
        &dim_parameters_1, &dim_parameters_2
    );
    
    return kernels;
}

PYBIND11_MODULE(ffchl_module, m) {
    m.doc() = "QMLlib FCHL representation functions";

    // Create a submodule for kernel types (constants)
    py::module_ kt = m.def_submodule("ffchl_kernel_types", "Kernel type constants");
    kt.attr("GAUSSIAN") = 1;
    kt.attr("LINEAR") = 2;
    kt.attr("POLYNOMIAL") = 3;
    kt.attr("SIGMOID") = 4;
    kt.attr("MULTIQUADRATIC") = 5;
    kt.attr("INV_MULTIQUADRATIC") = 6;
    kt.attr("BESSEL") = 7;
    kt.attr("L2") = 8;
    kt.attr("MATERN") = 9;
    kt.attr("CAUCHY") = 10;
    kt.attr("POLYNOMIAL2") = 11;

    m.def("fget_kernels_fchl", &get_kernels_fchl_wrapper,
        py::arg("x1"), py::arg("x2"), py::arg("verbose"),
        py::arg("n1"), py::arg("n2"),
        py::arg("nneigh1"), py::arg("nneigh2"),
        py::arg("nm1"), py::arg("nm2"), py::arg("nsigmas"),
        py::arg("t_width"), py::arg("d_width"),
        py::arg("cut_start"), py::arg("cut_distance"),
        py::arg("order"), py::arg("pd"),
        py::arg("distance_scale"), py::arg("angular_scale"),
        py::arg("alchemy"),
        py::arg("two_body_power"), py::arg("three_body_power"),
        py::arg("kernel_idx"), py::arg("parameters"),
        "FCHL kernel computation");
}
