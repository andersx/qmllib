! Wrapper subroutines with bind(C) for calling from pybind11
! These wrappers convert explicit-size arrays to assumed-shape arrays
! and call the original FCHL scalar kernel subroutines

! Wrapper 1: fget_kernels_fchl
subroutine fget_kernels_fchl_wrapper(x1, x2, verbose, n1, n2, nneigh1, nneigh2, nm1, nm2, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels, &
       & dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4, &
       & dim_x2_1, dim_x2_2, dim_x2_3, dim_x2_4, &
       & dim_n1, dim_n2, dim_nneigh1_1, dim_nneigh1_2, dim_nneigh2_1, dim_nneigh2_2, &
       & dim_pd_1, dim_pd_2, dim_parameters_1, dim_parameters_2) bind(C, name="fget_kernels_fchl_wrapper")

   use, intrinsic :: iso_c_binding
   implicit none

   integer(c_int), intent(in), value :: dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4
   integer(c_int), intent(in), value :: dim_x2_1, dim_x2_2, dim_x2_3, dim_x2_4
   integer(c_int), intent(in), value :: dim_n1, dim_n2, dim_nneigh1_1, dim_nneigh1_2
   integer(c_int), intent(in), value :: dim_nneigh2_1, dim_nneigh2_2, dim_pd_1, dim_pd_2
   integer(c_int), intent(in), value :: dim_parameters_1, dim_parameters_2

   double precision, dimension(dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4), intent(in) :: x1
   double precision, dimension(dim_x2_1, dim_x2_2, dim_x2_3, dim_x2_4), intent(in) :: x2
   integer(c_int), dimension(dim_n1), intent(in) :: n1
   integer(c_int), dimension(dim_n2), intent(in) :: n2
   integer(c_int), dimension(dim_nneigh1_1, dim_nneigh1_2), intent(in) :: nneigh1
   integer(c_int), dimension(dim_nneigh2_1, dim_nneigh2_2), intent(in) :: nneigh2
   double precision, dimension(dim_pd_1, dim_pd_2), intent(in) :: pd
   double precision, dimension(dim_parameters_1, dim_parameters_2), intent(in) :: parameters

   integer(c_int), intent(in), value :: verbose, nm1, nm2, nsigmas, order, kernel_idx, alchemy
   double precision, intent(in), value :: t_width, d_width, cut_start, cut_distance
   double precision, intent(in), value :: distance_scale, angular_scale, two_body_power, three_body_power

   double precision, dimension(nsigmas, nm1, nm2), intent(out) :: kernels

   logical :: verbose_log, alchemy_log

   verbose_log = (verbose /= 0)
   alchemy_log = (alchemy /= 0)

   call fget_kernels_fchl(x1, x2, verbose_log, n1, n2, nneigh1, nneigh2, nm1, nm2, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy_log, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels)

end subroutine fget_kernels_fchl_wrapper

! Wrapper 2: fget_symmetric_kernels_fchl
subroutine fget_symmetric_kernels_fchl_wrapper(x1, verbose, n1, nneigh1, nm1, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels, &
       & dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4, dim_n1, dim_nneigh1_1, dim_nneigh1_2, &
       & dim_pd_1, dim_pd_2, dim_parameters_1, dim_parameters_2) bind(C, name="fget_symmetric_kernels_fchl_wrapper")

   use, intrinsic :: iso_c_binding
   implicit none

   integer(c_int), intent(in), value :: dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4, dim_n1
   integer(c_int), intent(in), value :: dim_nneigh1_1, dim_nneigh1_2, dim_pd_1, dim_pd_2
   integer(c_int), intent(in), value :: dim_parameters_1, dim_parameters_2

   double precision, dimension(dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4), intent(in) :: x1
   integer(c_int), dimension(dim_n1), intent(in) :: n1
   integer(c_int), dimension(dim_nneigh1_1, dim_nneigh1_2), intent(in) :: nneigh1
   double precision, dimension(dim_pd_1, dim_pd_2), intent(in) :: pd
   double precision, dimension(dim_parameters_1, dim_parameters_2), intent(in) :: parameters

   integer(c_int), intent(in), value :: verbose, nm1, nsigmas, order, kernel_idx, alchemy
   double precision, intent(in), value :: t_width, d_width, cut_start, cut_distance
   double precision, intent(in), value :: distance_scale, angular_scale, two_body_power, three_body_power

   double precision, dimension(nsigmas, nm1, nm1), intent(out) :: kernels

   logical :: verbose_log, alchemy_log

   verbose_log = (verbose /= 0)
   alchemy_log = (alchemy /= 0)

   call fget_symmetric_kernels_fchl(x1, verbose_log, n1, nneigh1, nm1, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy_log, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels)

end subroutine fget_symmetric_kernels_fchl_wrapper

! Wrapper 3: fget_global_symmetric_kernels_fchl  
subroutine fget_global_symmetric_kernels_fchl_wrapper(x1, verbose, n1, nneigh1, nm1, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels, &
       & dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4, dim_n1, dim_nneigh1_1, dim_nneigh1_2, &
       & dim_pd_1, dim_pd_2, dim_parameters_1, dim_parameters_2) bind(C, name="fget_global_symmetric_kernels_fchl_wrapper")

   use, intrinsic :: iso_c_binding
   implicit none

   integer(c_int), intent(in), value :: dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4, dim_n1
   integer(c_int), intent(in), value :: dim_nneigh1_1, dim_nneigh1_2, dim_pd_1, dim_pd_2
   integer(c_int), intent(in), value :: dim_parameters_1, dim_parameters_2

   double precision, dimension(dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4), intent(in) :: x1
   integer(c_int), dimension(dim_n1), intent(in) :: n1
   integer(c_int), dimension(dim_nneigh1_1, dim_nneigh1_2), intent(in) :: nneigh1
   double precision, dimension(dim_pd_1, dim_pd_2), intent(in) :: pd
   double precision, dimension(dim_parameters_1, dim_parameters_2), intent(in) :: parameters

   integer(c_int), intent(in), value :: verbose, nm1, nsigmas, order, kernel_idx, alchemy
   double precision, intent(in), value :: t_width, d_width, cut_start, cut_distance
   double precision, intent(in), value :: distance_scale, angular_scale, two_body_power, three_body_power

   double precision, dimension(nsigmas, nm1, nm1), intent(out) :: kernels

   logical :: verbose_log, alchemy_log

   verbose_log = (verbose /= 0)
   alchemy_log = (alchemy /= 0)

   call fget_global_symmetric_kernels_fchl(x1, verbose_log, n1, nneigh1, nm1, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy_log, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels)

end subroutine fget_global_symmetric_kernels_fchl_wrapper

! Wrapper 4: fget_global_kernels_fchl
subroutine fget_global_kernels_fchl_wrapper(x1, x2, verbose, n1, n2, nneigh1, nneigh2, nm1, nm2, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels, &
       & dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4, dim_x2_1, dim_x2_2, dim_x2_3, dim_x2_4, &
       & dim_n1, dim_n2, dim_nneigh1_1, dim_nneigh1_2, dim_nneigh2_1, dim_nneigh2_2, &
       & dim_pd_1, dim_pd_2, dim_parameters_1, dim_parameters_2) bind(C, name="fget_global_kernels_fchl_wrapper")

   use, intrinsic :: iso_c_binding
   implicit none

   integer(c_int), intent(in), value :: dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4
   integer(c_int), intent(in), value :: dim_x2_1, dim_x2_2, dim_x2_3, dim_x2_4
   integer(c_int), intent(in), value :: dim_n1, dim_n2, dim_nneigh1_1, dim_nneigh1_2
   integer(c_int), intent(in), value :: dim_nneigh2_1, dim_nneigh2_2, dim_pd_1, dim_pd_2
   integer(c_int), intent(in), value :: dim_parameters_1, dim_parameters_2

   double precision, dimension(dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4), intent(in) :: x1
   double precision, dimension(dim_x2_1, dim_x2_2, dim_x2_3, dim_x2_4), intent(in) :: x2
   integer(c_int), dimension(dim_n1), intent(in) :: n1
   integer(c_int), dimension(dim_n2), intent(in) :: n2
   integer(c_int), dimension(dim_nneigh1_1, dim_nneigh1_2), intent(in) :: nneigh1
   integer(c_int), dimension(dim_nneigh2_1, dim_nneigh2_2), intent(in) :: nneigh2
   double precision, dimension(dim_pd_1, dim_pd_2), intent(in) :: pd
   double precision, dimension(dim_parameters_1, dim_parameters_2), intent(in) :: parameters

   integer(c_int), intent(in), value :: verbose, nm1, nm2, nsigmas, order, kernel_idx, alchemy
   double precision, intent(in), value :: t_width, d_width, cut_start, cut_distance
   double precision, intent(in), value :: distance_scale, angular_scale, two_body_power, three_body_power

   double precision, dimension(nsigmas, nm1, nm2), intent(out) :: kernels

   logical :: verbose_log, alchemy_log

   verbose_log = (verbose /= 0)
   alchemy_log = (alchemy /= 0)

   call fget_global_kernels_fchl(x1, x2, verbose_log, n1, n2, nneigh1, nneigh2, nm1, nm2, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy_log, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels)

end subroutine fget_global_kernels_fchl_wrapper

! Wrapper 5: fget_atomic_kernels_fchl (uses 3D arrays, not 4D)
subroutine fget_atomic_kernels_fchl_wrapper(x1, x2, verbose, nneigh1, nneigh2, na1, na2, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels, &
       & dim_x1_1, dim_x1_2, dim_x1_3, dim_x2_1, dim_x2_2, dim_x2_3, &
       & dim_nneigh1, dim_nneigh2, dim_pd_1, dim_pd_2, &
       & dim_parameters_1, dim_parameters_2) bind(C, name="fget_atomic_kernels_fchl_wrapper")

   use, intrinsic :: iso_c_binding
   implicit none

   integer(c_int), intent(in), value :: dim_x1_1, dim_x1_2, dim_x1_3
   integer(c_int), intent(in), value :: dim_x2_1, dim_x2_2, dim_x2_3
   integer(c_int), intent(in), value :: dim_nneigh1, dim_nneigh2, dim_pd_1, dim_pd_2
   integer(c_int), intent(in), value :: dim_parameters_1, dim_parameters_2

   double precision, dimension(dim_x1_1, dim_x1_2, dim_x1_3), intent(in) :: x1
   double precision, dimension(dim_x2_1, dim_x2_2, dim_x2_3), intent(in) :: x2
   integer(c_int), dimension(dim_nneigh1), intent(in) :: nneigh1
   integer(c_int), dimension(dim_nneigh2), intent(in) :: nneigh2
   double precision, dimension(dim_pd_1, dim_pd_2), intent(in) :: pd
   double precision, dimension(dim_parameters_1, dim_parameters_2), intent(in) :: parameters

   integer(c_int), intent(in), value :: verbose, na1, na2, nsigmas, order, kernel_idx, alchemy
   double precision, intent(in), value :: t_width, d_width, cut_start, cut_distance
   double precision, intent(in), value :: distance_scale, angular_scale, two_body_power, three_body_power

   double precision, dimension(nsigmas, na1, na2), intent(out) :: kernels

   logical :: verbose_log, alchemy_log

   verbose_log = (verbose /= 0)
   alchemy_log = (alchemy /= 0)

   call fget_atomic_kernels_fchl(x1, x2, verbose_log, nneigh1, nneigh2, na1, na2, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy_log, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels)

end subroutine fget_atomic_kernels_fchl_wrapper

! Wrapper 6: fget_atomic_symmetric_kernels_fchl (uses 3D arrays)
subroutine fget_atomic_symmetric_kernels_fchl_wrapper(x1, verbose, nneigh1, na1, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels, &
       & dim_x1_1, dim_x1_2, dim_x1_3, dim_nneigh1, dim_pd_1, dim_pd_2, &
       & dim_parameters_1, dim_parameters_2) bind(C, name="fget_atomic_symmetric_kernels_fchl_wrapper")

   use, intrinsic :: iso_c_binding
   implicit none

   integer(c_int), intent(in), value :: dim_x1_1, dim_x1_2, dim_x1_3, dim_nneigh1
   integer(c_int), intent(in), value :: dim_pd_1, dim_pd_2, dim_parameters_1, dim_parameters_2

   double precision, dimension(dim_x1_1, dim_x1_2, dim_x1_3), intent(in) :: x1
   integer(c_int), dimension(dim_nneigh1), intent(in) :: nneigh1
   double precision, dimension(dim_pd_1, dim_pd_2), intent(in) :: pd
   double precision, dimension(dim_parameters_1, dim_parameters_2), intent(in) :: parameters

   integer(c_int), intent(in), value :: verbose, na1, nsigmas, order, kernel_idx, alchemy
   double precision, intent(in), value :: t_width, d_width, cut_start, cut_distance
   double precision, intent(in), value :: distance_scale, angular_scale, two_body_power, three_body_power

   double precision, dimension(nsigmas, na1, na1), intent(out) :: kernels

   logical :: verbose_log, alchemy_log

   verbose_log = (verbose /= 0)
   alchemy_log = (alchemy /= 0)

   call fget_atomic_symmetric_kernels_fchl(x1, verbose_log, nneigh1, na1, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy_log, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels)

end subroutine fget_atomic_symmetric_kernels_fchl_wrapper

! Wrapper 7: fget_atomic_local_kernels_fchl (uses 4D arrays)
subroutine fget_atomic_local_kernels_fchl_wrapper(x1, x2, verbose, n1, n2, nneigh1, nneigh2, &
       & nm1, nm2, na1, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels, &
       & dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4, dim_x2_1, dim_x2_2, dim_x2_3, dim_x2_4, &
       & dim_n1, dim_n2, dim_nneigh1_1, dim_nneigh1_2, dim_nneigh2_1, dim_nneigh2_2, &
       & dim_pd_1, dim_pd_2, dim_parameters_1, dim_parameters_2) bind(C, name="fget_atomic_local_kernels_fchl_wrapper")

   use, intrinsic :: iso_c_binding
   implicit none

   integer(c_int), intent(in), value :: dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4
   integer(c_int), intent(in), value :: dim_x2_1, dim_x2_2, dim_x2_3, dim_x2_4
   integer(c_int), intent(in), value :: dim_n1, dim_n2, dim_nneigh1_1, dim_nneigh1_2
   integer(c_int), intent(in), value :: dim_nneigh2_1, dim_nneigh2_2, dim_pd_1, dim_pd_2
   integer(c_int), intent(in), value :: dim_parameters_1, dim_parameters_2

   double precision, dimension(dim_x1_1, dim_x1_2, dim_x1_3, dim_x1_4), intent(in) :: x1
   double precision, dimension(dim_x2_1, dim_x2_2, dim_x2_3, dim_x2_4), intent(in) :: x2
   integer(c_int), dimension(dim_n1), intent(in) :: n1
   integer(c_int), dimension(dim_n2), intent(in) :: n2
   integer(c_int), dimension(dim_nneigh1_1, dim_nneigh1_2), intent(in) :: nneigh1
   integer(c_int), dimension(dim_nneigh2_1, dim_nneigh2_2), intent(in) :: nneigh2
   double precision, dimension(dim_pd_1, dim_pd_2), intent(in) :: pd
   double precision, dimension(dim_parameters_1, dim_parameters_2), intent(in) :: parameters

   integer(c_int), intent(in), value :: verbose, nm1, nm2, na1, nsigmas, order, kernel_idx, alchemy
   double precision, intent(in), value :: t_width, d_width, cut_start, cut_distance
   double precision, intent(in), value :: distance_scale, angular_scale, two_body_power, three_body_power

   double precision, dimension(nsigmas, nm1, na1), intent(out) :: kernels

   logical :: verbose_log, alchemy_log

   verbose_log = (verbose /= 0)
   alchemy_log = (alchemy /= 0)

   call fget_atomic_local_kernels_fchl(x1, x2, verbose_log, n1, n2, nneigh1, nneigh2, &
       & nm1, nm2, na1, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy_log, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels)

end subroutine fget_atomic_local_kernels_fchl_wrapper
