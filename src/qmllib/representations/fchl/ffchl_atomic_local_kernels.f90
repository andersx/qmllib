subroutine fget_atomic_local_kernels_fchl(nm1, nm2, na1, nsigmas, n1_size, n2_size, &
       & nneigh1_size1, nneigh1_size2, nneigh2_size1, nneigh2_size2, &
       & x1_size1, x1_size2, x1_size3, x1_size4, &
       & x2_size1, x2_size2, x2_size3, x2_size4, &
       & pd_size1, pd_size2, parameters_size1, parameters_size2, &
       & x1, x2, verbose, n1, n2, nneigh1, nneigh2, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, &
       & kernel_idx, parameters, kernels) bind(C, name="fget_atomic_local_kernels_fchl")

   use iso_c_binding
   use ffchl_module, only: scalar, get_threebody_fourier, get_twobody_weights, &
       & get_angular_norm2, get_pmax, get_ksi, init_cosp_sinp, get_selfscalar
   use ffchl_kernels, only: kernel

   implicit none

   ! Dimensions (MUST be first with value attribute for bind(C))
   integer(c_int), intent(in), value :: nm1, nm2, na1, nsigmas
   integer(c_int), intent(in), value :: n1_size, n2_size
   integer(c_int), intent(in), value :: nneigh1_size1, nneigh1_size2
   integer(c_int), intent(in), value :: nneigh2_size1, nneigh2_size2
   integer(c_int), intent(in), value :: x1_size1, x1_size2, x1_size3, x1_size4
   integer(c_int), intent(in), value :: x2_size1, x2_size2, x2_size3, x2_size4
   integer(c_int), intent(in), value :: pd_size1, pd_size2
   integer(c_int), intent(in), value :: parameters_size1, parameters_size2

   ! fchl descriptors for the training set, format (nm1,maxatoms,5,maxneighbors)
   real(c_double), dimension(x1_size1, x1_size2, x1_size3, x1_size4), intent(in) :: x1
   real(c_double), dimension(x2_size1, x2_size2, x2_size3, x2_size4), intent(in) :: x2

   ! Whether to be verbose with output (C int, not logical)
   integer(c_int), intent(in), value :: verbose

   ! List of numbers of atoms in each molecule
   integer(c_int), dimension(n1_size), intent(in) :: n1
   integer(c_int), dimension(n2_size), intent(in) :: n2

   ! Number of neighbors for each atom in each compound
   integer(c_int), dimension(nneigh1_size1, nneigh1_size2), intent(in) :: nneigh1
   integer(c_int), dimension(nneigh2_size1, nneigh2_size2), intent(in) :: nneigh2

   real(c_double), intent(in), value :: two_body_power
   real(c_double), intent(in), value :: three_body_power

   real(c_double), intent(in), value :: t_width
   real(c_double), intent(in), value :: d_width
   real(c_double), intent(in), value :: cut_start
   real(c_double), intent(in), value :: cut_distance
   integer(c_int), intent(in), value :: order
   real(c_double), intent(in), value :: distance_scale
   real(c_double), intent(in), value :: angular_scale

   ! -1.0 / sigma^2 for use in the kernel
   real(c_double), dimension(pd_size1, pd_size2), intent(in) :: pd

   integer(c_int), intent(in), value :: kernel_idx
   real(c_double), dimension(parameters_size1, parameters_size2), intent(in) :: parameters

   ! Resulting kernel matrix
   real(c_double), dimension(nsigmas, na1, nm2), intent(out) :: kernels

   ! Convert C integer to Fortran logical
   logical :: verbose_logical
   logical :: alchemy_logical

   integer(c_int), intent(in), value :: alchemy

   integer :: idx1

   ! Internal counters
   integer :: i, j
   integer :: ni, nj
   integer :: a, b

   ! Temporary variables necessary for parallelization
   double precision :: s12

   ! Pre-computed terms in the full distance matrix
   double precision, allocatable, dimension(:, :) :: self_scalar1
   double precision, allocatable, dimension(:, :) :: self_scalar2

   ! Pre-computed terms
   double precision, allocatable, dimension(:, :, :) :: ksi1
   double precision, allocatable, dimension(:, :, :) :: ksi2

   double precision, allocatable, dimension(:, :, :, :, :) :: sinp1
   double precision, allocatable, dimension(:, :, :, :, :) :: sinp2
   double precision, allocatable, dimension(:, :, :, :, :) :: cosp1
   double precision, allocatable, dimension(:, :, :, :, :) :: cosp2

   ! Value of PI at full FORTRAN precision.
   double precision, parameter :: pi = 4.0d0*atan(1.0d0)

   ! counter for periodic distance
   integer :: pmax1
   integer :: pmax2

   double precision :: ang_norm2

   integer :: maxneigh1
   integer :: maxneigh2

   ! Work kernel
   double precision, allocatable, dimension(:) :: ktmp

   ! Convert C integers to Fortran logicals
   verbose_logical = (verbose /= 0)
   alchemy_logical = (alchemy /= 0)

   allocate (ktmp(size(parameters, dim=1)))

   maxneigh1 = maxval(nneigh1)
   maxneigh2 = maxval(nneigh2)

   ang_norm2 = get_angular_norm2(t_width)

   pmax1 = get_pmax(x1, n1)
   pmax2 = get_pmax(x2, n2)

   allocate (ksi1(size(x1, dim=1), maxval(n1), maxval(nneigh1)))
   allocate (ksi2(size(x2, dim=1), maxval(n2), maxval(nneigh2)))
   call get_ksi(x1, n1, nneigh1, two_body_power, cut_start, cut_distance, verbose_logical, ksi1)
   call get_ksi(x2, n2, nneigh2, two_body_power, cut_start, cut_distance, verbose_logical, ksi2)

   allocate (cosp1(nm1, maxval(n1), pmax1, order, maxval(nneigh1)))
   allocate (sinp1(nm1, maxval(n1), pmax1, order, maxval(nneigh1)))

   call init_cosp_sinp(x1, n1, nneigh1, three_body_power, order, cut_start, cut_distance, &
       & cosp1, sinp1, verbose_logical)

   allocate (cosp2(nm2, maxval(n2), pmax2, order, maxval(nneigh2)))
   allocate (sinp2(nm2, maxval(n2), pmax2, order, maxval(nneigh2)))

   call init_cosp_sinp(x2, n2, nneigh2, three_body_power, order, cut_start, cut_distance, &
       & cosp2, sinp2, verbose_logical)

   ! Pre-calculate self-scalar terms
   allocate (self_scalar1(nm1, maxval(n1)))
   allocate (self_scalar2(nm2, maxval(n2)))
   call get_selfscalar(x1, nm1, n1, nneigh1, ksi1, sinp1, cosp1, t_width, d_width, &
        & cut_distance, order, pd, ang_norm2, distance_scale, angular_scale, alchemy_logical, verbose_logical, self_scalar1)
   call get_selfscalar(x2, nm2, n2, nneigh2, ksi2, sinp2, cosp2, t_width, d_width, &
        & cut_distance, order, pd, ang_norm2, distance_scale, angular_scale, alchemy_logical, verbose_logical, self_scalar2)

   kernels(:, :, :) = 0.0d0

   !$OMP PARALLEL DO schedule(dynamic) PRIVATE(ni,nj,idx1,s12,ktmp)
   do a = 1, nm1
      ni = n1(a)
      do i = 1, ni

         idx1 = sum(n1(:a)) - ni + i

         do b = 1, nm2
            nj = n2(b)
            do j = 1, nj

               s12 = scalar(x1(a, i, :, :), x2(b, j, :, :), &
                   & nneigh1(a, i), nneigh2(b, j), ksi1(a, i, :), ksi2(b, j, :), &
                   & sinp1(a, i, :, :, :), sinp2(b, j, :, :, :), &
                   & cosp1(a, i, :, :, :), cosp2(b, j, :, :, :), &
                   & t_width, d_width, cut_distance, order, &
                   & pd, ang_norm2, distance_scale, angular_scale, alchemy_logical)

               ktmp = 0.0d0
               call kernel(self_scalar1(a, i), self_scalar2(b, j), s12, &
                   & kernel_idx, parameters, ktmp)
               kernels(:, idx1, b) = kernels(:, idx1, b) + ktmp

            end do
         end do

      end do
   end do
   !$OMP END PARALLEL DO

   deallocate (ktmp)
   deallocate (self_scalar1)
   deallocate (self_scalar2)
   deallocate (ksi1)
   deallocate (ksi2)
   deallocate (cosp1)
   deallocate (cosp2)
   deallocate (sinp1)
   deallocate (sinp2)

end subroutine fget_atomic_local_kernels_fchl
