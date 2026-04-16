module ggf100k_spline
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  integer, parameter :: cp = kind(1.0_real64)
  private
  public :: interv, bspline

contains

  subroutine interv(tknts, time, nspl, nleft)
    real(cp), intent(in) :: tknts(:)
    real(cp), intent(in) :: time
    integer, intent(in) :: nspl
    integer, intent(out) :: nleft
    integer :: n

    nleft = -1
    if (time < tknts(4) .or. time > tknts(nspl + 1)) return

    do n = 5, nspl + 1
      if (time <= tknts(n)) then
        nleft = n - 1
        return
      end if
    end do
  end subroutine interv

  subroutine bspline(tknts, t, jorder, nleft, spl)
    real(cp), intent(in) :: tknts(:)
    real(cp), intent(in) :: t
    integer, intent(in) :: jorder
    integer, intent(in) :: nleft
    real(cp), intent(out) :: spl(4)

    real(cp) :: deltal(4), deltar(4)
    real(cp) :: saved, term
    integer :: i, j

    spl = 0.0_cp
    deltal = 0.0_cp
    deltar = 0.0_cp
    spl(1) = 1.0_cp

    do j = 1, jorder - 1
      deltar(j) = tknts(nleft + j) - t
      deltal(j) = t - tknts(nleft + 1 - j)
      saved = 0.0_cp

      do i = 1, j
        term = spl(i) / (deltar(i) + deltal(j + 1 - i))
        spl(i) = saved + deltar(i) * term
        saved = deltal(j + 1 - i) * term
      end do

      spl(j + 1) = saved
    end do
  end subroutine bspline

end module ggf100k_spline


program get_coeffs
  use, intrinsic :: iso_fortran_env, only: real64
  use hdf5
  use ggf100k_spline, only: interv, bspline
  implicit none

  integer, parameter :: cp = kind(1.0_real64)
  integer, parameter :: jord = 4
  real(cp), parameter :: t_start = -98050.0_cp
  real(cp), parameter :: t_end   =   1950.0_cp

  integer :: lmax, n, nm, nspl, nrows, ntimes
  integer :: io, model_unit
  integer :: i, j, k, l, m, nleft, it, row
  real(cp) :: time, tstartin, tendin, dt
  real(cp), allocatable :: g(:), gt(:,:), tknts(:), spl(:)
  real(cp), allocatable :: times(:), g_out(:,:), h_out(:,:)
  integer, allocatable :: degrees(:), orders(:)
  character(len=*), parameter :: model_file = 'GGF100k'
  character(len=*), parameter :: out_file   = 'GGF100k_coefs.h5'

  ! HDF5 identifiers
  integer(hid_t) :: file_id, dspace_id, dset_id
  integer(hsize_t) :: dims1(1), dims2(2)
  integer :: hdferr

  ! ---- read model ----
  open (newunit=model_unit, file=model_file, status='old', action='read', iostat=io)
  if (io /= 0) error stop 'Could not open model file GGF100k.'
  read (model_unit, *) tstartin, tendin
  read (model_unit, *) lmax, nm, nspl

  n = lmax * (lmax + 2)
  allocate(g(n), gt(n, nspl), tknts(nspl + 4), spl(4))

  read (model_unit, *) (tknts(i), i = 1, nspl + 4)
  read (model_unit, *) gt
  close (model_unit)

  ! ---- build time array (500 equally spaced points) ----
  dt = 200.0_cp
  ntimes = (t_end - t_start) / dt + 1
  allocate(times(ntimes))
  do it = 1, ntimes
    times(it) = t_start + real(it - 1, cp) * dt
  end do

  ! ---- build degree / order index (one entry per l,m pair) ----
  nrows = lmax * (lmax + 3) / 2          ! sum_{l=1}^{lmax} (l+1)
  allocate(degrees(nrows), orders(nrows))
  allocate(g_out(nrows, ntimes), h_out(nrows, ntimes))

  row = 0
  do l = 1, lmax
    do m = 0, l
      row = row + 1
      degrees(row) = l
      orders(row)  = m
    end do
  end do

  ! ---- evaluate coefficients at every time step ----
  do it = 1, ntimes
    time = times(it)

    call interv(tknts, time, nspl, nleft)
    if (nleft < 4 .or. nleft > nspl) then
      write(*, '(a,f10.1,a)') 'Warning: skipping time ', time, ' (outside knot range)'
      g_out(:, it) = 0.0_cp
      h_out(:, it) = 0.0_cp
      cycle
    end if

    call bspline(tknts, time, jord, nleft, spl)

    do k = 1, n
      g(k) = 0.0_cp
      do j = 1, 4
        g(k) = g(k) + spl(j) * gt(k, j + nleft - 4)
      end do
    end do

    ! pack into output arrays following the original (l,m) layout
    row = 0
    do l = 1, lmax
      row = row + 1
      g_out(row, it) = g(l * l)
      h_out(row, it) = 0.0_cp
      do m = 1, l
        row = row + 1
        k = l * l + 2 * m - 1
        g_out(row, it) = g(k)
        h_out(row, it) = g(k + 1)
      end do
    end do
  end do

  ! ---- write HDF5 file ----
  call h5open_f(hdferr)
  call h5fcreate_f(out_file, H5F_ACC_TRUNC_F, file_id, hdferr)

  ! times (ntimes)
  dims1(1) = ntimes
  call h5screate_simple_f(1, dims1, dspace_id, hdferr)
  call h5dcreate_f(file_id, 'times', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, times, dims1, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5sclose_f(dspace_id, hdferr)

  ! degree (nrows)
  dims1(1) = nrows
  call h5screate_simple_f(1, dims1, dspace_id, hdferr)
  call h5dcreate_f(file_id, 'degree', H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferr)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, degrees, dims1, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5sclose_f(dspace_id, hdferr)

  ! order (nrows)
  dims1(1) = nrows
  call h5screate_simple_f(1, dims1, dspace_id, hdferr)
  call h5dcreate_f(file_id, 'order', H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferr)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, orders, dims1, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5sclose_f(dspace_id, hdferr)

  ! g_coeffs (nrows x ntimes)
  dims2 = [int(nrows, hsize_t), int(ntimes, hsize_t)]
  call h5screate_simple_f(2, dims2, dspace_id, hdferr)
  call h5dcreate_f(file_id, 'g_coeffs', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, g_out, dims2, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5sclose_f(dspace_id, hdferr)

  ! h_coeffs (nrows x ntimes)
  dims2 = [int(nrows, hsize_t), int(ntimes, hsize_t)]
  call h5screate_simple_f(2, dims2, dspace_id, hdferr)
  call h5dcreate_f(file_id, 'h_coeffs', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, h_out, dims2, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5sclose_f(dspace_id, hdferr)

  call h5fclose_f(file_id, hdferr)
  call h5close_f(hdferr)

  write(*, '(a,i0,a,a)') 'Wrote ', ntimes, ' time steps to ', out_file

end program get_coeffs
