module hdf5_interface

  use constants
  use error,           only: fatal_error
  use global
  use output,          only: write_message, time_stamp
  use string,          only: to_str

#ifdef HDF5
  use hdf5
  use h5lt
  use, intrinsic :: ISO_C_BINDING
#endif

  implicit none

#ifdef HDF5

contains

!===============================================================================
! HDF5_INITIALIZE
!===============================================================================

  subroutine hdf5_initialize()

    ! Initialize FORTRAN interface.
    call h5open_f(hdf5_err)

    ! Determine type for integer(8)
    hdf5_integer8_t = h5kind_to_type(8, H5_INTEGER_KIND)

  end subroutine hdf5_initialize

!===============================================================================
! HDF5_FINALIZE
!===============================================================================

  subroutine hdf5_finalize()
  
    ! Close FORTRAN interface.
    call h5close_f(hdf5_err)

  end subroutine hdf5_finalize

!===============================================================================
! HDF5_WRITE_HEADER
!===============================================================================

  subroutine hdf5_write_header()

    ! Write version information
    call hdf5_write_integer(hdf5_output_file, "version_major", VERSION_MAJOR)
    call hdf5_write_integer(hdf5_output_file, "version_minor", VERSION_MINOR)
    call hdf5_write_integer(hdf5_output_file, "version_release", VERSION_RELEASE)

    ! Write current date and time
    call h5ltmake_dataset_string_f(hdf5_output_file, "date_and_time", &
         time_stamp(), hdf5_err)

!     ! Write MPI information
!     call hdf5_write_integer(hdf5_output_file, "n_procs", n_procs)
!     call h5ltset_attribute_string_f(hdf5_output_file, "n_procs", &
!          "description", "Number of MPI processes", hdf5_err)

  end subroutine hdf5_write_header

  !===============================================================================
! HDF5_WRITE_INTEGER
!===============================================================================

  subroutine hdf5_write_integer(group, name, buffer)

    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    integer,        intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltmake_dataset_int_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_write_integer

!===============================================================================
! HDF5_WRITE_LONG
!===============================================================================

  subroutine hdf5_write_long(group, name, buffer)

    integer(HID_T),     intent(in) :: group
    character(*),       intent(in) :: name
    integer(8), target, intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
    integer(HID_T)   :: dspace
    integer(HID_T)   :: dset
    type(c_ptr)      :: f_ptr

    ! Create dataspace and dataset
    call h5screate_simple_f(rank, dims, dspace, hdf5_err)
    call h5dcreate_f(group, name, hdf5_integer8_t, dspace, dset, hdf5_err)

    ! Write eight-byte integer
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, hdf5_integer8_t, f_ptr, hdf5_err)

    ! Close dataspace and dataset for long integer
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

  end subroutine hdf5_write_long

!===============================================================================
! HDF5_WRITE_DOUBLE
!===============================================================================

  subroutine hdf5_write_double(group, name, buffer)

    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    real(8),        intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltmake_dataset_double_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_write_double

!===============================================================================
! HDF5_READ_INTEGER
!===============================================================================

  subroutine hdf5_read_integer(group, name, buffer)

    integer(HID_T), intent(in)    :: group
    character(*),   intent(in)    :: name
    integer,        intent(inout) :: buffer

    integer          :: buffer_copy(1)
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltread_dataset_int_f(group, name, buffer_copy, dims, hdf5_err)
    buffer = buffer_copy(1)

  end subroutine hdf5_read_integer

!===============================================================================
! HDF5_READ_LONG
!===============================================================================

  subroutine hdf5_read_long(group, name, buffer)

    integer(HID_T),     intent(in)  :: group
    character(*),       intent(in)  :: name
    integer(8), target, intent(out) :: buffer

    integer(HID_T) :: dset
    type(c_ptr)    :: f_ptr

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Get pointer to buffer
    f_ptr = c_loc(buffer)

    ! Read data from dataset
    call h5dread_f(dset, hdf5_integer8_t, f_ptr, hdf5_err)

    ! Close dataset
    call h5dclose_f(dset, hdf5_err)

  end subroutine hdf5_read_long

!===============================================================================
! HDF5_READ_DOUBLE
!===============================================================================

  subroutine hdf5_read_double(group, name, buffer)

    integer(HID_T), intent(in)  :: group
    character(*),   intent(in)  :: name
    real(8),        intent(out) :: buffer

    real(8)          :: buffer_copy(1)
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltread_dataset_double_f(group, name, buffer_copy, dims, hdf5_err)
    buffer = buffer_copy(1)

  end subroutine hdf5_read_double

#endif

end module hdf5_interface
