module array_merge

  implicit none

  contains

!===============================================================================
! MERGE combines two arrays in to one longer array, maintaining the sorted order
!===============================================================================

    subroutine merge(a, b, result)
      real(8), target, intent(in)    :: a(:)
      real(8), target, intent(in)    :: b(:)
      real(8), allocatable, intent(inout) :: result(:)

      real(8), allocatable :: merged(:)
      integer :: ndata1, ndata2, nab
      integer :: idata1, idata2, ires
      logical :: no_exit = .true.
      real(8), pointer :: data1(:), data2(:)

      if (a(size(a)) > b(size(b))) then
        data1 => b
        data2 => a
      else
        data1 => a
        data2 => b
      end if

      ndata1 = size(data1)
      ndata2 = size(data2)
      nab = ndata1 + ndata2
      allocate(merged(nab))

      idata1 = 1
      idata2 = 1
      do ires = 1, nab
        if (idata1 <= ndata1 .and. idata2 <= ndata2) then
          if (data1(idata1) < data2(idata2)) then
            ! take data2 info, unless it is zero
            ! a zero Ein results in zero scattering, which is not a useful
            ! point to interpolate to.  MC codes will extrapolate in this case
            ! from the bottom-two points. This is more desirable than interpolating
            ! between 0 and the next highest Ein.
            ! Even more desirable, perhaps, would be interpolating between
            ! a suitably low, but non-zero Ein, and the next highest Ein
            if (data1(idata1) == 0.0_8) then
              merged(ires) = data2(idata2) * 1.0E-3_8
            else
              merged(ires) = data1(idata1)
            end if
            idata1 = idata1 + 1
          else if (data1(idata1) == data2(idata2)) then
            ! take data1 info, but increment both
            merged(ires) = data1(idata1)
            idata1 = idata1 + 1
            idata2 = idata2 + 1
          else
            ! take data2 info, unless it is zero
            ! a zero Ein results in zero scattering, which is not a useful
            ! point to interpolate to.  MC codes will extrapolate in this case
            ! from the bottom-two points. This is more desirable than interpolating
            ! between 0 and the next highest Ein.
            ! Even more desirable, perhaps, would be interpolating between
            ! a suitably low, but non-zero Ein, and the next highest Ein
            if (data2(idata2) == 0.0_8) then
              merged(ires) = data1(idata1) * 1.0E-3_8
            else
              merged(ires) = data2(idata2)
            end if
            idata2 = idata2 + 1
          end if
        else if (idata1 <= ndata1) then
          ! There are more data1 data than data2s
          ! Take an data1 and then stop building array
          merged(ires) = data1(idata1)
          idata1 = idata1 + 1
          no_exit = .false.
          exit
        else if (idata2 <= ndata2) then
          ! There are more data2s than data1 data
          merged(ires) = data2(idata2)
          idata2 = idata2 + 1
        else
          no_exit = .false.
          exit
        end if
      end do

      ! Clear result if it has values (as it will when it gets here)
      if (allocated(result)) then
        deallocate(result)
      end if

      ! Adjust ires
      if ((.not. no_exit) .or. (ires > nab)) then
        ires = ires - 1
      end if

      ! Store our results
      allocate(result(ires))
      result = merged(1: ires)
      ! Clean up
      deallocate(merged)
    end subroutine merge

    subroutine extend_grid(a)
      real(8), allocatable, intent(inout) :: a(:)
      real(8), allocatable :: temp(:)
      integer :: i, j

      allocate(temp(4*size(a) - 1))
      j = 1
      do i = 1, size(a) - 1
        temp(j) = a(i)
        temp(j + 1) = 0.25_8 * (a(i + 1) + a(i))
        temp(j + 2) = 0.5_8 * (a(i + 1) + a(i))
        temp(j + 3) = 0.75_8 * (a(i + 1) + a(i))
        j = j + 4
      end do
      temp(size(temp)) = a(size(a))

      deallocate(a)
      allocate(a(size(temp)))
      a = temp
      deallocate(temp)

    end subroutine extend_grid

end module array_merge
