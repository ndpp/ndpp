module thin

  implicit none

  contains

!===============================================================================
! THIN_GRID Removes the number of f(x) points needed to reduce memory footprint
! while maintaining linear interpolation accuracy to a supplied tolerance.
! This is currently designed to work with the 3D scatt_mat arrays from ndpp.F90
! as this is the only real need at the moment.  This routine writes the results
! in place.  Thin_grid_one works on simply scatt_mat while thin_grid_two works
! on both scatt_mat and nuscatt_mat.  Thin_grid calls the correct function,
! depending on the arguments.
!===============================================================================

    subroutine thin_grid(xout, yout, tokeep, tol, compression, maxerr, yout2, yout3)
      real(8), allocatable, intent(inout) :: xout(:)      ! Resultant x grid
      real(8), allocatable, intent(inout) :: yout(:,:,:)  ! Resultant y values
      real(8), allocatable, intent(in)    :: tokeep(:)    ! x points to keep
      real(8), intent(in)                 :: tol          ! Desired fractional error to maintain
      real(8), intent(inout)              :: compression  ! Data reduction fraction)
      real(8), intent(inout)              :: maxerr       ! Maximum error due to compression
      ! Following are y grids which are also to be thinned the same as
      ! yout
      real(8), allocatable, optional, intent(inout) :: yout2(:,:,:)
      real(8), allocatable, optional, intent(inout) :: yout3(:)

      if (present(yout2)) then
        if (allocated(yout2)) then
          if (present(yout3)) then
            if (allocated(yout3)) then
              call thin_grid_three(xout, yout, yout2, yout3, tokeep, tol, compression, maxerr)
            end if
          else
            call thin_grid_two(xout, yout, yout2, tokeep, tol, compression, maxerr)
          end if
        else
          call thin_grid_one(xout, yout, tokeep, tol, compression, maxerr)
        end if
      else
        call thin_grid_one(xout, yout, tokeep, tol, compression, maxerr)
      end if

    end subroutine thin_grid

!===============================================================================
! THIN_GRID_ONE implements thin-grid for only scatt_mat
!===============================================================================

    subroutine thin_grid_one(xout, yout, tokeep, tol, compression, maxerr)
      real(8), allocatable, intent(inout) :: xout(:)     ! Resultant x grid
      real(8), allocatable, intent(inout) :: yout(:,:,:) ! Resultant y values
      real(8), allocatable, intent(in)    :: tokeep(:)    ! x points to keep
      real(8), intent(in)                 :: tol         ! Desired fractional error to maintain
      real(8), intent(out)                :: compression ! Data reduction fraction
      real(8), intent(inout)              :: maxerr      ! Maximum error due to compression

      real(8), allocatable :: xin(:)      ! Incoming x grid
      real(8), allocatable :: yin(:,:,:)  ! Incoming y values
      integer :: i, j, k, klo, khi
      integer :: all_ok
      real(8) :: x1, y1, x2, y2, x, y, testval
      integer :: num_keep, remove_it
      real(8) :: initial_size
      real(8) :: error
      real(8) :: x_frac

      initial_size = real(size(xout), 8)

      allocate(xin(size(xout)))
      xin = xout
      allocate(yin(size(yout,dim=1),size(yout,dim=2),size(yout,dim=3)))
      yin = yout

      all_ok = size(yin,dim=1) * size(yin,dim=2)
      maxerr = 0.0_8

      ! This loop will step through each entry in dim==3 and check to see if
      ! all of the values in other 2 dims can be replaced with linear interp.
      ! If not, the value will be saved to a new array, if so, it will be
      ! skipped.

      xout = 0.0_8
      yout = 0.0_8

      ! Keep first point's data
      xout(1) = xin(1)
      yout(:,:,1) = yin(:,:,1)

      ! Initialize data
      num_keep = 1
      klo = 1
      khi = 3
      k = 2
      do while (khi <= size(xin))
        remove_it = 0
        x1 = xin(klo)
        x2 = xin(khi)
        x  = xin(k)
        ! x_frac = 1.0_8 / (x2 - x1) * (x - x1)  ! Linear interp.
        x_frac = 1.0_8 / log(x2 / x1) * log(x / x1)  ! Log interp.

        ! Check to see if this x is in the tokeep list, if its not, then check
        ! it for removal. Otherwise, it stays. This is accomplished by leaving
        ! remove_it as 0, entering the else portion of if(remove_it==all_ok)
        if (.not. any(tokeep == x)) then
          do i = 1, size(yin, dim=1)
            do j = 1, size(yin, dim=2)
              y1 = yin(i,j,klo)
              y2 = yin(i,j,khi)
              y  = yin(i,j,k)

              testval = y1 + (y2 - y1) * x_frac
              error = abs(testval - y)
              if (y /= 0.0_8) then
                error = error / y
              end if
              if (error <= tol) then
                remove_it = remove_it + 1
                if (error > maxerr) then
                  maxerr = abs(testval - y)
                end if
              end if
            end do
          end do
        end if
        ! Now place the point in to the proper bin and advance iterators.
        if (remove_it == all_ok) then
          ! Then don't put it in the new grid but advance iterators
          k = k + 1
          khi = khi + 1
        else
          ! Put it in new grid and advance iterators accordingly
          num_keep = num_keep + 1
          xout(num_keep) = xin(k)
          yout(:,:,num_keep) = yin(:,:,k)
          klo = k
          k = k + 1
          khi = khi + 1
        end if
      end do
      ! Save the last point's data
      num_keep = num_keep + 1
      xout(num_keep) = xin(size(xin))
      yout(:,:,num_keep) = yin(:,:,size(xin))

      ! Finally, xout and yout were sized to match xin and yin since we knew
      ! they would be no larger than those.  Now we must resize these arrays
      ! and copy only the useful data in. Will use xin/yin for temp arrays.
      xin = xout(1:num_keep)
      yin = yout(:,:,1:num_keep)

      deallocate(xout)
      deallocate(yout)
      allocate(xout(num_keep))
      allocate(yout(size(yin,dim=1),size(yin,dim=2),num_keep))

      xout = xin(1:num_keep)
      yout = yin(:,:,1:num_keep)

      ! Clean up
      deallocate(xin)
      deallocate(yin)


      compression = (initial_size - real(size(xout),8)) / initial_size

    end subroutine thin_grid_one

!===============================================================================
! THIN_GRID_TWO implements thin-grid for scatt_mat and nuscatt_mat
!===============================================================================

    subroutine thin_grid_two(xout, yout, yout2, tokeep, tol, compression, maxerr)
      real(8), allocatable, intent(inout) :: xout(:)      ! Resultant x grid
      real(8), allocatable, intent(inout) :: yout(:,:,:)  ! Resultant y values
      real(8), allocatable, intent(inout) :: yout2(:,:,:) ! Resultant y values
      real(8), allocatable, intent(in)    :: tokeep(:)    ! x points to keep
      real(8), intent(in)                 :: tol          ! Desired fractional error to maintain
      real(8), intent(out)                :: compression  ! Data reduction fraction
      real(8), intent(inout)              :: maxerr       ! Maximum error due to compression

      real(8), allocatable :: xin(:)       ! Incoming x grid
      real(8), allocatable :: yin(:,:,:)   ! Incoming y values
      real(8), allocatable :: yin2(:,:,:)  ! Incoming y values
      integer :: i, j, k, klo, khi
      integer :: all_ok
      real(8) :: x1, y1, x2, y2, x, y, testval
      integer :: num_keep, remove_it
      real(8) :: initial_size
      real(8) :: error
      real(8) :: x_frac

      initial_size = real(size(xout), 8)

      allocate(xin(size(xout)))
      xin = xout
      allocate(yin(size(yout,dim=1),size(yout,dim=2),size(yout,dim=3)))
      yin = yout
      allocate(yin2(size(yout2,dim=1),size(yout2,dim=2),size(yout2,dim=3)))
      yin2 = yout2

      all_ok = size(yin,dim=1) * size(yin,dim=2) + &
        size(yin2,dim=1) * size(yin2,dim=2)
      maxerr = 0.0_8

      ! This loop will step through each entry in dim==3 and check to see if
      ! all of the values in other 2 dims can be replaced with linear interp.
      ! If not, the value will be saved to a new array, if so, it will be
      ! skipped.

      xout = 0.0_8
      yout = 0.0_8
      yout2 = 0.0_8

      ! Keep first point's data
      xout(1) = xin(1)
      yout(:,:,1) = yin(:,:,1)
      yout2(:,:,1) = yin2(:,:,1)

      ! Initialize data
      num_keep = 1
      klo = 1
      khi = 3
      k = 2
      do while (khi <= size(xin))
        remove_it = 0
        x1 = xin(klo)
        x2 = xin(khi)
        x  = xin(k)
        ! x_frac = 1.0_8 / (x2 - x1) * (x - x1)  ! Linear interp.
        x_frac = 1.0_8 / log(x2 / x1) * log(x / x1)  ! Log interp.

        ! Check to see if this x is in the tokeep list, if its not, then check
        ! it for removal. Otherwise, it stays. This is accomplished by leaving
        ! remove_it as 0, entering the else portion of if(remove_it==all_ok)
        if (.not. any(tokeep == x)) then
          do i = 1, size(yin, dim=1)
            do j = 1, size(yin, dim=2)
              ! First check yin
              y1 = yin(i,j,klo)
              y2 = yin(i,j,khi)
              y  = yin(i,j,k)
              testval = y1 + (y2-y1) * x_frac
              error = abs(testval - y)
              if (y /= 0.0_8) then
                error = error / y
              end if
              if (error <= tol) then
                remove_it = remove_it + 1
                if (error > maxerr) then
                  maxerr = abs(testval - y)
                end if
              end if
              ! And now check yin2
              y1 = yin2(i,j,klo)
              y2 = yin2(i,j,khi)
              y  = yin2(i,j,k)
              testval = y1 + (y2-y1) * x_frac
              error = abs(testval - y)
              if (y /= 0.0_8) then
                error = error / y
              end if
              if (error <= tol) then
                remove_it = remove_it + 1
                if (error > maxerr) then
                  maxerr = abs(testval - y)
                end if
              end if
            end do
          end do
        end if
        ! Now place the point in to the proper bin and advance iterators.
        if (remove_it == all_ok) then
          ! Then don't put it in the new grid but advance iterators
          k = k + 1
          khi = khi + 1
        else
          ! Put it in new grid and advance iterators accordingly
          num_keep = num_keep + 1
          xout(num_keep) = xin(k)
          yout(:,:,num_keep) = yin(:,:,k)
          yout2(:,:,num_keep) = yin2(:,:,k)
          klo = k
          k = k + 1
          khi = khi + 1
        end if
      end do
      ! Save the last point's data
      num_keep = num_keep + 1
      xout(num_keep) = xin(size(xin))
      yout(:,:,num_keep) = yin(:,:,size(xin))
      yout2(:,:,num_keep) = yin2(:,:,size(xin))

      ! Finally, xout and yout were sized to match xin and yin since we knew
      ! they would be no larger than those.  Now we must resize these arrays
      ! and copy only the useful data in. Will use xin/yin for temp arrays.
      xin = xout(1:num_keep)
      yin = yout(:,:,1:num_keep)
      yin2 = yout2(:,:,1:num_keep)

      deallocate(xout)
      deallocate(yout)
      deallocate(yout2)
      allocate(xout(num_keep))
      allocate(yout(size(yin,dim=1),size(yin,dim=2),num_keep))
      allocate(yout2(size(yin2,dim=1),size(yin2,dim=2),num_keep))

      xout = xin(1:num_keep)
      yout = yin(:,:,1:num_keep)
      yout2 = yin2(:,:,1:num_keep)

      ! Clean up
      deallocate(xin)
      deallocate(yin)
      deallocate(yin2)

      compression = (initial_size - real(size(xout),8)) / initial_size

    end subroutine thin_grid_two

!===============================================================================
! THIN_GRID_THREE implements thin-grid for scatt_mat and nuscatt_mat and the
! normalization terms
!===============================================================================

    subroutine thin_grid_three(xout, yout, yout2, yout3, tokeep, tol, compression, maxerr)
      real(8), allocatable, intent(inout) :: xout(:)      ! Resultant x grid
      real(8), allocatable, intent(inout) :: yout(:,:,:)  ! Resultant y values
      real(8), allocatable, intent(inout) :: yout2(:,:,:) ! Resultant y values
      real(8), allocatable, intent(inout) :: yout3(:)     ! Resultant y values
      real(8), allocatable, intent(in)    :: tokeep(:)    ! x points to keep
      real(8), intent(in)                 :: tol          ! Desired fractional error to maintain
      real(8), intent(out)                :: compression  ! Data reduction fraction
      real(8), intent(inout)              :: maxerr       ! Maximum error due to compression

      real(8), allocatable :: xin(:)       ! Incoming x grid
      real(8), allocatable :: yin(:,:,:)   ! Incoming y values
      real(8), allocatable :: yin2(:,:,:)  ! Incoming y values
      real(8), allocatable :: yin3(:)  ! Incoming y values
      integer :: i, j, k, klo, khi
      integer :: all_ok
      real(8) :: x1, y1, x2, y2, x, y, testval
      integer :: num_keep, remove_it
      real(8) :: initial_size
      real(8) :: error
      real(8) :: x_frac

      initial_size = real(size(xout), 8)

      allocate(xin(size(xout)))
      xin = xout
      allocate(yin(size(yout,dim=1),size(yout,dim=2),size(yout,dim=3)))
      yin = yout
      allocate(yin2(size(yout2,dim=1),size(yout2,dim=2),size(yout2,dim=3)))
      yin2 = yout2
      allocate(yin3(size(yout3,dim=1)))
      yin3 = yout3

      all_ok = size(yin,dim=1) * size(yin,dim=2) + &
        size(yin2,dim=1) * size(yin2,dim=2) + 1
      maxerr = 0.0_8

      ! This loop will step through each entry in dim==3 and check to see if
      ! all of the values in other 2 dims can be replaced with linear interp.
      ! If not, the value will be saved to a new array, if so, it will be
      ! skipped.

      xout = 0.0_8
      yout = 0.0_8
      yout2 = 0.0_8

      ! Keep first point's data
      xout(1) = xin(1)
      yout(:,:,1) = yin(:,:,1)
      yout2(:,:,1) = yin2(:,:,1)
      yout3(1) = yin3(1)

      ! Initialize data
      num_keep = 1
      klo = 1
      khi = 3
      k = 2
      do while (khi <= size(xin))
        remove_it = 0
        x1 = xin(klo)
        x2 = xin(khi)
        x  = xin(k)
        ! x_frac = 1.0_8 / (x2 - x1) * (x - x1)  ! Linear interp.
        x_frac = 1.0_8 / log(x2 / x1) * log(x / x1)  ! Log interp.

        ! Check to see if this x is in the tokeep list, if its not, then check
        ! it for removal. Otherwise, it stays. This is accomplished by leaving
        ! remove_it as 0, entering the else portion of if(remove_it==all_ok)
        if (.not. any(tokeep == x)) then
          do i = 1, size(yin, dim=1)
            do j = 1, size(yin, dim=2)
              ! First check yin
              y1 = yin(i,j,klo)
              y2 = yin(i,j,khi)
              y  = yin(i,j,k)
              testval = y1 + (y2-y1) * x_frac
              error = abs(testval - y)
              if (y /= 0.0_8) then
                error = error / y
              end if
              if (error <= tol) then
                remove_it = remove_it + 1
                if (error > maxerr) then
                  maxerr = abs(testval - y)
                end if
              end if
              ! And now check yin2
              y1 = yin2(i,j,klo)
              y2 = yin2(i,j,khi)
              y  = yin2(i,j,k)
              testval = y1 + (y2-y1) * x_frac
              error = abs(testval - y)
              if (y /= 0.0_8) then
                error = error / y
              end if
              if (error <= tol) then
                remove_it = remove_it + 1
                if (error > maxerr) then
                  maxerr = abs(testval - y)
                end if
              end if
            end do
          end do
          ! And now check yin3, which only has one dimension, not 3 like others
          y1 = yin3(klo)
          y2 = yin3(khi)
          y  = yin3(k)
          testval = y1 + (y2-y1) * x_frac
          error = abs(testval - y)
          if (y /= 0.0_8) then
            error = error / y
          end if
          if (error <= tol) then
            remove_it = remove_it + 1
            if (error > maxerr) then
              maxerr = abs(testval - y)
            end if
          end if
        end if
        ! Now place the point in to the proper bin and advance iterators.
        if (remove_it == all_ok) then
          ! Then don't put it in the new grid but advance iterators
          k = k + 1
          khi = khi + 1
        else
          ! Put it in new grid and advance iterators accordingly
          num_keep = num_keep + 1
          xout(num_keep) = xin(k)
          yout(:,:,num_keep) = yin(:,:,k)
          yout2(:,:,num_keep) = yin2(:,:,k)
          yout3(num_keep) = yin3(k)
          klo = k
          k = k + 1
          khi = khi + 1
        end if
      end do
      ! Save the last point's data
      num_keep = num_keep + 1
      xout(num_keep) = xin(size(xin))
      yout(:,:,num_keep) = yin(:,:,size(xin))
      yout2(:,:,num_keep) = yin2(:,:,size(xin))
      yout3(num_keep) = yin3(size(xin))

      ! Finally, xout and yout were sized to match xin and yin since we knew
      ! they would be no larger than those.  Now we must resize these arrays
      ! and copy only the useful data in. Will use xin/yin for temp arrays.
      xin = xout(1:num_keep)
      yin = yout(:,:,1:num_keep)
      yin2 = yout2(:,:,1:num_keep)
      yin3 = yout3(1:num_keep)

      deallocate(xout)
      deallocate(yout)
      deallocate(yout2)
      deallocate(yout3)
      allocate(xout(num_keep))
      allocate(yout(size(yin,dim=1),size(yin,dim=2),num_keep))
      allocate(yout2(size(yin2,dim=1),size(yin2,dim=2),num_keep))
      allocate(yout3(num_keep))

      xout = xin(1:num_keep)
      yout = yin(:,:,1:num_keep)
      yout2 = yin2(:,:,1:num_keep)
      yout3 = yin3(1:num_keep)

      ! Clean up
      deallocate(xin)
      deallocate(yin)
      deallocate(yin2)
      deallocate(yin3)

      compression = (initial_size - real(size(xout),8)) / initial_size

    end subroutine thin_grid_three


end module thin