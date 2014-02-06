module chi

  use ace_header
  use array_merge,    only: merge
  use chidata_header, only: ChiData
  use constants
  use error,          only: fatal_error, warning
  use global,         only: message
  use output,         only: write_message, print_ascii_array
  use string,         only: to_str

#ifdef HDF5
  use hdf5_interface
#endif

  implicit none

contains

!===============================================================================
! CALC_CHI Calculates the group-wise chi values for a given nuclide.
!===============================================================================
  subroutine calc_chi(nuc, E_bins, E_grid, chi_total, chi_prompt, chi_delay)
    type(Nuclide),pointer,intent(in)  :: nuc        ! Working Nuclide
    real(8),intent(in)                :: E_bins(:)  ! Energy group boundaries
    real(8),allocatable,intent(inout) :: E_grid(:)  ! Incoming energy grid
    real(8),allocatable,intent(inout) :: chi_total(:,:)   ! Total Chi Values
    real(8),allocatable,intent(inout) :: chi_prompt(:,:)  ! Prompt Chi Values
    real(8),allocatable,intent(inout) :: chi_delay(:,:,:) ! Delayed Chi Values

    type(DistEnergy),  pointer :: edist => null()
    type(Reaction), pointer    :: rxn => null()
    integer :: i_rxn, i, iE, num_fiss, groups
    real(8) :: Ein
    type(ChiData), allocatable :: prompt_data(:)
    type(ChiData), allocatable :: delay_data(:)
    real(8)              :: beta, prob, norm
    real(8), allocatable :: temp_E_grid(:), chi_p(:)
    real(8), pointer :: sigma(:) => null()

    groups = size(E_bins) - 1
    allocate(chi_p(groups))

    ! This routine will calculate chi on each energy grid for all the reactions
    ! and then combine them on to a union grid at the end.

    ! First we have to find the maximum number of dimensions and allocate the grids
    ! To do this, loop through the n_fission reactions, and count the number of nested
    ! energy distributions for each. This is an issue for at least Pu-240 in ENDF-7.
    num_fiss = 0
    do i = 1, nuc % n_fission
      rxn => nuc % reactions(nuc % index_fission(i))
      edist => rxn % edist
      num_fiss = num_fiss + 1
      do
        if (associated(edist % next)) then
          num_fiss = num_fiss + 1
          edist => edist % next
        else
          exit
        end if
      end do
    end do
    allocate(prompt_data(num_fiss))
    allocate(delay_data(nuc % n_precursor))

    ! Now we can go through and initialize each of these.
    ! Start with prompt fissions
    i_rxn = 0
    do i = 1, nuc % n_fission
      i_rxn = i_rxn + 1
      rxn => nuc % reactions(nuc % index_fission(i))
      edist => rxn % edist
      if (rxn % MT == N_FISSION) then
        sigma => nuc % fission
      else
        sigma => rxn % sigma
      end if
      call prompt_data(i_rxn) % init(nuc, rxn, sigma, edist, E_bins, .false.)
      do
        if (associated(edist % next)) then
          i_rxn = i_rxn + 1
          edist => edist % next
          call prompt_data(i_rxn) % init(nuc, rxn, sigma, edist, E_bins, .false.)
        else
          exit
        end if
      end do
    end do

    ! Move to delayed
    do i = 1, nuc % n_precursor
      edist => nuc % nu_d_edist(i)
      call delay_data(i) % init(nuc, rxn, sigma, edist, E_bins, .true., i)
    end do

    ! Build the combined energy grid
    ! First lets do for the prompt data
    allocate(E_grid(prompt_data(1) % NE))
    E_grid = prompt_data(1) % E_grid
    do i = 2, size(prompt_data)
      call merge(E_grid, prompt_data(i) % E_grid, temp_E_grid)
      deallocate(E_grid)
      allocate(E_grid(size(temp_E_grid)))
      E_grid = temp_E_grid
    end do

    ! Move on to the delayed data, incorporating them in to this grid too
    do i = 1, size(delay_data)
      call merge(E_grid, delay_data(i) % E_grid, temp_E_grid)
      deallocate(E_grid)
      allocate(E_grid(size(temp_E_grid)))
      E_grid = temp_E_grid
    end do
    deallocate(temp_E_grid)

    ! Now we can step through each incoming energy point in E_grid,
    ! and calculate the chi, beta, and probability values for each
    ! And init our final values to pass back
    allocate(chi_total(groups, size(E_grid)))
    allocate(chi_prompt(groups, size(E_grid)))
    allocate(chi_delay(groups, size(E_grid), nuc % n_precursor))
    chi_total = ZERO
    chi_prompt = ZERO

    do iE = 1, size(E_grid)
      Ein = E_grid(iE)
      ! Prompt
      do i = 1, size(prompt_data)
        chi_p = prompt_data(i) % integrate(Ein)
        prob  = prompt_data(i) % prob(Ein)
        beta  = prompt_data(i) % beta(Ein)
        chi_total(:,iE) = chi_total(:,iE) + prob * beta * chi_p
        chi_prompt(:,iE) = chi_prompt(:,iE) + prob * chi_p
      end do

      ! Delayed
      do i = 1, size(delay_data)
        chi_delay(:,iE,i) = delay_data(i) % integrate(Ein)
        prob = delay_data(i) % prob(Ein)
        beta = delay_data(i) % beta(Ein)
        chi_total(:,iE) = chi_total(:,iE) + prob * beta * chi_delay(:,iE,i)
      end do

      ! Finally, normalize each of these to ensure the groups PDFs sum to 1.0
      norm = sum(chi_total(:,iE))
      if (norm > ZERO) then
        chi_total(:,iE) = chi_total(:,iE) / norm
      end if
      norm = sum(chi_prompt(:,iE))
      if (norm > ZERO) then
        chi_prompt(:,iE) = chi_prompt(:,iE) / norm
      end if
      do i = 1, size(delay_data)
        norm = sum(chi_delay(:,iE,i))
        if (norm > ZERO) then
          chi_delay(:,iE,i) = chi_delay(:,iE,i) / norm
        end if
      end do
    end do

    ! Clean up
    do i = 1, size(prompt_data)
      call prompt_data(i) % clear()
    end do

    do i = 1, size(delay_data)
      call delay_data(i) % clear()
    end do
  end subroutine calc_chi

!===============================================================================
! PRINT_CHI prints the chi data to the specified output file
! in the specified format.
!===============================================================================

  subroutine print_chi(name, lib_format, E_grid, chi_t, chi_p, chi_d)
    character(len=*),     intent(in) :: name         ! (hdf5 specific) name of group
    integer,              intent(in) :: lib_format   ! Library output type
    real(8), allocatable, intent(in) :: E_grid(:)    ! Energy grid
    real(8), allocatable, intent(in) :: chi_t(:,:)   ! Unionized Total Chi
    real(8), allocatable, intent(in) :: chi_p(:,:)   ! Unionized Prompt Chi
    real(8), allocatable, intent(in) :: chi_d(:,:,:) ! Unionized Delayed Chi

    if (lib_format == ASCII) then
      call print_chi_ascii(E_grid, chi_t, chi_p, chi_d)
    else if (lib_format == BINARY) then
      call print_chi_bin(E_grid, chi_t, chi_p, chi_d)
    else if (lib_format == HUMAN) then
      call print_chi_human(E_grid, chi_t, chi_p, chi_d)
    else if (lib_format == H5) then
      call print_chi_hdf5(name, E_grid, chi_t, chi_p, chi_d)
    end if

  end subroutine print_chi

!===============================================================================
! PRINT_CHI_ASCII prints the chi data to the specified output file
! in an ASCII format.
!===============================================================================

  subroutine print_chi_ascii(E_grid, chi_t, chi_p, chi_d)
    real(8), allocatable, intent(in) :: E_grid(:)    ! Energy grid
    real(8), allocatable, intent(in) :: chi_t(:,:)   ! Unionized Total Chi
    real(8), allocatable, intent(in) :: chi_p(:,:)   ! Unionized Prompt Chi
    real(8), allocatable, intent(in) :: chi_d(:,:,:) ! Unionized Delayed Chi

    character(MAX_LINE_LEN) :: line
    integer :: i

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format with at max 4 entries per line:
    ! <NE, num_precursors>
    ! <1:NE energies>
    ! <1:NE chi_t(g,E)
    ! <1:NE chi_p(g,E)
    !  For each precursor:
    ! <1:NE chi_d(g,precursors,E)

    ! Begin writing:

    ! <NE, num_precursors>
    line = ''
    write(line,'(I20,I20)') size(E_grid), size(chi_d, DIM=3)
    write(UNIT_NUC,'(A)') trim(line)

    ! <1:NE energies>
    call print_ascii_array(E_grid, UNIT_NUC)

    ! <1:NE chi_t(g,E)
    call print_ascii_array(pack(chi_t, .true.), UNIT_NUC)

    ! <1:NE chi_p(g,E)
    call print_ascii_array(pack(chi_p, .true.), UNIT_NUC)

    !  For each precursor:
    ! <1:NE chi_d(g,precursors,E)
    do i = 1, size(chi_d, DIM=3)
      call print_ascii_array(pack(chi_d(:,:,i), .true.), UNIT_NUC)
    end do

  end subroutine print_chi_ascii

!===============================================================================
! PRINT_CHI_HUMAN prints the chi data to the specified output file
! in a human-readable ASCII format.
!===============================================================================

  subroutine print_chi_human(E_grid, chi_t, chi_p, chi_d)
    real(8), allocatable, intent(in) :: E_grid(:)    ! Energy grid
    real(8), allocatable, intent(in) :: chi_t(:,:)   ! Unionized Total Chi
    real(8), allocatable, intent(in) :: chi_p(:,:)   ! Unionized Prompt Chi
    real(8), allocatable, intent(in) :: chi_d(:,:,:) ! Unionized Delayed Chi

    integer :: iE

    character(MAX_LINE_LEN) :: line
    integer :: i

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format with at max 4 entries per line:
    ! <NE, num_precursors>
    ! <1:NE energies>
    ! <1:NE chi_t(g,E)
    ! <1:NE chi_p(g,E)
    !  For each precursor:
    ! <1:NE chi_d(g,precursors,E)

    ! Begin writing:

    ! <NE, num_precursors>
    write(UNIT_NUC,*) 'Chi Data'
    line = ''
    write(line,'(I20,I20)') size(E_grid), size(chi_d, DIM=3)
    write(UNIT_NUC,'(A)') trim(line)

    ! <1:NE energies>
    write(UNIT_NUC,*) 'Chi Energy Grid'
    call print_ascii_array(E_grid, UNIT_NUC)

    ! <1:NE chi_t(g,E)
    write(UNIT_NUC,*) 'Chi Total'
    do iE = 1, size(E_grid)
      write(line, '(A,1PE20.12)') 'Ein=',E_grid(iE)
      write(UNIT_NUC, '(A)') trim(line)
      call print_ascii_array(chi_t(:, iE), UNIT_NUC)
    end do

    ! <1:NE chi_p(g,E)
    write(UNIT_NUC,*) 'Chi Prompt'
    do iE = 1, size(E_grid)
      write(line, '(A,1PE20.12)') 'Ein=',E_grid(iE)
      write(UNIT_NUC, '(A)') trim(line)
      call print_ascii_array(chi_p(:, iE), UNIT_NUC)
    end do

    !  For each precursor:
    ! <1:NE chi_d(g,precursors,E)
    write(UNIT_NUC,*) 'Chi Delayed'
    do i = 1, size(chi_d, DIM=3)
      write(UNIT_NUC,*) 'Precursor', i
      do iE = 1, size(E_grid)
        write(line, '(A,1PE20.12)') 'Ein=',E_grid(iE)
        write(UNIT_NUC, '(A)') trim(line)
        call print_ascii_array(chi_d(:, iE, i), UNIT_NUC)
      end do
    end do

  end subroutine print_chi_human

!===============================================================================
! PRINT_CHI_BIN prints the chi data to the specified output file
! in Fortran binary (stream) format.
!===============================================================================

  subroutine print_chi_bin(E_grid, chi_t, chi_p, chi_d)
    real(8), allocatable, intent(in) :: E_grid(:)    ! Energy grid
    real(8), allocatable, intent(in) :: chi_t(:,:)   ! Unionized Total Chi
    real(8), allocatable, intent(in) :: chi_p(:,:)   ! Unionized Prompt Chi
    real(8), allocatable, intent(in) :: chi_d(:,:,:) ! Unionized Delayed Chi

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format
    ! <NE, num_precursors>
    ! <1:NE energies>
    ! <1:NE chi_t(g,E)
    ! <1:NE chi_p(g,E)
    !  For each precursor:
    ! <1:NE chi_d(g,precursors,E)

    ! Begin writing:

    ! <NE, num_precursors>
    write(UNIT_NUC) size(E_grid)

    ! <1:NE energies>
    write(UNIT_NUC) E_grid

    ! <1:NE chi_t(g,E)
    write(UNIT_NUC) chi_t

    ! <1:NE chi_p(g,E)
    write(UNIT_NUC) chi_p

    !  For each precursor:
    ! <1:NE chi_d(g,precursors,E)
    write(UNIT_NUC) chi_d

  end subroutine print_chi_bin

!===============================================================================
! PRINT_CHI_BIN prints the chi data to the specified output file
! in Fortran binary (stream) format.
!===============================================================================

  subroutine print_chi_hdf5(name, E_grid, chi_t, chi_p, chi_d)
    character(len=*),     intent(in) :: name       ! name of group
    real(8), allocatable, intent(in) :: E_grid(:)    ! Energy grid
    real(8), allocatable, intent(in) :: chi_t(:,:)   ! Unionized Total Chi
    real(8), allocatable, intent(in) :: chi_p(:,:)   ! Unionized Prompt Chi
    real(8), allocatable, intent(in) :: chi_d(:,:,:) ! Unionized Delayed Chi

#ifdef HDF5
    ! Will follow this format
    ! <NE, num_precursors>
    ! <1:NE energies>
    ! <1:NE chi_t(g,E)
    ! <1:NE chi_p(g,E)
    !  For each precursor:
    ! <1:NE chi_d(g,precursors,E)
#endif
  end subroutine print_chi_hdf5

end module chi
