module ndpp_chi
  
  use ace_header
  use constants
  use dict_header
  use error,            only: fatal_error, warning
  use fission,          only: nu_total, nu_delayed
  use global,           only: nuclides, message
  use interpolation,    only: interpolate_tab1
  use output,           only: write_message, header, print_ascii_array
  use search,           only: binary_search
  use string,           only: to_str

  implicit none
  
contains
 
!===============================================================================
! CALC_CHIS Calculates the group-wise chi values for a given nuclide.
!===============================================================================
  subroutine calc_chis(nuc, energy_bins, energy_groups, chi_t_union, chi_p_union, &
    chi_d_union, E_t_union, E_p_union, E_d_union, thin_tol)
    type(Nuclide), pointer, intent(in) :: nuc            ! Nuclide
    real(8), intent(in)                :: energy_bins(:) ! Energy groups
    integer, intent(in)                :: energy_groups  ! Number of groups
    real(8), allocatable,intent(inout) :: chi_t_union(:,:) ! Unionized Total Chi Values
    real(8), allocatable,intent(inout) :: chi_p_union(:,:) ! Unionized Prompt Chi Values
    real(8), allocatable,intent(inout) :: chi_d_union(:,:) ! Unionized Delayed Chi Values
    real(8), allocatable,intent(inout) :: E_t_union(:)     ! Unionized Total Energy Values
    real(8), allocatable,intent(inout) :: E_p_union(:)     ! Unionized Prompt Energy Values
    real(8), allocatable,intent(inout) :: E_d_union(:)     ! Unionized Delayed Energy Values
    real(8), intent(in)                :: thin_tol         ! Thinning tolerance
    
    real(8) :: beta
    real(8) :: yield
    real(8) :: probability
    real(8), allocatable :: E_grid(:,:)     ! The energy grid for each fission rxn
    real(8), allocatable :: chi_grid(:,:,:) ! Chis for the global grid E values
    real(8), allocatable :: chi_temp(:)     ! Temporary Chi array
    integer, allocatable :: NE_grid(:)      ! NE values for each reaction
    integer, allocatable :: INT_grid(:,:)   ! The interpolation type grid for each fission rxn
    real(8), allocatable :: betas(:,:)        ! Delayed-n Beta at each E_in
    real(8), allocatable :: probs(:,:)      ! Probabilities of each fission reaction for each E_in
    
    type(DistEnergy),  pointer :: edist
    type(Reaction), pointer    :: rxn
    integer :: NE, NR, NE_yield, NR_yield
    integer :: lc, lc_yield
    integer :: NE_max
    integer :: num_rxn
    integer :: i_rxn
    integer :: i_E
    integer :: i, j
    real(8) :: max_err
    integer :: num_fission_tot
    
    
    ! This routine will calculate chi on each energy grid for all the reactions
    ! and then combine them on to a union grid at the end.
    
    ! First we have to find the maximum number of dimensions and allocate the grids
    ! To do this, loop through the n_fission reactions, and count the number of nested
    ! energy distributions for each. This is an issue for at least Pu-240 in ENDF-7.
    num_rxn = 0
    do i = 1, nuc % n_fission
      rxn => nuc % reactions(nuc % index_fission(i))
      edist => rxn % edist
      num_rxn = num_rxn + 1
      do
        if (associated(edist % next)) then
          num_rxn = num_rxn + 1
          edist => edist % next
        else
          exit
        end if
      end do
    end do
    num_fission_tot = num_rxn
    ! Now add in the standard fission channels and delayed channels
    num_rxn = num_rxn + nuc % n_precursor
    allocate(NE_grid(num_rxn))
    
    NE_max = 0
    i_rxn = 0
    ! Check prompt
    do i = 1, nuc % n_fission
      i_rxn = i_rxn + 1
      ! Check prompt rxn for the total number of points
      rxn => nuc % reactions(nuc % index_fission(i))
      edist => rxn % edist
      NE = get_e_grid_count(edist)
      if (NE > NE_max) NE_max = NE
      NE_grid(i_rxn) = NE
      ! Do the nested rxn channels
      do
        if (associated(edist % next)) then
          i_rxn = i_rxn + 1
          edist => edist % next
          NE = get_e_grid_count(edist)
          if (NE > NE_max) NE_max = NE
          NE_grid(i_rxn) = NE
        else
          exit
        end if
      end do
    end do
    ! Check delayed reactions for the total number of points
    do i = 1, nuc % n_precursor
      i_rxn = i_rxn + 1
      edist => nuc % nu_d_edist(i)
      NE = get_e_grid_count(edist)
      if (NE > NE_max) NE_max = NE
      NE_grid(i_rxn) = NE
    end do
    
    allocate(E_grid(NE_max,num_rxn))
    allocate(INT_grid(NE_max,num_rxn))
    allocate(chi_grid(energy_groups, NE_max, num_rxn))
    chi_grid = ZERO
    allocate(betas(NE_max, num_rxn))
    allocate(probs(NE_max, num_rxn))
    allocate(chi_temp(energy_groups))
    
    ! Calculate Chi for the prompt reaction
    ! Chi_prompt = integral[Eg-1,Eg](Chi_prompt(E)*(1-beta)*prob_of_rxn)dE
    i = 0
    do i_rxn = 1, nuc % n_fission
      i = i + 1
      rxn => nuc % reactions(nuc % index_fission(i_rxn))
      edist => rxn % edist
      ! Get the energy grid
      NR = int(edist % data(1))
      !!! This seems to not be important.
!       if (NR /= 0) write(*,*) 'NR /= 0!', i, NR, edist%law
      lc = 2 + 2*NR
      E_grid(1:NE_grid(i),i) = edist % data(lc + 1 : lc + NE_grid(i))
      
      ! Loop over energy grid and calculate
      do j = 1, NE_grid(i)
        ! Calculate beta
        beta = nu_delayed(nuc, E_grid(j, i)) / nu_total(nuc, E_grid(j, i))
        
        ! calculate probabilty of this rxn over all fission rxns of this nuclide
        if (rxn % MT == N_FISSION) then
          probability = ONE
        else
          if (E_grid(j, i) < nuc % energy(1)) then
            i_E = 1
          elseif (E_grid(j, i) >= nuc % energy(nuc % n_grid)) then
            i_E = nuc % n_grid
          else
            i_E = binary_search(nuc % energy, nuc % n_grid, E_grid(j, i))
          end if
          ! If block needed because of Pu-240 issue.
          if (i_E < rxn % threshold) then
            probability = ZERO
          else
            probability = rxn % sigma(i_E - rxn % threshold + 1) / nuc % fission(i_E)
          end if
          if (associated(edist % next) .and. edist % p_valid % n_regions > 0) then
            probability = probability * interpolate_tab1(edist % p_valid, E_grid(j, i))
          end if
        end if
        
        ! Now perform the integration
        call integrate_chi(edist, E_grid(j, i), INT_grid(j, i), energy_groups, &
          energy_bins, chi_temp)
        chi_grid(:, j, i) = chi_grid(:, j, i) + chi_temp(:)

        ! Store beta and probs
        betas(j, i) = (ONE - beta)
        probs(j, i) = probability
      end do
      
      ! Repeat for the nested loops
      do
        if (associated(edist % next)) then
          i = i + 1
          edist => edist % next
          ! Get the energy grid
          NR = int(edist % data(1))
          !!! This seems to not be important.
          if (NR /= 0) write(*,*) 'NR /= 0!', i, NR, edist%law, 'nested'
          lc = 2 + 2*NR
          E_grid(1:NE_grid(i),i) = edist % data(lc + 1 : lc + NE_grid(i))
          
          ! Loop over energy grid and calculate
          do j = 1, NE_grid(i)
            ! Calculate beta
            beta = nu_delayed(nuc, E_grid(j, i)) / nu_total(nuc, E_grid(j, i))
            
            ! calculate probabilty of this nested distribution over the others
            if (E_grid(j, i) < nuc % energy(1)) then
              i_E = 1
            elseif (E_grid(j, i) >= nuc % energy(nuc % n_grid)) then
              i_E = nuc % n_grid
            else
              i_E = binary_search(nuc % energy, nuc % n_grid, E_grid(j, i))
            end if
            ! Get the probability of law validity. If block needed because of Pu-240 issue.
            if (i_E < rxn % threshold + 1) then
              probability = ZERO
            else
              probability = rxn % sigma(i_E - rxn % threshold + 1) / nuc % fission(i_E)
            end if
            if (associated(edist % next) .and. edist % p_valid % n_regions > 0) then
              probability = probability * interpolate_tab1(edist % p_valid, E_grid(j, i))
            end if
            
            ! Now perform the integration
            call integrate_chi(edist, E_grid(j, i), INT_grid(j, i), energy_groups, &
              energy_bins, chi_temp)
            chi_grid(:, j, i) = chi_grid(:, j, i) + chi_temp(:)

            ! Store beta and probs
            betas(j, i) = (ONE - beta)
            probs(j, i) = probability
          end do
        else
          exit
        end if
      end do
    end do
    
    ! Calculate Chi for the delayed neutrons
    ! Chi_delay= integral[Eg-1,Eg](Chi_delay(E)*(beta)*prob_of_rxn)dE
    i = num_fission_tot
    lc_yield = 1
    rxn => nuc % reactions(nuc % index_fission(1))
    do i_rxn = 1, nuc % n_precursor
      i = i + 1
      ! Get the energy grid
      edist => nuc % nu_d_edist(i_rxn)
      NR = int(edist % data(1))
      if (NR /= 0) write(*,*) 'NR /= 0!'
      lc = 2 + 2*NR
      E_grid(1:NE_grid(i),i) = edist % data(lc + 1 : lc + NE_grid(i))
      ! Loop over energy grid and calculate
      do j = 1, NE_grid(i)
        ! Calculate beta
        beta = nu_delayed(nuc, E_grid(j, i)) / nu_total(nuc, E_grid(j, i))
        
        ! Get the yield data
        ! determine number of interpolation regions and energies
        NR_yield = int(nuc % nu_d_precursor_data(lc + 1))
        NE_yield = int(nuc % nu_d_precursor_data(lc + 2 + 2*NR))
        ! determine delayed neutron precursor yield for group j
        yield = interpolate_tab1(nuc % nu_d_precursor_data( &
              lc_yield + 1 :lc_yield + 2 + 2 * NR_yield + 2 * NE_yield), &
              E_grid(j, i))
        probability = yield
        
        ! Now perform the integration
        call integrate_chi(edist, E_grid(j, i), INT_grid(j, i), energy_groups, &
          energy_bins, chi_temp)
        chi_grid(:, j, i) = chi_grid(:, j, i) + chi_temp(:)
          
        ! Store beta and probs
        betas(j, i) = beta
        probs(j, i) = probability
      end do
      ! Advance lc pointer
      lc_yield = lc_yield + 2 + 2*NR_yield + 2*NE_yield + 1
    end do
    
    deallocate(chi_temp)
    
    ! Combine the total chi values to a  to a union grid
    call unionize_chi(chi_t_union, E_t_union, E_grid, chi_grid, NE_grid, &
      INT_grid, betas, probs, 1, size(NE_grid))
    betas = ONE ! betas is not needed any more.
    ! Combine the prompt chi values to a  to a union grid
    call unionize_chi(chi_p_union, E_p_union, E_grid, chi_grid, NE_grid, &
      INT_grid, betas, probs, 1, num_fission_tot)
    ! Create the delayed neutron grid, if there is delayed data
    if (num_fission_tot < size(NE_grid)) then
      call unionize_chi(chi_d_union, E_d_union, E_grid, chi_grid, NE_grid, &
        INT_grid, betas, probs, num_fission_tot + 1, size(NE_grid))
    else
       ! There were no chi values, just set to 0.0s for chi at the energy group boundaries
      allocate(chi_d_union(energy_groups, 2))
      chi_d_union = ZERO
      allocate(E_d_union(2))
      E_d_union(1) = energy_bins(1)
      E_d_union(2) = energy_bins(energy_groups + 1)
    end if
    
!     call thin_union_chi(chi_t_union, E_t_union, thin_tol, max_err)
!     call thin_union_chi(chi_p_union, E_p_union, thin_tol, max_err_p)
!     call thin_union_chi(chi_d_union, E_d_union, thin_tol, max_err_d)
    
    ! Free memory
    deallocate(E_grid)
    deallocate(chi_grid)
    deallocate(NE_grid)
    deallocate(INT_grid)
    deallocate(betas)
    deallocate(probs)
    
  end subroutine calc_chis
  
  function get_e_grid_count(edist) result(NE)
    type(DistEnergy), pointer, intent(in) :: edist
    integer                               :: NE
    integer                               :: NR
    
    ! read number of interpolation regions and incoming energies
    NR = int(edist % data(1))
    NE = int(edist % data(2 + 2*NR))
    
  end function get_e_grid_count
  
  subroutine integrate_chi(edist, E_in, INTT_in, energy_groups, energy_bins, chi)
    type(DistEnergy), pointer, intent(in) :: edist
    real(8), intent(in)                   :: E_in
    integer, intent(inout)                :: INTT_in
    integer, intent(in)                   :: energy_groups
    real(8), intent(in)                   :: energy_bins(:)
    real(8), intent(inout)                :: chi(:)
    
    integer :: NR, NE, NP, i_E, INTTp, INTT, ND
    integer :: lEout_min
    real(8) :: T              ! Fission Spectra Theta Value
    real(8) :: U              ! Restriction Energy
    real(8) :: I              ! Non-dimensional Fiss Spectra param.
    real(8) :: x, x0          ! Spectra constants
    integer :: lc             ! location in the data array
    integer :: g, g2          ! E group indices
    real(8) :: Egp1           ! upper bound of integral
    real(8) :: Eg             ! lower bound of integral
    real(8) :: Watt_a, Watt_b ! Watt spectrum values
    real(8) :: interp         ! interpolation value of the energy point
    
    chi = ZERO
    
    ! Determine which secondary energy distribution law to use
    select case (edist % law)
    case (1)
      ! =======================================================================
      ! TABULAR EQUIPROBABLE ENERGY BINS
      message = "Energy Distribution Type " // trim(to_str(edist % law)) // &
        " Not Yet Supported."
      call warning()
    case (3)
      ! =======================================================================
      ! INELASTIC LEVEL SCATTERING
      message = "Energy Distribution Type " // trim(to_str(edist % law)) // &
        " Not Yet Supported."
      call warning()
    case (4, 61)
      ! =======================================================================
      ! CONTINUOUS TABULAR DISTRIBUTION AND 
      ! CORRELATED ENERGY AND ANGLE DISTRIBUTION
      
      ! read number of interpolation regions and incoming energies
      NR  = int(edist % data(1))
      NE  = int(edist % data(2 + 2*NR))
      
      ! find energy bin and calculate interpolation factor -- if the energy is
      ! outside the range of the tabulated energies, choose the first or last
      ! bins
      lc = 2 + 2*NR
      if (E_in < edist % data(lc+1)) then
        i_E = 1
      elseif (E_in >= edist % data(lc+NE)) then
        i_E = NE
      else
        i_E = binary_search(edist % data(lc+1:lc+NE), NE, E_in)
      end if
      
      ! Get interpolation type
      if (NR == 0) then
        INTT_in = LINEAR_LINEAR
      else
        !!! If I come across this in chi, then I need to account for it here.
        !!! If not in ENDF-7, then I should be okay.
        write(*,*) 'LAW 4 Error, INTT /= LINEAR_LINEAR', edist % law
      end if
      
      ! determine location of outgoing energies, pdf, cdf for E(l)
      lc = int(edist % data(2 + 2*NR + NE + i_E))

      ! determine type of interpolation and number of discrete lines
      INTTp = int(edist % data(lc + 1))
      NP    = int(edist % data(lc + 2))
      if (INTTp > 10) then
        INTT = mod(INTTp,10)
        ND = (INTTp - INTT)/10
      else
        INTT = INTTp
        ND = 0
      end if
      if (ND > 0) then
        ! discrete lines present
        message = "Discrete lines in continuous tabular distributed not &
             &yet supported"
        call fatal_error()
      end if
      
      ! Loop through energy groups
      lc = lc + 3
      lEout_min = lc
      do g = 1, size(energy_bins) - 1
        ! Find the location in edist%data corresponding to the upper
        ! bound of this energy group.
        ! The lower bound is already known: lEout_min
        do i_E = lEout_min, NP + lc - 2
          if (edist % data(i_E + 1) > energy_bins(g + 1)) exit
        end do
        if (i_E == NP + lc - 1) i_E = i_E - 1
        ! Since we have the CDF, the integral of chi(g) is simply:
        ! cdf_chi(E_high) - cdf_chi(E_low)
        ! Calculate cdf_chi(E_high),set equal to chi(g)
        if (INTT == LINEAR_LINEAR) then
          interp = (energy_bins(g + 1) - edist % data(i_E)) / &
            (edist % data(i_E + 1) - edist % data(i_E))
        elseif (INTT == HISTOGRAM) then
          interp = ZERO
        end if
        chi(g) = (edist % data(i_E + 2 * NP) + interp * &
          (edist % data(i_E + 1 + 2 * NP) - edist % data(i_E + 2 * NP)))
        ! Calculate cdf_chi(E_low), subtract from cdf_chi(E_high), which is chi
        ! The term to subtract off is the chi(g2<g) values.
        do g2 = 1, g - 1
          chi(g) = chi(g) - chi(g2)
        end do
        lEout_min = i_E
      end do
      
    case (5)
      ! =======================================================================
      ! GENERAL EVAPORATION SPECTRUM
      message = "Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported."
      call warning()
    case (7)
      ! =======================================================================
      ! MAXWELL FISSION SPECTRUM
      ! read number of interpolation regions and incoming energies 
      NR  = int(edist % data(1))
      NE  = int(edist % data(2 + 2*NR))
      
      ! determine nuclear temperature from tabulated function
      T = interpolate_tab1(edist % data, E_in)
      
      ! Get interpolation scheme
      lc = 2 + 2*NR
      if (E_in < edist % data(lc+1)) then
        i_E = 1
      elseif (E_in >= edist % data(lc+NE)) then
        i_E = NE
      else
        i_E = binary_search(edist % data(lc+1:lc+NE), NE, E_in)
      end if
      if (NR == 0) then
        INTT_in = LINEAR_LINEAR
      elseif (NR == 1) then
        INTT_in = int(edist % data(1 + NR + 1))
      elseif (NR > 1) then
        do g = 1, NR
          if (i_E < edist % data(1 + g)) then
            INTT_in = int(edist % data(1 + NR + g))
            exit
          end if
        end do
      end if
      
      ! determine restriction energy
      lc = 2 + 2*NR + 2*NE
      U = edist % data(lc + 1)
      x = (E_in - U) / T
      I = sqrt(T*T*T) * (sqrt(0.25_8*PI) * erf(x) - sqrt(x) * exp(-x))
      if (U < ZERO) U = E_in - U
      do g = 1, energy_groups
        ! The upper energy boundary
        Egp1 = energy_bins(g+1)
        if (Egp1 > U) Egp1 = U
        chi(g) = 0.5_8 * sqrt(PI * T * T * T) * erf( sqrt(Egp1 / T)) - &
          T * sqrt(Egp1) * exp(-Egp1 / T)
        ! The lower energy boundary
        Eg = energy_bins(g)
        if (Eg > U) Eg = U
        chi(g) = chi(g) - (0.25_8*sqrt(T*T*T*PI) * erf(sqrt(Eg / T)) - &
          T * sqrt(Eg) * exp(-Eg / T))
        chi(g) = chi(g) / I
      end do
      
    case (9)
      ! =======================================================================
      ! EVAPORATION SPECTRUM
      ! read number of interpolation regions and incoming energies 
      NR  = int(edist % data(1))
      NE  = int(edist % data(2 + 2*NR))
      
      ! determine nuclear temperature from tabulated function
      T = interpolate_tab1(edist % data, E_in)
      
      ! Get interpolation scheme
      lc = 2 + 2*NR
      if (E_in < edist % data(lc+1)) then
        i_E = 1
      elseif (E_in >= edist % data(lc+NE)) then
        i_E = NE
      else
        i_E = binary_search(edist % data(lc+1:lc+NE), NE, E_in)
      end if
      if (NR == 0) then
        INTT_in = LINEAR_LINEAR
      elseif (NR == 1) then
        INTT_in = int(edist % data(1 + NR + 1))
      elseif (NR > 1) then
        do g = 1, NR
          if (i_E < edist % data(1 + g)) then
            INTT_in = int(edist % data(1 + NR + g))
            exit
          end if
        end do
      end if
      
      ! determine restriction energy
      lc = 2 + 2*NR + 2*NE
      U = edist % data(lc + 1)
      x = (E_in - U) / T
      I = T * T * (ONE - exp(-x) * (ONE + x))
      if (U < ZERO) U = E_in - U
      do g = 1, energy_groups
        Egp1 = energy_bins(g+1)
        if (Egp1 > U) Egp1 = U
        chi(g) = T * exp(-Egp1 / T) * (T + Egp1)
        Eg = energy_bins(g)
        if (Eg > U) Eg = U
        chi(g) = chi(g) - (T * exp(-Eg / T) * (T + Eg))
        chi(g) = - chi(g) / I
      end do
    case (11)
      ! =======================================================================
      ! ENERGY-DEPENDENT WATT SPECTRUM
      ! read number of interpolation regions and incoming energies 
      NR  = int(edist % data(1))
      NE  = int(edist % data(2 + 2*NR))
      
      ! Get interpolation type
      if (NR == 0) then
        INTT_in = LINEAR_LINEAR
      else
        write(*,*) 'Error, INTT /= LINEAR_LINEAR', edist % law
      end if
      
      ! determine Watt parameter 'a' from tabulated function
      Watt_a = interpolate_tab1(edist % data, E_in)
      ! determine Watt parameter 'b' from tabulated function
      lc = 2 + 2*(NR + NE)
      Watt_b = interpolate_tab1(edist % data, E_in, lc + 1)
      ! read number of interpolation regions and incoming energies for
      ! parameter 'a'
      NR = int(edist % data(lc + 1))
      NE = int(edist % data(lc + 2 + 2*NR))
      ! determine restriction energy
      lc = lc + 2 + 2*(NR + NE)
      U = edist % data(lc + 1)
      x = (E_in - U) / Watt_a
      x0 = Watt_a * Watt_b * 0.25_8
      I = 0.25_8 * sqrt(PI * Watt_a ** 3 * Watt_b) * exp(x0) * &
        (erf(sqrt(x) - sqrt(x0)) + erf(sqrt(x) + sqrt(x0))) - &
        Watt_a * exp(-x * sinh(Watt_a * Watt_b * x))
      !Reuse Watt_b, and x
      Watt_b = sqrt(Watt_b)
      x = sqrt(PI * Watt_a) * Watt_b * exp(0.25_8 * Watt_a * Watt_b**2)
      if (U < ZERO) U = E_in - U
      do g = 1, energy_groups
        Egp1 = energy_bins(g+1)
        if (Egp1 > U) Egp1 = U
        chi(g) = &
          (-x * erf((Watt_a * Watt_b - TWO * sqrt(Egp1)/(TWO * Watt_a))) + &
          x * erf((Watt_a * Watt_b + TWO * sqrt(Egp1)/(TWO * Watt_a))) - &
          TWO * (exp(TWO * Watt_b * sqrt(Egp1)) * exp(-(Watt_a * Watt_b * sqrt(Egp1))/Watt_a)))
        Eg = energy_bins(g)
        if (Eg > U) Eg = U
        chi(g) = chi(g) - &
          (-x * erf((Watt_a * Watt_b - TWO * sqrt(Eg)/(TWO * Watt_a))) + &
          x * erf((Watt_a * Watt_b + TWO * sqrt(Eg)/(TWO * Watt_a))) - &
          TWO * (exp(TWO * Watt_b * sqrt(Eg)) * exp(-(Watt_a * Watt_b * sqrt(Eg))/Watt_a)))
        chi(g) = 0.25_8 * Watt_a *chi(g) / I
      end do
    case (12)
      ! =======================================================================
      ! Madland-Nix Fission Spectrum
      message = "Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported."
      call warning()
    case (44)
      ! =======================================================================
      ! KALBACH-MANN CORRELATED SCATTERING
      message = "Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported."
      call warning()
!     case (61)
!       ! =======================================================================
!       ! CORRELATED ENERGY AND ANGLE DISTRIBUTION
!       !!! Needed for 90232, 91231, 91233
!       message = "Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported."
!       call warning()
    case (66)
      ! =======================================================================
      ! N-BODY PHASE SPACE DISTRIBUTION
      message = "Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported."
      call warning()
    case (67)
      ! =======================================================================
      ! LABORATORY ENERGY-ANGLE LAW
      message = "Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported."
      call warning()
    end select
    
    ! Normalize chi(g) in case the interpolation rule causes it to be > 1.0
    I = ZERO
    do g = 1, energy_groups
      I = I + chi(g)
    end do
    if (I /= ONE) then
      I = ONE / I
      do g = 1, energy_groups
        chi(g) = chi(g) * I
      end do
    end if     
    
  end subroutine integrate_chi
  
  subroutine unionize_chi(chi_union, E_union, E_grid, chi_grid, NE_grid, &
    INT_grid, betas, probs, rxn_low, rxn_high)
    real(8), allocatable, intent(inout) :: chi_union(:,:)  ! Unioninzed Chi grid
    real(8), allocatable, intent(inout) :: E_union(:)      ! Unioninzed Energy grid
    real(8), intent(in)    :: E_grid(:,:)     ! The energy grid for each fission rxn
    real(8), intent(in)    :: chi_grid(:,:,:) ! Chis for the global grid E values
    integer, intent(in)    :: NE_grid(:)      ! NE values for each reaction
    integer, intent(in)    :: INT_grid(:,:)   ! The interpolation grid for each fission rxn
    real(8), intent(in)    :: betas(:,:)      ! Delayed-n Beta at each E_in
    real(8), intent(in)    :: probs(:,:)      ! Probabilities of each fission reaction for each E_in
    integer, intent(in)    :: rxn_low         ! Low bound of reaction index
    integer, intent(in)    :: rxn_high        ! High bound of reaction index
    
    real(8), allocatable :: chi_temp(:,:) ! Temporary array for the unionized chi
    real(8), allocatable :: E_temp(:)     ! Temporary array for the unionized E_grid
    integer, allocatable :: chi_loc(:)      ! Location of last-used energy
    logical, allocatable :: chi_done(:)     ! Flag to state whether a rxn's grid is complete
    logical, allocatable :: chi_hits(:)     ! Flag for when a chi value is the lowest
    real(8)              :: r               ! Interpolation Factor
    integer              :: i_energy, i_rxn ! Loop indices
    
    ! Chi_temp is sized to allow for all of this nuclide's reactions
    ! to not have the same energy points. One of the last steps here
    ! will be to allocate chi_union and copy in only the needed values
    allocate(chi_temp(size(chi_grid, dim = 1), SUM(NE_grid(rxn_low : rxn_high))))
    chi_temp = ZERO
    allocate(E_temp(size(chi_temp, dim = 2)))
    E_temp = ZERO
    allocate(chi_loc(rxn_low : rxn_high))
    chi_loc = 1
    allocate(chi_done(rxn_low : rxn_high))
    chi_done = .false.
    allocate(chi_hits(rxn_low : rxn_high))
    chi_hits = .false.
    
    do i_energy = 1, size(E_temp)
      ! Find the reactions which have the smallest values, stored in chi_hits
      call min_value_locs(chi_hits, chi_loc, chi_done, E_grid, rxn_low, rxn_high)
      ! Get the hits first
      do i_rxn = rxn_low, rxn_high
        if (.not. chi_done(i_rxn)) then
          if (chi_hits(i_rxn)) then
            ! Add the energy point to the union grid
            E_temp(i_energy) = E_grid(chi_loc(i_rxn), i_rxn)
            ! Sum this chi value to the union grid at the point it is defined at
            chi_temp(:, i_energy) = chi_temp(:, i_energy) + betas(chi_loc(i_rxn), i_rxn) * &
              probs(chi_loc(i_rxn), i_rxn) * chi_grid(:, chi_loc(i_rxn), i_rxn)
            chi_loc(i_rxn) = chi_loc(i_rxn) + 1
            if (chi_loc(i_rxn) > NE_grid(i_rxn)) chi_done(i_rxn) = .true.
          end if
        end if
      end do
      ! Now interpolate those which are not the next data point
      do i_rxn = rxn_low, rxn_high
        if (.not. chi_done(i_rxn)) then
          if (.not. chi_hits(i_rxn)) then
            ! If the i_energy point is outside the range of this reaction, then 
            ! do not score it.
            if (E_temp(i_energy) < E_grid(1, i_rxn)) cycle
            if (INT_grid(chi_loc(i_rxn), i_rxn) == HISTOGRAM) then
              ! If its a histogram, use the low energy value
              chi_temp(:, i_energy) = chi_temp(:, i_energy) + betas(chi_loc(i_rxn), i_rxn) * &
                probs(chi_loc(i_rxn), i_rxn) * chi_grid(:, chi_loc(i_rxn), i_rxn)
            elseif (INT_grid(chi_loc(i_rxn), i_rxn) == LINEAR_LINEAR) then
              ! For Linear-Linear, I need to find an interpolation factor first
              r = (E_temp(i_energy) - E_grid(chi_loc(i_rxn) - 1, i_rxn)) / &
                (E_grid(chi_loc(i_rxn), i_rxn) - E_grid(chi_loc(i_rxn) - 1, i_rxn))
              ! Now perform linear interpolation and sum to chi_temp
              !!! Check to see if this separate betas probs interpolation is worth anything
              chi_temp(:, i_energy) = chi_temp(:, i_energy) + &
                (betas(chi_loc(i_rxn) - 1, i_rxn) + r * &
                (betas(chi_loc(i_rxn), i_rxn) - betas(chi_loc(i_rxn) - 1, i_rxn)))* &
                (probs(chi_loc(i_rxn) - 1, i_rxn) + r * &
                (probs(chi_loc(i_rxn), i_rxn) - probs(chi_loc(i_rxn) - 1, i_rxn))) * &
                (chi_grid(:, chi_loc(i_rxn) - 1, i_rxn) + &
                r * (chi_grid(:, chi_loc(i_rxn), i_rxn) - &
                chi_grid(:, chi_loc(i_rxn) - 1, i_rxn)))
            else
              message = "Invalid Interpolation Grid Type: " // &
                to_str(INT_GRID(chi_loc(i_rxn), i_rxn))
              call warning()
              ! But still just do linear-linear.
              ! For Linear-Linear, I need to find an interpolation factor first
              r = (E_temp(i_energy) - E_grid(chi_loc(i_rxn) - 1, i_rxn)) / &
                (E_grid(chi_loc(i_rxn), i_rxn) - E_grid(chi_loc(i_rxn) - 1, i_rxn))
              ! Now perform linear interpolation and sum to chi_temp
              !!! Check to see if this separate betas probs interpolation is worth anything
              chi_temp(:, i_energy) = chi_temp(:, i_energy) + &
                (betas(chi_loc(i_rxn) - 1, i_rxn) + r * &
                (betas(chi_loc(i_rxn), i_rxn) - betas(chi_loc(i_rxn) - 1, i_rxn)))* &
                (probs(chi_loc(i_rxn) - 1, i_rxn) + r * &
                (probs(chi_loc(i_rxn), i_rxn) - probs(chi_loc(i_rxn) - 1, i_rxn))) * &
                (chi_grid(:, chi_loc(i_rxn) - 1, i_rxn) + &
                r * (chi_grid(:, chi_loc(i_rxn), i_rxn) - &
                chi_grid(:, chi_loc(i_rxn) - 1, i_rxn)))
            end if
          end if
        end if
      end do
      ! Now that all data points are combined, re-check the normalization
      ! This is necessary since the interpolation likely has caused chi
      ! to become "un-normalized" again.
      r = sum(chi_temp(:, i_energy))
      if (r > 0.0) then
        chi_temp(:, i_energy) = chi_temp(:, i_energy) / r
      end if
      
      ! Check to see if all chi values have reached their last data point
      ! If so, then we can complete the unionizing
      if (all(chi_done)) exit
    end do
    
    ! Allocate the final union grid spaces
    allocate(chi_union(size(chi_temp, dim = 1), i_energy))
    allocate(E_union(i_energy))
    
    chi_union(:, :) = chi_temp(:, 1 : i_energy)
    E_union(:) = E_temp(1 : i_energy)
    
    ! Free memory
    deallocate(chi_temp)
    deallocate(E_temp)
    deallocate(chi_loc)
    deallocate(chi_done)
    deallocate(chi_hits)
  end subroutine unionize_chi
  
  subroutine thin_union_chi(chi, E, tol, max_err)
    real(8), allocatable, intent(inout) :: chi(:,:)  ! Unioninzed Chi grid
    real(8), allocatable, intent(inout) :: E(:)      ! Unioninzed Energy grid
    real(8), intent(in)    :: tol       ! The fractional thinning tolerance
    real(8), intent(out)   :: max_err   ! The maximum error of thinning
    
    real(8) :: err
    integer :: i_low, i_high ! low and high bounds of the grid
    integer :: i_test        ! value to test
    integer :: g             ! energy group index
    integer :: to_strip(size(E))
    integer :: num_strip
    logical :: cont_loop, all_pass
    real(8) :: interp
    real(8) :: test_val(size(chi,dim = 1))
    real(8), allocatable :: thin(:), E_thin(:)
    
    to_strip = 0
    num_strip = 0
    cont_loop = .true.
    i_low = 1
    max_err = ZERO
    i_test = i_low
    do while (cont_loop)
      i_high = i_low + 2
      i_test = i_test + 1
      interp = (E(i_test) - E(i_low)) / (E(i_test) - E(i_high))
      do g = 1, size(chi, dim = 1)
        test_val(g) = chi(g, i_low) + interp * (chi(g, i_high) - chi(g, i_low))
      end do
      all_pass = .true.
      do g = 1, size(chi, dim = 1)
        if (test_val(g) > tol) then
          all_pass = .false.
        end if
      end do
      write(*,*) 'test_val =',test_val
      if (all_pass) then
        i_high = i_high + 1
        num_strip = num_strip + 1
        to_strip(num_strip) = i_test
        i_test = i_test + 1
      else
        i_low = i_high - 1
        i_test = i_low
      end if
      if (i_low > (size(E) - 2)) cont_loop = .false.
    end do
    write(*,*) 'to_strip: ',to_strip
    ! Reconstruct the grid, but using the thinned grid
    allocate(thin(num_strip))
    allocate(E_thin(num_strip))
    i_low = 0
    i_high = 0
    
  end subroutine thin_union_chi
  
  subroutine min_value_locs(hits, chi_loc, chi_done, E_grid, low, high)
    integer, intent(in)    :: low         ! Low index of chi to check
    integer, intent(in)    :: high        ! High index of chi to check
    logical, intent(inout) :: hits(low:high)     ! The reactions which are the lowest E
    integer, intent(in)    :: chi_loc(low:high)  ! Location of last-used energy
    logical, intent(in)    :: chi_done(low:high) ! Flag to state whether a rxn's grid is complete
    real(8), intent(in)    :: E_grid(:,:) ! The energy grid for each fission rxn
    
    integer :: i
    real(8) :: min_val
    
    hits = .false.
    min_val = INFINITY
    
    ! Find the next smallest energy value
    do i = low, high
      if (.not. chi_done(i)) then
        if (E_grid(chi_loc(i), i) < min_val) &
          min_val = E_grid(chi_loc(i), i)
      end if
    end do
    
    ! Find all locs  that are within a certain tolerance of the min value
    do i = low, high
      if (.not. chi_done(i)) then
        if (abs(E_grid(chi_loc(i), i) - min_val) < FP_PRECISION) then
          hits(i) = .true.
        end if
      end if
    end do
  
  end subroutine min_value_locs
  
!===============================================================================
! PRINT_CHI prints the chi data to the specified output file
! in the specified format.
!===============================================================================  
  
  subroutine print_chi(lib_format, chi_t, chi_p, chi_d, E_t, E_p, E_d)
    integer,              intent(in) :: lib_format ! Library output type
    real(8), allocatable, intent(in) :: chi_t(:,:) ! Unionized Total Chi
    real(8), allocatable, intent(in) :: chi_p(:,:) ! Unionized Prompt Chi
    real(8), allocatable, intent(in) :: chi_d(:,:) ! Unionized Delayed Chi
    real(8), allocatable, intent(in) :: E_t(:)     ! Unionized Total Energy
    real(8), allocatable, intent(in) :: E_p(:)     ! Unionized Prompt Energy
    real(8), allocatable, intent(in) :: E_d(:)     ! Unionized Delayed Energy
    
    if (lib_format == ASCII) then
      call print_chi_ascii(chi_t, chi_p, chi_d, E_t, E_p, E_d)
    else if (lib_format == BINARY) then
      call print_chi_bin(chi_t, chi_p, chi_d, E_t, E_p, E_d)
    else if (lib_format == HUMAN) then
      call print_chi_human(chi_t, chi_p, chi_d, E_t, E_p, E_d)
    else if (lib_format == HDF5) then
      ! TBI
    end if
    
  end subroutine print_chi

!===============================================================================
! PRINT_CHI_ASCII prints the chi data to the specified output file
! in an ASCII format.
!===============================================================================
  
  subroutine print_chi_ascii(chi_t, chi_p, chi_d, E_t, E_p, E_d)
    real(8), allocatable, intent(in) :: chi_t(:,:) ! Unionized Total Chi
    real(8), allocatable, intent(in) :: chi_p(:,:) ! Unionized Prompt Chi
    real(8), allocatable, intent(in) :: chi_d(:,:) ! Unionized Delayed Chi
    real(8), allocatable, intent(in) :: E_t(:)     ! Unionized Total Energy
    real(8), allocatable, intent(in) :: E_p(:)     ! Unionized Prompt Energy
    real(8), allocatable, intent(in) :: E_d(:)     ! Unionized Delayed Energy
    
    character(MAX_LINE_LEN) :: line
  
    ! Assumes that the file and header information is already printed 
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format with at max 4 entries per line: 
    ! <n_E_t>
    ! <1:n_E_t energies>
    ! <1:n_E_t chi_t(g,E)
    ! <n_E_p>
    ! <1:n_E_p energies>
    ! <1:n_E_p chi_p(g,E)
    ! <n_E_d>
    ! <1:n_E_d energies>
    ! <1:n_E_d chi_d(g,E)
    
    ! Begin writing:
    
    ! <n_E_t>
    line = ''
    write(line,'(I20)') size(E_t)
    write(UNIT_NUC,'(A)') trim(line)
    
    ! <1:n_E_t energies>
    call print_ascii_array(E_t, UNIT_NUC)
    
    ! <1:n_E_t chi_t>
    call print_ascii_array(pack(chi_t, .true.), UNIT_NUC)
    
    ! <n_E_p>
    line = ''
    write(line,'(I20)') size(E_p)
    write(UNIT_NUC,'(A)') trim(line)
    
    ! <1:n_E_p energies>
    call print_ascii_array(E_p, UNIT_NUC)
    
    ! <1:n_E_p chi_p>
    call print_ascii_array(pack(chi_p, .true.), UNIT_NUC)
    
    ! <n_E_d>
    line = ''
    write(line,'(I20)') size(E_d)
    write(UNIT_NUC,'(A)') trim(line)
    
    ! <1:n_E_d energies>
    call print_ascii_array(E_d, UNIT_NUC)
    
      ! <1:n_E_d chi_d>
    call print_ascii_array(pack(chi_d, .true.), UNIT_NUC)
    
  end subroutine print_chi_ascii
  
!===============================================================================
! PRINT_CHI_HUMAN prints the chi data to the specified output file
! in a human-readable ASCII format.
!===============================================================================
  
  subroutine print_chi_human(chi_t, chi_p, chi_d, E_t, E_p, E_d)
    real(8), allocatable, intent(in) :: chi_t(:,:) ! Unionized Total Chi
    real(8), allocatable, intent(in) :: chi_p(:,:) ! Unionized Prompt Chi
    real(8), allocatable, intent(in) :: chi_d(:,:) ! Unionized Delayed Chi
    real(8), allocatable, intent(in) :: E_t(:)     ! Unionized Total Energy
    real(8), allocatable, intent(in) :: E_p(:)     ! Unionized Prompt Energy
    real(8), allocatable, intent(in) :: E_d(:)     ! Unionized Delayed Energy
    
    integer :: iE
    
    character(MAX_LINE_LEN) :: line
  
    ! Assumes that the file and header information is already printed 
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format with at max 4 entries per line: 
    ! <n_E_t>
    ! <1:n_E_t energies>
    ! <1:n_E_t E_t, chi_t(g,E)
    ! <n_E_p>
    ! <1:n_E_p energies>
    ! <1:n_E_p E_p, chi_p(g,E)
    ! <n_E_d>
    ! <1:n_E_d energies>
    ! <1:n_E_d E_d, chi_d(g,E)
    
    ! Begin writing:
    
    ! <n_E_t>
    line = ''
    write(line,'(I20)') size(E_t)
    write(UNIT_NUC,'(A)') trim(line)
    
    ! <1:n_E_t energies>
    call print_ascii_array(E_t, UNIT_NUC)
    
    ! <1:n_E_t E_t, chi_t>
    do iE = 1, size(E_t)
      write(line, '(1PE20.12)') E_t(iE)
      write(UNIT_NUC, '(A)') trim(line)
      call print_ascii_array(chi_t(:, iE), UNIT_NUC)
    end do
    
    ! <n_E_p>
    line = ''
    write(line,'(I20)') size(E_p)
    write(UNIT_NUC,'(A)') trim(line)
    
    ! <1:n_E_p energies>
    call print_ascii_array(E_p, UNIT_NUC)
    
    ! <1:n_E_p E_p, chi_p>
    do iE = 1, size(E_p)
      write(line, '(1PE20.12)') E_p(iE)
      write(UNIT_NUC, '(A)') trim(line)
      call print_ascii_array(chi_p(:, iE), UNIT_NUC)
    end do
    
    ! <n_E_d>
    line = ''
    write(line,'(I20)') size(E_d)
    write(UNIT_NUC,'(A)') trim(line)
    
    ! <1:n_E_d energies>
    call print_ascii_array(E_d, UNIT_NUC)
    
    ! <1:n_E_d E_d, chi_p>
    do iE = 1, size(E_d)
      write(line, '(1PE20.12)') E_d(iE)
      write(UNIT_NUC, '(A)') trim(line)
      call print_ascii_array(chi_d(:, iE), UNIT_NUC)
    end do
    
  end subroutine print_chi_human

!===============================================================================
! PRINT_CHI_BIN prints the chi data to the specified output file
! in Fortran binary (stream) format.
!===============================================================================
  
  subroutine print_chi_bin(chi_t, chi_p, chi_d, E_t, E_p, E_d)
    real(8), allocatable, intent(in) :: chi_t(:,:) ! Unionized Total Chi
    real(8), allocatable, intent(in) :: chi_p(:,:) ! Unionized Prompt Chi
    real(8), allocatable, intent(in) :: chi_d(:,:) ! Unionized Delayed Chi
    real(8), allocatable, intent(in) :: E_t(:)     ! Unionized Total Energy
    real(8), allocatable, intent(in) :: E_p(:)     ! Unionized Prompt Energy
    real(8), allocatable, intent(in) :: E_d(:)     ! Unionized Delayed Energy
  
    ! Assumes that the file and header information is already printed 
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format:
    ! <n_E_t>
    ! <1:n_E_t energies>
    ! <1:n_E_t chi_t(g,E)
    ! <n_E_p>
    ! <1:n_E_p energies>
    ! <1:n_E_p chi_p(g,E)
    ! <n_E_d>
    ! <1:n_E_d energies>
    ! <1:n_E_d chi_d(g,E)
    
    ! Begin writing:
    
    ! <n_E_t>
    write(UNIT_NUC) size(E_t)
    
    ! <1:n_E_t energies>
    write(UNIT_NUC) E_t
    
    ! <1:n_E_t chi_t>
    write(UNIT_NUC) chi_t
    
    ! <n_E_p>
    write(UNIT_NUC) size(E_p)
    
    ! <1:n_E_p energies>
    write(UNIT_NUC) E_p
    
    ! <1:n_E_p chi_p>
    write(UNIT_NUC) chi_p
    
    ! <n_E_d>
    write(UNIT_NUC) size(E_d)
    
    ! <1:n_E_d energies>
    write(UNIT_NUC) E_d
    
    ! <1:n_E_d chi_d>
    write(UNIT_NUC) chi_d
    
  end subroutine print_chi_bin
  
end module ndpp_chi
